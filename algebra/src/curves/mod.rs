use crate::UniformRand;
use crate::{
    bits::{FromCompressedBits, ToCompressedBits},
    bytes::{FromBytes, ToBytes},
    fields::{Field, PrimeField, SquareRootField},
    groups::Group,
    CanonicalDeserialize, CanonicalSerialize, Error, FromBytesChecked, SemanticallyValid,
};
use serde::{Deserialize, Serialize};
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};

pub mod models;

#[cfg(feature = "bls12_377")]
pub mod bls12_377;

#[cfg(feature = "bls12_381")]
pub mod bls12_381;

#[cfg(feature = "bn_382")]
pub mod bn_382;

#[cfg(feature = "edwards_bls12")]
pub mod edwards_bls12;

#[cfg(feature = "edwards_sw6")]
pub mod edwards_sw6;

#[cfg(feature = "jubjub")]
pub mod jubjub;

#[cfg(feature = "mnt4_753")]
pub mod mnt4753;

#[cfg(feature = "mnt6_753")]
pub mod mnt6753;

#[cfg(feature = "mnt6")]
pub mod mnt6;

#[cfg(feature = "sw6")]
pub mod sw6;

#[cfg(feature = "tweedle")]
pub mod tweedle;

#[cfg(test)]
pub mod tests;

pub use self::models::*;

pub trait PairingEngine: Sized + 'static + Copy + Debug + Sync + Send + Eq + PartialEq {
    /// This is the scalar field of the G1/G2 groups.
    type Fr: PrimeField + SquareRootField + Into<<Self::Fr as PrimeField>::BigInt>;

    /// The projective representation of an element in G1.
    type G1Projective: ProjectiveCurve<BaseField = Self::Fq, ScalarField = Self::Fr, Affine = Self::G1Affine>
        + From<Self::G1Affine>
        + Into<Self::G1Affine>;

    /// The affine representation of an element in G1.
    type G1Affine: AffineCurve<BaseField = Self::Fq, ScalarField = Self::Fr, Projective = Self::G1Projective>
        + From<Self::G1Projective>
        + Into<Self::G1Projective>
        + Into<Self::G1Prepared>;

    /// A G1 element that has been preprocessed for use in a pairing.
    type G1Prepared: ToBytes
        + FromBytes
        + Serialize
        + for<'a> Deserialize<'a>
        + Default
        + Clone
        + Send
        + Sync
        + Debug
        + From<Self::G1Affine>;

    /// The projective representation of an element in G2.
    type G2Projective: ProjectiveCurve<BaseField = Self::Fqe, ScalarField = Self::Fr, Affine = Self::G2Affine>
        + From<Self::G2Affine>
        + Into<Self::G2Affine>;

    /// The affine representation of an element in G2.
    type G2Affine: AffineCurve<BaseField = Self::Fqe, ScalarField = Self::Fr, Projective = Self::G2Projective>
        + From<Self::G2Projective>
        + Into<Self::G2Projective>
        + Into<Self::G2Prepared>;

    /// A G2 element that has been preprocessed for use in a pairing.
    type G2Prepared: ToBytes
        + FromBytes
        + Serialize
        + for<'a> Deserialize<'a>
        + Default
        + Eq
        + PartialEq
        + Clone
        + Send
        + Sync
        + Debug
        + From<Self::G2Affine>;

    /// The base field that hosts G1.
    type Fq: PrimeField + SquareRootField;

    /// The extension field that hosts G2.
    type Fqe: SquareRootField;

    /// The extension field that hosts the target group of the pairing.
    type Fqk: Field;

    /// Perform a miller loop with some number of (G1, G2) pairs.
    #[must_use]
    fn miller_loop<'a, I>(i: I) -> Result<Self::Fqk, Error>
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>;

    /// Perform final exponentiation of the result of a miller loop.
    #[must_use]
    fn final_exponentiation(_: &Self::Fqk) -> Result<Self::Fqk, Error>;

    /// Computes a product of pairings.
    #[must_use]
    fn product_of_pairings<'a, I>(i: I) -> Result<Self::Fqk, Error>
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        Self::final_exponentiation(&Self::miller_loop(i)?)
    }

    /// Performs multiple pairing operations
    #[must_use]
    fn pairing<G1, G2>(p: G1, q: G2) -> Result<Self::Fqk, Error>
    where
        G1: Into<Self::G1Affine>,
        G2: Into<Self::G2Affine>,
    {
        let g1_prep = Self::G1Prepared::from(p.into());
        let g2_prep = Self::G2Prepared::from(q.into());
        Self::product_of_pairings(std::iter::once(&(g1_prep, g2_prep)))
    }
}

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait ProjectiveCurve:
    Eq
    + Sized
    + ToBytes
    + FromBytes
    + Serialize
    + for<'a> Deserialize<'a>
    + CanonicalSerialize
    + CanonicalDeserialize
    + SemanticallyValid
    + FromBytesChecked
    + Copy
    + Clone
    + Default
    + Send
    + Sync
    + Hash
    + Debug
    + Display
    + UniformRand
    + 'static
    + Neg<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
{
    type ScalarField: PrimeField + SquareRootField + Into<<Self::ScalarField as PrimeField>::BigInt>;
    type BaseField: Field;
    type Affine: AffineCurve<Projective = Self, ScalarField = Self::ScalarField>;

    /// Returns the additive identity.
    #[must_use]
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    #[must_use]
    fn prime_subgroup_generator() -> Self;

    /// Determines if this point is the point at infinity.
    #[must_use]
    fn is_zero(&self) -> bool;

    /// Checks that the current point is on curve and is in the
    /// prime order subgroup
    #[must_use]
    fn group_membership_test(&self) -> bool;

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

    fn batch_normalization_into_affine(mut v: Vec<Self>) -> Vec<Self::Affine> {
        Self::batch_normalization(v.as_mut_slice());
        v.into_iter().map(|p| p.into_affine()).collect()
    }

    /// Checks if the point is already "normalized" so that
    /// cheap affine conversion is possible.
    #[must_use]
    fn is_normalized(&self) -> bool;

    /// Doubles this element.
    #[must_use]
    fn double(&self) -> Self {
        let mut copy = *self;
        copy.double_in_place();
        copy
    }

    fn double_in_place(&mut self) -> &mut Self;

    /// Adds an affine element to this element.
    fn add_assign_mixed(&mut self, other: &Self::Affine);

    /// Performs scalar multiplication of this element.
    fn mul_assign<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&mut self, other: S);

    /// Converts this element into its affine representation.
    #[must_use]
    fn into_affine(&self) -> Self::Affine;

    /// Recommends a wNAF window table size given a scalar. Always returns a
    /// number between 2 and 22, inclusive.
    #[must_use]
    fn recommended_wnaf_for_scalar(scalar: <Self::ScalarField as PrimeField>::BigInt) -> usize;

    /// Recommends a wNAF window size given the number of scalars you intend to
    /// multiply a base by. Always returns a number between 2 and 22,
    /// inclusive.
    #[must_use]
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize;
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait AffineCurve:
    Eq
    + Sized
    + ToBytes
    + FromBytes
    + Serialize
    + for<'a> Deserialize<'a>
    + CanonicalSerialize
    + CanonicalDeserialize
    + SemanticallyValid
    + FromBytesChecked
    + ToCompressedBits
    + FromCompressedBits
    + Copy
    + Clone
    + Default
    + Send
    + Sync
    + Hash
    + Debug
    + Display
    + Neg<Output = Self>
    + 'static
{
    type ScalarField: PrimeField + SquareRootField + Into<<Self::ScalarField as PrimeField>::BigInt>;
    type BaseField: Field;
    type Projective: ProjectiveCurve<Affine = Self, ScalarField = Self::ScalarField>;

    /// Returns the additive identity.
    #[must_use]
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    #[must_use]
    fn prime_subgroup_generator() -> Self;

    /// Determines if this point represents the point at infinity; the
    /// additive identity.
    #[must_use]
    fn is_zero(&self) -> bool;

    /// Returns a group element if the set of bytes forms a valid group element,
    /// otherwise returns None. This function is primarily intended for sampling
    /// random group elements from a hash-function or RNG output.
    fn from_random_bytes(bytes: &[u8]) -> Option<Self>;

    /// Checks that the current point is on curve and is in the
    /// prime order subgroup
    #[must_use]
    fn group_membership_test(&self) -> bool;

    /// Adds, for each vector in 'to_add', its elements together
    /// using Affine point arithmetic
    fn add_points(to_add: &mut [Vec<Self>]);

    /// Performs scalar multiplication of this element with mixed addition.
    #[must_use]
    fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&self, other: S)
        -> Self::Projective;

    /// Converts this element into its projective representation.
    #[must_use]
    fn into_projective(&self) -> Self::Projective;

    /// Multiply this element by the cofactor.
    #[must_use]
    fn mul_by_cofactor(&self) -> Self;

    /// Multiply this element by the inverse of the cofactor modulo the size of
    /// `Self::ScalarField`.
    #[must_use]
    fn mul_by_cofactor_inv(&self) -> Self;
}

impl<C: ProjectiveCurve> Group for C {
    type ScalarField = C::ScalarField;
    #[must_use]
    fn zero() -> Self {
        <C as ProjectiveCurve>::zero()
    }

    #[must_use]
    fn is_zero(&self) -> bool {
        <C as ProjectiveCurve>::is_zero(&self)
    }

    #[inline]
    #[must_use]
    fn double(&self) -> Self {
        let mut tmp = *self;
        tmp += self;
        tmp
    }

    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        <C as ProjectiveCurve>::double_in_place(self)
    }
}

/// Preprocess a G1 element for use in a pairing.
pub fn prepare_g1<E: PairingEngine>(g: impl Into<E::G1Affine>) -> E::G1Prepared {
    let g: E::G1Affine = g.into();
    E::G1Prepared::from(g)
}

/// Preprocess a G2 element for use in a pairing.
pub fn prepare_g2<E: PairingEngine>(g: impl Into<E::G2Affine>) -> E::G2Prepared {
    let g: E::G2Affine = g.into();
    E::G2Prepared::from(g)
}
