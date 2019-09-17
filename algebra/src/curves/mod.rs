use crate::{
    bytes::{FromBytes, ToBytes},
    fields::{Field, PrimeField, SquareRootField},
    groups::Group,
};
use crate::UniformRand;
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};

pub mod bls12_377;
pub mod bls12_381;
pub mod edwards_bls12;
pub mod edwards_sw6;
pub mod jubjub;
pub mod mnt6;
pub mod models;
pub mod sw6;

#[cfg(test)]
pub mod tests;

pub use self::models::*;

pub trait PairingEngine: Sized + 'static + Copy + Debug + Sync + Send {
    /// This is the scalar field of the G1/G2 groups.
    type Fr: PrimeField + SquareRootField + Into<<Self::Fr as PrimeField>::BigInt>;

    /// The projective representation of an element in G1.
    type G1Projective: ProjectiveCurve<
            BaseField = Self::Fq,
            ScalarField = Self::Fr,
            Affine = Self::G1Affine,
        > + From<Self::G1Affine>;

    /// The affine representation of an element in G1.
    type G1Affine: AffineCurve<
            BaseField = Self::Fq,
            ScalarField = Self::Fr,
            Projective = Self::G1Projective,
        > + PairingCurve<PairWith = Self::G2Affine, PairingResult = Self::Fqk>
        + From<Self::G1Projective>;

    /// The projective representation of an element in G2.
    type G2Projective: ProjectiveCurve<
            BaseField = Self::Fqe,
            ScalarField = Self::Fr,
            Affine = Self::G2Affine,
        > + From<Self::G2Affine>;

    /// The affine representation of an element in G2.
    type G2Affine: AffineCurve<
            BaseField = Self::Fqe,
            ScalarField = Self::Fr,
            Projective = Self::G2Projective,
        > + PairingCurve<PairWith = Self::G1Affine, PairingResult = Self::Fqk>
        + From<Self::G2Projective>;

    /// The base field that hosts G1.
    type Fq: PrimeField + SquareRootField;

    /// The extension field that hosts G2.
    type Fqe: SquareRootField;

    /// The extension field that hosts the target group of the pairing.
    type Fqk: Field;

    /// Perform a miller loop with some number of (G1, G2) pairs.
    #[must_use]
    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<
            Item = &'a (
                &'a <Self::G1Affine as PairingCurve>::Prepared,
                &'a <Self::G2Affine as PairingCurve>::Prepared,
            ),
        >;

    /// Perform final exponentiation of the result of a miller loop.
    #[must_use]
    fn final_exponentiation(_: &Self::Fqk) -> Option<Self::Fqk>;

    /// Computes a product of pairings.
    #[must_use]
    fn product_of_pairings<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<
            Item = &'a (
                &'a <Self::G1Affine as PairingCurve>::Prepared,
                &'a <Self::G2Affine as PairingCurve>::Prepared,
            ),
        >,
    {
        Self::final_exponentiation(&Self::miller_loop(i)).unwrap()
    }

    /// Performs multiple pairing operations
    #[must_use]
    fn pairing<G1, G2>(p: G1, q: G2) -> Self::Fqk
    where
        G1: Into<Self::G1Affine>,
        G2: Into<Self::G2Affine>,
    {
        Self::final_exponentiation(&Self::miller_loop(
            [(&(p.into().prepare()), &(q.into().prepare()))].into_iter(),
        ))
        .unwrap()
    }
}

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait ProjectiveCurve:
    Eq
    + Sized
    + ToBytes
    + FromBytes
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

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

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

pub trait PairingCurve: AffineCurve {
    type Engine: PairingEngine<Fr = Self::ScalarField>;
    type Prepared: ToBytes + Default + Clone + Send + Sync + Debug + 'static;
    type PairWith: PairingCurve<PairWith = Self>;
    type PairingResult: Field;

    /// Prepares this element for pairing purposes.
    #[must_use]
    fn prepare(&self) -> Self::Prepared;

    /// Perform a pairing
    #[must_use]
    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult;
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
