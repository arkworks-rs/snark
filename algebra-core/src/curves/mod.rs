use crate::{
    biginteger::BigInteger,
    bytes::{FromBytes, ToBytes},
    fields::{Field, PrimeField, SquareRootField},
    groups::Group,
    CanonicalDeserialize, CanonicalSerialize, ConstantSerializedSize, UniformRand, Vec,
};
use core::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, MulAssign, Neg, Sub, SubAssign},
};
use num_traits::Zero;

pub mod batch_verify;
pub use self::batch_verify::*;

pub mod models;

pub use self::models::*;

pub trait PairingEngine: Sized + 'static + Copy + Debug + Sync + Send {
    /// This is the scalar field of the G1/G2 groups.
    type Fr: PrimeField + SquareRootField;

    /// The projective representation of an element in G1.
    type G1Projective: ProjectiveCurve<BaseField = Self::Fq, ScalarField = Self::Fr, Affine = Self::G1Affine>
        + From<Self::G1Affine>
        + Into<Self::G1Affine>
        + MulAssign<Self::Fr>; // needed due to https://github.com/rust-lang/rust/issues/69640

    /// The affine representation of an element in G1.
    type G1Affine: AffineCurve<BaseField = Self::Fq, ScalarField = Self::Fr, Projective = Self::G1Projective>
        + From<Self::G1Projective>
        + Into<Self::G1Projective>
        + Into<Self::G1Prepared>;

    /// A G1 element that has been preprocessed for use in a pairing.
    type G1Prepared: ToBytes + Default + Clone + Send + Sync + Debug + From<Self::G1Affine>;

    /// The projective representation of an element in G2.
    type G2Projective: ProjectiveCurve<BaseField = Self::Fqe, ScalarField = Self::Fr, Affine = Self::G2Affine>
        + From<Self::G2Affine>
        + Into<Self::G2Affine>
        + MulAssign<Self::Fr>; // needed due to https://github.com/rust-lang/rust/issues/69640

    /// The affine representation of an element in G2.
    type G2Affine: AffineCurve<BaseField = Self::Fqe, ScalarField = Self::Fr, Projective = Self::G2Projective>
        + From<Self::G2Projective>
        + Into<Self::G2Projective>
        + Into<Self::G2Prepared>;

    /// A G2 element that has been preprocessed for use in a pairing.
    type G2Prepared: ToBytes + Default + Clone + Send + Sync + Debug + From<Self::G2Affine>;

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
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>;

    /// Perform final exponentiation of the result of a miller loop.
    #[must_use]
    fn final_exponentiation(_: &Self::Fqk) -> Option<Self::Fqk>;

    /// Computes a product of pairings.
    #[must_use]
    fn product_of_pairings<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
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
        let g1_prep = Self::G1Prepared::from(p.into());
        let g2_prep = Self::G2Prepared::from(q.into());
        Self::product_of_pairings(core::iter::once(&(g1_prep, g2_prep)))
    }
}

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait ProjectiveCurve:
    Eq
    + 'static
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
    + Zero
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<<Self as ProjectiveCurve>::ScalarField>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + core::iter::Sum<Self>
    + for<'a> core::iter::Sum<&'a Self>
    + From<<Self as ProjectiveCurve>::Affine>
{
    const COFACTOR: &'static [u64];
    type ScalarField: PrimeField + SquareRootField;
    type BaseField: Field;
    type Affine: AffineCurve<Projective = Self, ScalarField = Self::ScalarField, BaseField = Self::BaseField>
        + From<Self>
        + Into<Self>;

    /// Returns a fixed generator of unknown exponent.
    #[must_use]
    fn prime_subgroup_generator() -> Self;

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

    /// Normalizes a slice of projective elements and outputs a vector
    /// containing the affine equivalents.
    fn batch_normalization_into_affine(v: &[Self]) -> Vec<Self::Affine> {
        let mut v = v.to_vec();
        Self::batch_normalization(&mut v);
        v.into_iter().map(|v| v.into()).collect()
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

    /// Doubles this element in place.
    fn double_in_place(&mut self) -> &mut Self;

    /// Converts self into the affine representation.
    fn into_affine(&self) -> Self::Affine {
        (*self).into()
    }

    /// Set `self` to be `self + other`, where `other: Self::Affine`.
    /// This is usually faster than adding `other` in projective form.
    fn add_mixed(mut self, other: &Self::Affine) -> Self {
        self.add_assign_mixed(other);
        self
    }

    /// Set `self` to be `self + other`, where `other: Self::Affine`.
    /// This is usually faster than adding `other` in projective form.
    fn add_assign_mixed(&mut self, other: &Self::Affine);

    /// Performs scalar multiplication of this element.
    fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(mut self, other: S) -> Self {
        let mut res = Self::zero();

        let mut found_one = false;

        for i in crate::fields::BitIterator::new(other.into()) {
            if found_one {
                res.double_in_place();
            } else {
                found_one = i;
            }

            if i {
                res += self;
            }
        }

        self = res;
        self
    }
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait AffineCurve:
    Eq
    + 'static
    + Sized
    + ToBytes
    + FromBytes
    + CanonicalSerialize
    + ConstantSerializedSize
    + CanonicalDeserialize
    + Copy
    + Clone
    + Default
    + Send
    + Sync
    + Hash
    + Debug
    + Display
    + Zero
    + Neg<Output = Self>
    + From<<Self as AffineCurve>::Projective>
    + BatchGroupArithmetic
{
    const COFACTOR: &'static [u64];
    type ScalarField: PrimeField + SquareRootField + Into<<Self::ScalarField as PrimeField>::BigInt>;
    type BaseField: Field;
    type Projective: ProjectiveCurve<Affine = Self, ScalarField = Self::ScalarField, BaseField = Self::BaseField>
        + From<Self>
        + Into<Self>
        + MulAssign<Self::ScalarField>; // needed due to https://github.com/rust-lang/rust/issues/69640

    /// Returns a fixed generator of unknown exponent.
    #[must_use]
    fn prime_subgroup_generator() -> Self;

    /// Converts self into the projective representation.
    fn into_projective(&self) -> Self::Projective {
        (*self).into()
    }

    /// Returns a group element if the set of bytes forms a valid group element,
    /// otherwise returns None. This function is primarily intended for sampling
    /// random group elements from a hash-function or RNG output.
    fn from_random_bytes(bytes: &[u8]) -> Option<Self>;

    /// Performs scalar multiplication of this element with mixed addition.
    #[must_use]
    fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&self, other: S)
        -> Self::Projective;

    /// Multiply this element by the cofactor and output the
    /// resulting projective element.
    #[must_use]
    fn mul_by_cofactor_to_projective(&self) -> Self::Projective;

    /// Multiply this element by the cofactor.
    #[must_use]
    fn mul_by_cofactor(&self) -> Self {
        self.mul_by_cofactor_to_projective().into()
    }

    /// Multiply this element by the inverse of the cofactor in
    /// `Self::ScalarField`.
    #[must_use]
    fn mul_by_cofactor_inv(&self) -> Self;
}

impl<C: ProjectiveCurve> Group for C {
    type ScalarField = C::ScalarField;

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

/// A cycle of pairing-friendly elliptic curves.
pub trait CycleEngine: Sized + 'static + Copy + Debug + Sync + Send
where
    <Self::E2 as PairingEngine>::G1Projective: MulAssign<<Self::E1 as PairingEngine>::Fq>,
    <Self::E2 as PairingEngine>::G2Projective: MulAssign<<Self::E1 as PairingEngine>::Fq>,
{
    type E1: PairingEngine;
    type E2: PairingEngine<
        Fr = <Self::E1 as PairingEngine>::Fq,
        Fq = <Self::E1 as PairingEngine>::Fr,
    >;
}

pub trait BatchGroupArithmetic
where
    Self: Sized + Clone + Copy + Zero + Neg<Output = Self>,
{
    // This function consumes the scalars
    // We can make this more generic in the future to use other than u16.

    // TODO: Generalise to A != 0
    // Computes [-p, p, -3p, 3p, ..., -2^wp, 2^wp]
    fn batch_wnaf_tables(bases: &[Self], w: usize) -> Vec<Vec<Self>> {
        let half_size = 1 << w;
        let batch_size = bases.len();

        let mut tables = vec![Vec::<Self>::with_capacity(half_size); batch_size];

        let mut a_2 = bases[..].to_vec();
        let mut tmp = bases[..].to_vec();

        let instr = (0..batch_size).collect::<Vec<usize>>();
        Self::batch_double_in_place(&mut a_2, &instr[..]);

        for i in 0..half_size {
            if i != 0 {
                let instr = (0..batch_size)
                    .map(|x| (x, x))
                    .collect::<Vec<(usize, usize)>>();
                Self::batch_add_in_place(&mut tmp, &mut a_2.to_vec()[..], &instr[..]);
            }

            for (table, p) in tables.iter_mut().zip(&tmp) {
                table.push(p.clone());
            }
        }
        tables
    }

    // This function mutates the scalars in place
    // We can make this more generic in the future to use other than u16.
    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>> {
        assert!(w > 0);
        let batch_size = scalars.len();
        let window_size: i16 = 1 << (w + 1);
        let half_window_size: i16 = 1 << w;

        let mut op_code_vectorised =
            Vec::<Vec<Option<i16>>>::with_capacity(scalars[0].as_ref().len() * 64);

        let mut all_none = false;
        while !all_none {
            let mut opcode_row = Vec::with_capacity(batch_size);

            for s in scalars.iter_mut() {
                if s.is_zero() {
                    opcode_row.push(None);
                } else {
                    let op = if s.is_odd() {
                        let mut z: i16 = (s.as_ref()[0] % (1 << (w + 1))) as i16;

                        if z < half_window_size {
                            s.sub_noborrow(&BigInt::from(z as u64));
                        } else {
                            z = z - window_size;
                            s.add_nocarry(&BigInt::from((-z) as u64));
                        }
                        z
                    } else {
                        0
                    };
                    opcode_row.push(Some(op));
                    s.div2();
                }
            }

            all_none = opcode_row.iter().all(|x| x.is_none());
            if !all_none {
                op_code_vectorised.push(opcode_row);
            }
        }
        op_code_vectorised
    }

    // This function consumes the second op as it mutates it in place
    // to prevent memory allocation
    fn batch_double_in_place(bases: &mut [Self], index: &[usize]);

    fn batch_add_in_place_same_slice(bases: &mut [Self], index: &[(usize, usize)]);

    fn batch_add_in_place(bases: &mut [Self], other: &mut [Self], index: &[(usize, usize)]);

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(
        mut bases: &mut [Self],
        scalars: &mut [BigInt],
        w: usize,
    ) {
        let opcode_vectorised = Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w);
        let tables = Self::batch_wnaf_tables(bases, w);

        // Set all points to 0;
        let zero = Self::zero();
        for p in bases.iter_mut() {
            *p = zero;
        }

        for opcode_row in opcode_vectorised.iter().rev() {
            let index_double: Vec<usize> = opcode_row
                .iter()
                .enumerate()
                .filter(|x| x.1.is_some())
                .map(|x| x.0)
                .collect();

            Self::batch_double_in_place(&mut bases, &index_double[..]);

            let mut add_ops: Vec<Self> = tables
                .iter()
                .zip(opcode_row)
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|(t, op)| {
                    let idx = op.unwrap();
                    if idx > 0 {
                        t[(idx as usize) / 2].clone()
                    } else {
                        t[((-idx) as usize) / 2].clone().neg()
                    }
                })
                .collect();

            let index_add: Vec<(usize, usize)> = opcode_row
                .iter()
                .enumerate()
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|x| x.0)
                .enumerate()
                .map(|(x, y)| (y, x))
                .collect();

            Self::batch_add_in_place(&mut bases, &mut add_ops[..], &index_add[..]);
        }
    }

    fn get_chunked_instr<T: Clone>(instr: &[T], batch_size: usize) -> Vec<Vec<T>> {
        let mut res = Vec::new();

        let rem = instr.chunks_exact(batch_size).remainder();
        let mut chunks = instr.chunks_exact(batch_size).peekable();

        if chunks.len() == 0 {
            res.push(rem.to_vec());
        }

        while let Some(chunk) = chunks.next() {
            let chunk = if chunks.peek().is_none() {
                [chunk, rem].concat()
            } else {
                chunk.to_vec()
            };
            res.push(chunk);
        }
        res
    }
}

// We make the syntax cleaner by defining corresponding trait and impl for [G]
pub trait BatchGroupArithmeticSlice<G: AffineCurve> {
    fn batch_wnaf_tables(&self, w: usize) -> Vec<Vec<G>>;

    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>>;

    fn batch_double_in_place(&mut self, index: &[usize]);

    fn batch_add_in_place_same_slice(&mut self, index: &[(usize, usize)]);

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(usize, usize)]);

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize);
}

impl<G: AffineCurve> BatchGroupArithmeticSlice<G> for [G] {
    fn batch_wnaf_tables(&self, w: usize) -> Vec<Vec<G>> {
        G::batch_wnaf_tables(self, w)
    }

    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>> {
        G::batch_wnaf_opcode_recoding::<BigInt>(scalars, w)
    }

    fn batch_double_in_place(&mut self, index: &[usize]) {
        G::batch_double_in_place(self, index);
    }

    fn batch_add_in_place_same_slice(&mut self, index: &[(usize, usize)]) {
        G::batch_add_in_place_same_slice(self, index);
    }

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(usize, usize)]) {
        G::batch_add_in_place(self, other, index);
    }

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize) {
        G::batch_scalar_mul_in_place(self, scalars, w);
    }
}

trait GLV: AffineCurve {
    fn glv_scalar_decomposition<BigInt: BigInteger, SmallBigInt: BigInteger>(
        k: BigInt,
    ) -> (SmallBigInt, SmallBigInt);

    fn glv_endomorphism_in_place(&mut self);

    fn batch_scalar_mul_in_place_glv<BigInt: BigInteger, SmallBigInt: BigInteger>(
        w: usize,
        points: &mut [Self],
        scalars: &mut [BigInt],
    ) {
        assert_eq!(points.len(), scalars.len());
        let batch_size = points.len();
        let glv_scalars: Vec<(SmallBigInt, SmallBigInt)> = scalars
            .iter()
            .map(|&s| Self::glv_scalar_decomposition::<BigInt, SmallBigInt>(s))
            .collect();
        let (mut k1, mut k2): (Vec<SmallBigInt>, Vec<SmallBigInt>) = (
            glv_scalars.iter().map(|x| x.0).collect(),
            glv_scalars.iter().map(|x| x.1).collect(),
        );

        let mut p2 = points.to_vec();
        p2.iter_mut().for_each(|p| p.glv_endomorphism_in_place());
        Self::batch_scalar_mul_in_place::<SmallBigInt>(points, &mut k1[..], w);
        Self::batch_scalar_mul_in_place::<SmallBigInt>(&mut p2[..], &mut k2[..], w);
        Self::batch_add_in_place(
            points,
            &mut p2,
            &(0..batch_size)
                .map(|x| (x, x))
                .collect::<Vec<(usize, usize)>>()[..],
        );
    }
}
