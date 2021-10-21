use super::quadratic_extension::*;
use std::marker::PhantomData;
use std::ops::{MulAssign, Neg};

use crate::{
    bits::{FromBits, FromCompressedBits, ToBits, ToCompressedBits},
    fields::{Field, Fp3, Fp3Parameters, SquareRootField},
    BitSerializationError, Error,
};

/// Model for quadratic extension field F6 as towered extension
///
//     F6 = F2[Y]/(Y^2-X),
//     F3 = Fp[X]/(X^3-alpha),
///
/// using a "non-residue" alpha mod p such that (X^6-alpha) is irreducible over Fp.
/// Its arithmetics includes pairing-relevant operations such as exponentiation and
/// squaring on the r-th unit roots of F6 (cyclotomic exp. and squ.).
pub trait Fp6Parameters: 'static + Send + Sync {
    type Fp3Params: Fp3Parameters;

    const NONRESIDUE: Fp3<Self::Fp3Params>;

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP6_C1: &'static [<Self::Fp3Params as Fp3Parameters>::Fp];

    #[inline(always)]
    fn mul_fp3_by_nonresidue(fe: &Fp3<Self::Fp3Params>) -> Fp3<Self::Fp3Params> {
        let mut res = *fe;
        res.c0 = fe.c2;
        res.c1 = fe.c0;
        res.c2 = fe.c1;
        res.c0 = <Self::Fp3Params as Fp3Parameters>::mul_fp_by_nonresidue(&res.c0);
        res
    }
}

pub struct Fp6ParamsWrapper<P: Fp6Parameters>(PhantomData<P>);

impl<P: Fp6Parameters> QuadExtParameters for Fp6ParamsWrapper<P> {
    type BasePrimeField = <P::Fp3Params as Fp3Parameters>::Fp;
    type BaseField = Fp3<P::Fp3Params>;
    type FrobCoeff = Self::BasePrimeField;

    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 6;

    const NONRESIDUE: Self::BaseField = P::NONRESIDUE;

    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff] = P::FROBENIUS_COEFF_FP6_C1;

    #[inline(always)]
    fn mul_base_field_by_nonresidue(fe: &Self::BaseField) -> Self::BaseField {
        P::mul_fp3_by_nonresidue(fe)
    }

    fn mul_base_field_by_frob_coeff(fe: &mut Self::BaseField, power: usize) {
        fe.mul_assign_by_fp(&Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD]);
    }
}

pub type Fp6<P> = QuadExtField<Fp6ParamsWrapper<P>>;

impl<P: Fp6Parameters> Fp6<P> {
    pub fn mul_by_034(
        &mut self,
        c0: &<P::Fp3Params as Fp3Parameters>::Fp,
        c3: &<P::Fp3Params as Fp3Parameters>::Fp,
        c4: &<P::Fp3Params as Fp3Parameters>::Fp,
    ) {
        let z0 = self.c0.c0;
        let z1 = self.c0.c1;
        let z2 = self.c0.c2;
        let z3 = self.c1.c0;
        let z4 = self.c1.c1;
        let z5 = self.c1.c2;

        let x0 = *c0;
        let x3 = *c3;
        let x4 = *c4;

        let mut tmp1 = x3;
        tmp1.mul_assign(&<P::Fp3Params as Fp3Parameters>::NONRESIDUE);
        let mut tmp2 = x4;
        tmp2.mul_assign(&<P::Fp3Params as Fp3Parameters>::NONRESIDUE);

        self.c0.c0 = x0 * &z0 + &(tmp1 * &z5) + &(tmp2 * &z4);
        self.c0.c1 = x0 * &z1 + &(x3 * &z3) + &(tmp2 * &z5);
        self.c0.c2 = x0 * &z2 + &(x3 * &z4) + &(x4 * &z3);
        self.c1.c0 = x0 * &z3 + &(x3 * &z0) + &(tmp2 * &z2);
        self.c1.c1 = x0 * &z4 + &(x3 * &z1) + &(x4 * &z0);
        self.c1.c2 = x0 * &z5 + &(x3 * &z2) + &(x4 * &z1);
    }

    pub fn mul_by_014(
        &mut self,
        c0: &<P::Fp3Params as Fp3Parameters>::Fp,
        c1: &<P::Fp3Params as Fp3Parameters>::Fp,
        c4: &<P::Fp3Params as Fp3Parameters>::Fp,
    ) {
        let z0 = self.c0.c0;
        let z1 = self.c0.c1;
        let z2 = self.c0.c2;
        let z3 = self.c1.c0;
        let z4 = self.c1.c1;
        let z5 = self.c1.c2;

        let x0 = *c0;
        let x1 = *c1;
        let x4 = *c4;

        let mut tmp1 = x1;
        tmp1.mul_assign(&<P::Fp3Params as Fp3Parameters>::NONRESIDUE);
        let mut tmp2 = x4;
        tmp2.mul_assign(&<P::Fp3Params as Fp3Parameters>::NONRESIDUE);

        self.c0.c0 = x0 * &z0 + &(tmp1 * &z2) + &(tmp2 * &z4);
        self.c0.c1 = x0 * &z1 + &(x1 * &z0) + &(tmp2 * &z5);
        self.c0.c2 = x0 * &z2 + &(x1 * &z1) + &(x4 * &z3);
        self.c1.c0 = x0 * &z3 + &(tmp1 * &z5) + &(tmp2 * &z2);
        self.c1.c1 = x0 * &z4 + &(x1 * &z3) + &(x4 * &z0);
        self.c1.c2 = x0 * &z5 + &(x1 * &z4) + &(x4 * &z1);
    }

    //Mul by an element of the form [c0: (0, 0, a), c1: (b, c, d)]
    pub fn mul_by_2345(self, other: &Self) -> Self
/* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */
    {
        let v0 = {
            let t = other.c0.c2 * &<P::Fp3Params as Fp3Parameters>::NONRESIDUE;
            Fp3::<P::Fp3Params>::new(self.c0.c1 * &t, self.c0.c2 * &t, self.c0.c0 * &other.c0.c2)
        };
        let v1 = self.c1 * &other.c1;
        let beta_v1 = P::mul_fp3_by_nonresidue(&v1);
        let c0 = v0 + &beta_v1;
        let c1 = (self.c0 + &self.c1) * &(other.c0 + &other.c1) - &v0 - &v1;
        Self::new(c0, c1)
    }
}

/// Note: compression and decompression of a Fqk element is possible thanks to a property of Ate pairing.
/// if c0 + i*c1 is the output of an Ate pairing, then holds that c0^2 - nr * c1^2 = 1.
/// Therefore, we can save c1 and compute c0 as sqrt(1 + nr*c1^2), dedicating a bit also for the sign
/// of the result.

impl<P: Fp6Parameters> ToCompressedBits for Fp6<P> {
    #[inline]
    fn compress(&self) -> Vec<bool> {
        //Serialize c1
        let mut res = self.c1.write_bits();

        //Set the MSB to indicate the parity of c0
        let parity = self.c0.is_odd();
        res.push(parity);

        res
    }
}

impl<P: Fp6Parameters> FromCompressedBits for Fp6<P> {
    #[inline]
    fn decompress(compressed: Vec<bool>) -> Result<Self, Error> {
        let len = compressed.len() - 1;
        let parity_flag_set = compressed[len];

        //Mask away the flag bits and try to get the c1 component
        let c1 = Fp3::read_bits(compressed[..len].to_vec())?;

        //Compute c0
        let c0 = {
            let t = Fp3::one() + &P::mul_fp3_by_nonresidue(&(c1.square()));
            t.sqrt()
        };

        match c0 {
            //Estabilish c0 parity
            Some(c0_u) => {
                let neg_c0u = c0_u.neg();
                let c0_s = if c0_u.is_odd() ^ parity_flag_set {
                    neg_c0u
                } else {
                    c0_u
                };
                Ok(Self::new(c0_s, c1))
            }

            //sqrt(1 + nr*c1^2) doesn't exists in the field
            _ => Err(Box::new(BitSerializationError::UndefinedSqrt)),
        }
    }
}
