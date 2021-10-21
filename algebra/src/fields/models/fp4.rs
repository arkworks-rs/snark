use super::quadratic_extension::*;
use std::marker::PhantomData;

use crate::{
    bits::{FromBits, FromCompressedBits, ToBits, ToCompressedBits},
    fields::{Field, Fp2, Fp2Parameters, SquareRootField},
    BitSerializationError, Error,
};

/// Model for quadratic extension field F4 as towered extension
///
//     F4 = F2[Y]/(Y^2-X),
//     F2 = Fp[X]/(X^2-alpha),
///
/// using a "non-residue" alpha mod p such that (X^4-alpha) is irreducible over Fp.
/// Its arithmetics includes pairing-relevant operations such as exponentiation and
/// squaring on the r-th unit roots of F4 (cyclotomic exp. and squ.).
pub trait Fp4Parameters: 'static + Send + Sync {
    type Fp2Params: Fp2Parameters;

    /// This *must* equal (0, 1);
    /// see [[DESD06, Section 5.1]](https://eprint.iacr.org/2006/471.pdf).
    const NONRESIDUE: Fp2<Self::Fp2Params>;

    /// Coefficients for the Frobenius automorphism.
    /// non_residue^((modulus^i-1)/4) for i=0,1,2,3
    const FROBENIUS_COEFF_FP4_C1: &'static [<Self::Fp2Params as Fp2Parameters>::Fp];

    #[inline(always)]
    fn mul_fp2_by_nonresidue(fe: &Fp2<Self::Fp2Params>) -> Fp2<Self::Fp2Params> {
        // see [[DESD06, Section 5.1]](https://eprint.iacr.org/2006/471.pdf).
        Fp2::new(
            <Self::Fp2Params as Fp2Parameters>::NONRESIDUE * &fe.c1,
            fe.c0,
        )
    }
}

pub struct Fp4ParamsWrapper<P: Fp4Parameters>(PhantomData<P>);

impl<P: Fp4Parameters> QuadExtParameters for Fp4ParamsWrapper<P> {
    type BasePrimeField = <P::Fp2Params as Fp2Parameters>::Fp;
    type BaseField = Fp2<P::Fp2Params>;
    type FrobCoeff = Self::BasePrimeField;

    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 4;

    const NONRESIDUE: Self::BaseField = P::NONRESIDUE;

    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff] = P::FROBENIUS_COEFF_FP4_C1;

    #[inline(always)]
    fn mul_base_field_by_nonresidue(fe: &Self::BaseField) -> Self::BaseField {
        P::mul_fp2_by_nonresidue(fe)
    }

    fn mul_base_field_by_frob_coeff(fe: &mut Self::BaseField, power: usize) {
        fe.mul_assign_by_fp(&Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD]);
    }

    fn cyclotomic_square(fe: &QuadExtField<Self>) -> QuadExtField<Self> {
        let a = fe.c1.square();
        let b = fe.c1 + &fe.c0;
        let c = b.square() - &a;
        let d = Self::mul_base_field_by_nonresidue(&a);
        let e = c - &d;
        QuadExtField::<Self>::new(
            d.double() + &Self::BaseField::one(),
            e - &Self::BaseField::one(),
        )
    }
}

pub type Fp4<P> = QuadExtField<Fp4ParamsWrapper<P>>;

impl<P: Fp4Parameters> Fp4<P> {
    pub fn mul_by_fp(&mut self, element: &<P::Fp2Params as Fp2Parameters>::Fp) {
        self.c0.mul_assign_by_fp(element);
        self.c1.mul_assign_by_fp(element);
    }

    pub fn mul_by_fp2(&mut self, element: &Fp2<P::Fp2Params>) {
        self.c0 *= element;
        self.c1 *= element;
    }

    //Mul by an element of the form (c0: [c0, 0] c1: [c2, c3])
    pub fn mul_by_023(self, other: &Self) -> Self {
        let v0 = {
            let v0_c0 = self.c0.c0 * &other.c0.c0;
            let v0_c1 = self.c0.c1 * &other.c0.c0;
            Fp2::new(v0_c0, v0_c1)
        };
        let v1 = self.c1 * &other.c1;

        let c0 = v0 + &P::mul_fp2_by_nonresidue(&v1);
        let c1 = (self.c0 + &self.c1) * &(other.c0 + &other.c1) - &v0 - &v1;

        Self::new(c0, c1)
    }
}

/// Note: compression and decompression of a Fqk element is possible thanks to a property of Ate pairing.
/// if c0 + i*c1 is the output of an Ate pairing, then holds that c0^2 - nr * c1^2 = 1.
/// Therefore, we can save c1 and compute c0 as sqrt(1 + nr*c1^2), dedicating a bit also for the sign
/// of the result.

impl<P: Fp4Parameters> ToCompressedBits for Fp4<P> {
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

impl<P: Fp4Parameters> FromCompressedBits for Fp4<P> {
    #[inline]
    fn decompress(compressed: Vec<bool>) -> Result<Self, Error> {
        let len = compressed.len() - 1;
        let parity_flag_set = compressed[len];

        //Mask away the flag bits and try to get the c1 component
        let c1 = Fp2::read_bits(compressed[..len].to_vec())?;

        //Compute c0
        let c0 = {
            let t = Fp2::one() + &P::mul_fp2_by_nonresidue(&(c1.square()));
            t.sqrt()
        };

        match c0 {
            //Estabilish c0 parity
            Some(c0_u) => {
                let c0_s = if c0_u.is_odd() ^ parity_flag_set {
                    -c0_u
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
