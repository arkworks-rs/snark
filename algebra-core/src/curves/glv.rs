use crate::{biginteger::BigInteger, ModelParameters, PrimeField};
use core::ops::Neg;

/// The GLV parameters here require the following conditions to be satisfied:
/// 1. MODULUS_BITS < NUM_LIMBS * 64 - 1. So 2 * n < 1 << (64 * NUM_LIMBS)
/// We also assume that |b1| * |b2| < 2 * n
/// We also know that either B1 is neg or B2 is.
pub trait GLVParameters: Send + Sync + 'static + ModelParameters {
    type WideBigInt: BigInteger;

    const LAMBDA: Self::ScalarField; // lambda in ZZ s.t. phi(P) = lambda*P for all P
    const OMEGA: Self::BaseField; // phi((x, y)) = (\omega x, y)
    const Q1: <Self::ScalarField as PrimeField>::BigInt; // round(R*|b2|/n)
    const Q2: <Self::ScalarField as PrimeField>::BigInt; // round(R*|b1|/n)
    const B1: <Self::ScalarField as PrimeField>::BigInt; // |b1|
    const B2: <Self::ScalarField as PrimeField>::BigInt; // |b2|
    const B1_IS_NEG: bool;
    const R_BITS: u32;

    #[inline]
    fn glv_scalar_decomposition_inner(
        k: <Self::ScalarField as PrimeField>::BigInt,
    ) -> (
        (bool, <Self::ScalarField as PrimeField>::BigInt),
        (bool, <Self::ScalarField as PrimeField>::BigInt),
    ) {
        let limbs = <Self::ScalarField as PrimeField>::BigInt::NUM_LIMBS;
        let modulus = Self::ScalarField::modulus();

        // If we are doing a subgroup check, we should multiply by the original scalar
        // since the GLV decomposition does not guarantee that we would not be
        // adding and subtracting back to zero
        if k == modulus {
            return (
                (false, k),
                (false, <Self::ScalarField as PrimeField>::BigInt::from(0)),
            );
        }

        let mut half = Self::WideBigInt::from(1);
        half.muln(Self::R_BITS - 1);

        let mut c1_wide = Self::WideBigInt::mul_no_reduce(k.as_ref(), Self::Q1.as_ref());
        // add half to achieve rounding rather than flooring
        c1_wide.add_nocarry(&half);
        // Approximation to round(|b2|*k/n)
        c1_wide.divn(Self::R_BITS);
        let c1 = &c1_wide.as_ref()[..limbs];

        let mut c2_wide = Self::WideBigInt::mul_no_reduce(k.as_ref(), Self::Q2.as_ref());
        c2_wide.add_nocarry(&half);
        c2_wide.divn(Self::R_BITS);
        let c2 = &c2_wide.as_ref()[..limbs];

        // We first assume that the final 2 bits of the representation for the modulus
        // is not set, so that 2 * n < R = 1 << (64 * NUM_LIMBS). Then, since we
        // know that |b_i| < \sqrt{2n}, wlog k|b1|/n * |b2| < 2 * k <  2 * n <
        // R.
        let d1 =
            <Self::ScalarField as PrimeField>::BigInt::mul_no_reduce_lo(&c1, Self::B1.as_ref());
        let d2 =
            <Self::ScalarField as PrimeField>::BigInt::mul_no_reduce_lo(&c2, Self::B2.as_ref());

        // We check if they have the same sign. If they do, we must do a subtraction.
        // Else, we must do an addition. Then, we will conditionally add or
        // subtract the product of this with lambda from k. We do this to obtain
        // the result k_2 = -(c1.b1 + c1.b1) = sign(b1)*(c2|b2| - c1|b1|) = sign(b1)(d2
        // - d1)
        let mut k2_field = if Self::B1_IS_NEG {
            Self::ScalarField::from(d2)
        } else {
            Self::ScalarField::from(d1)
        };
        if Self::B1_IS_NEG {
            k2_field -= &Self::ScalarField::from(d1);
        } else {
            k2_field -= &Self::ScalarField::from(d2);
        }

        let k1 = (Self::ScalarField::from(k) - &(k2_field * &Self::LAMBDA)).into_repr();

        let k2 = k2_field.into_repr();

        let (neg2, k2) = if k2.num_bits() > Self::R_BITS / 2 + 1 {
            (true, k2_field.neg().into_repr())
        } else {
            (false, k2)
        };

        let (neg1, k1) = if k1.num_bits() > Self::R_BITS / 2 + 1 {
            (true, Self::ScalarField::from(k1).neg().into_repr())
        } else {
            (false, k1)
        };

        ((neg1, k1), (neg2, k2))
    }
}
