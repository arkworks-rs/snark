use crate::{biginteger::BigInteger, PrimeField, ModelParameters};

// TODO: Make GLV override slower mul
pub trait GLVParameters: Send + Sync + 'static + ModelParameters {
    type SmallBigInt: BigInteger;
    type WideBigInt: BigInteger;

    const LAMBDA: Self::ScalarField;                         // lambda in ZZ s.t. phi(P) = lambda*P for all P
    const OMEGA: Self::BaseField;                            // phi((x, y)) = (\omega x, y)
    const Q1: <Self::ScalarField as PrimeField>::BigInt;     // round(R*|b2|/n)
    const Q2: <Self::ScalarField as PrimeField>::BigInt;     // round(R*|b1|/n)
    const B1: <Self::ScalarField as PrimeField>::BigInt;     // |b1|
    const B2: <Self::ScalarField as PrimeField>::BigInt;     // |b2|
    const B1_IS_NEG: bool;

    // Not sure if all the data copying due to `from_slice` would result in a very inefficient implementation
    fn glv_scalar_decomposition(
        k: <Self::ScalarField as PrimeField>::BigInt,
    ) -> ((bool, Self::SmallBigInt), (bool, Self::SmallBigInt)) {
        let limbs = <Self::ScalarField as PrimeField>::BigInt::NUM_LIMBS;
        let modulus = Self::ScalarField::modulus();

        // We set R = 2^(NUM_LIMBS * 64)
        let mut half = Self::WideBigInt::from(1);
        half.muln((limbs as u32 * 64) - 1);

        let mut c1_wide = Self::WideBigInt::mul_no_reduce(k.as_ref(), Self::Q1.as_ref());
        // add half to achieve rounding rather than flooring
        c1_wide.add_nocarry(&half);
        // Approximation to round(|b2|*k/n)
        let c1 = &c1_wide.as_ref()[limbs..];

        let mut c2_wide = Self::WideBigInt::mul_no_reduce(k.as_ref(), Self::Q2.as_ref());
        c2_wide.add_nocarry(&half);
        let c2 = &c2_wide.as_ref()[limbs..];

        let d1 = <Self::ScalarField as PrimeField>::BigInt::mul_no_reduce_lo(&c1, Self::B1.as_ref());
        let d2 = <Self::ScalarField as PrimeField>::BigInt::mul_no_reduce_lo(&c2, Self::B2.as_ref());

        // Exactly one of B1, B2 is neg. Their
        let mut k2 = if Self::B1_IS_NEG { d2.clone() } else { d1.clone() };
        let borrow = if Self::B1_IS_NEG {
            k2.sub_noborrow(&d1)
        } else {
            k2.sub_noborrow(&d2)
        };
        let neg2 = !borrow;
        if borrow {
            k2.add_nocarry(&modulus);
        } else if k2 > modulus {
            k2.sub_noborrow(&modulus);
        }

        let mut k1 = k;
        let borrow = k2.sub_noborrow(&(Self::ScalarField::from(k1) * &Self::LAMBDA).into_repr());
        let neg1 = borrow;
        if borrow {
            k1.add_nocarry(&modulus);
        }

        let s_limbs = Self::SmallBigInt::NUM_LIMBS;

        // We should really return field elements and then let the next part of the process determine if
        let k1 = Self::SmallBigInt::from_slice(&k1.as_ref()[..s_limbs]);
        let k2 = Self::SmallBigInt::from_slice(&k2.as_ref()[..s_limbs]);

        ((neg1, k1), (neg2, k2))
    }
}

    // fn mul_glv(&self, ) {
    //
    // }

    // fn batch_scalar_mul_in_place_glv(
    //     w: usize,
    //     points: &mut [Self],
    //     scalars: &mut [<Self::Fr as PrimeField>::BigInt],
    // ) {
    //     assert_eq!(points.len(), scalars.len());
    //     let batch_size = points.len();
    //     let glv_scalars: Vec<(Self::SmallBigInt, Self::SmallBigInt)> = scalars
    //         .iter()
    //         .map(|&s| Self::glv_scalar_decomposition(s))
    //         .collect();
    //     let (mut k1, mut k2): (Vec<Self::SmallBigInt>, Vec<Self::SmallBigInt>) = (
    //         glv_scalars.iter().map(|x| x.0).collect(),
    //         glv_scalars.iter().map(|x| x.1).collect(),
    //     );
    //
    //     let mut p2 = points.to_vec();
    //     p2.iter_mut().for_each(|p| p.glv_endomorphism_in_place());
    //
    //     // THIS IS WRONG and does not achieve the savings hoped for
    //     Self::batch_scalar_mul_in_place::<Self::SmallBigInt>(points, &mut k1[..], w);
    //     Self::batch_scalar_mul_in_place::<Self::SmallBigInt>(&mut p2[..], &mut k2[..], w);
    //     Self::batch_add_in_place(
    //         points,
    //         &mut p2,
    //         &(0..batch_size)
    //             .map(|x| (x, x))
    //             .collect::<Vec<(usize, usize)>>()[..],
    //     );
    // }
