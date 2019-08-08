use crate::{
    AffineCurve, BigInteger, Field, FpParameters, PairingCurve, PairingEngine, PrimeField,
    ProjectiveCurve,
};
use rayon::prelude::*;

pub struct VariableBaseMSM;

impl VariableBaseMSM {
    fn msm_inner<G: PairingCurve>(
        bases: &[G],
        scalars: &[<<G::Engine as PairingEngine>::Fr as PrimeField>::BigInt],
    ) -> <G as AffineCurve>::Projective {
        let c = if scalars.len() < 32 {
            3
        } else {
            (f64::from(scalars.len() as u32)).ln().ceil() as usize
        };

        let num_bits =
            <<G::Engine as PairingEngine>::Fr as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = <G::Engine as PairingEngine>::Fr::one().into_repr();

        let zero = G::zero().into_projective();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts
            .into_par_iter()
            .map(|w_start| {
                let mut res = zero;
                // We don't need the "zero" bucket, so we only have 2^c - 1 buckets
                let mut buckets = vec![zero; (1 << c) - 1];
                for (&scalar, base) in scalars.iter().zip(bases) {
                    if scalar.is_zero() {
                        continue;
                    }
                    if scalar == fr_one {
                        // We only process unit scalars once in the first window.
                        if w_start == 0 {
                            res.add_assign_mixed(base);
                        }
                    } else {
                        let mut scalar = scalar;

                        // We right-shift by w_start, thus getting rid of the
                        // lower bits.
                        scalar.divn(w_start as u32);

                        // We mod the remaining bits by the window size.
                        let scalar = scalar.as_ref()[0] % (1 << c);

                        // If the scalar is non-zero, we update the corresponding
                        // bucket.
                        // (Recall that `buckets` doesn't have a zero bucket.)
                        if scalar != 0 {
                            buckets[(scalar - 1) as usize].add_assign_mixed(&base);
                        } else {
                            continue;
                        }
                    }
                }
                let mut running_sum = G::Projective::zero();
                for b in buckets.iter().rev() {
                    running_sum += &b;
                    res += &running_sum;
                }

                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = *window_sums.first().unwrap();

        // We're traversing windows from high to low.
        let num_windows = window_sums.len();
        let mut total = zero;
        for i in (1..num_windows).rev() {
            total += &window_sums[i];
            for _ in 0..c {
                total.double_in_place();
            }
        }
        total + &lowest
    }

    pub fn multi_scalar_mul<G: PairingCurve>(
        bases: &[G],
        scalars: &[<<G::Engine as PairingEngine>::Fr as PrimeField>::BigInt],
    ) -> <G as AffineCurve>::Projective {
        Self::msm_inner(bases, scalars)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::curves::bls12_381::Bls12_381;
    use rand::{self, Rand, SeedableRng, XorShiftRng};

    fn naive_var_base_msm<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let mut acc = G::Projective::zero();

        for (base, scalar) in bases.iter().zip(scalars.iter()) {
            acc += &base.mul(*scalar);
        }
        acc
    }

    #[test]
    fn test_with_bls12() {
        const SAMPLES: usize = 1 << 10;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v = (0..SAMPLES)
            .map(|_| <Bls12_381 as PairingEngine>::Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| <Bls12_381 as PairingEngine>::G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

        assert_eq!(naive.into_affine(), fast.into_affine());
    }


    #[test]
    fn test_with_bls12_unequal_numbers() {
        const SAMPLES: usize = 1 << 10;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v = (0..SAMPLES-1)
            .map(|_| <Bls12_381 as PairingEngine>::Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| <Bls12_381 as PairingEngine>::G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

        assert_eq!(naive.into_affine(), fast.into_affine());
    }
}
