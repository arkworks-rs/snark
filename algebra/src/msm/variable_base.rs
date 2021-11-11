use crate::{AffineCurve, BigInteger, Error, Field, FpParameters, PrimeField, ProjectiveCurve};
use rayon::prelude::*;

pub struct VariableBaseMSM;

impl VariableBaseMSM {
    /// WARNING: This function allows scalars and bases to have different length
    /// (as long as scalars.len() <= bases.len()): internally, bases are trimmed
    /// to have the same length of the scalars; this may lead to potential message
    /// malleability issue: e.g. MSM([s1, s2], [b1, b2]) == MSM([s1, s2], [b1, b2, b3]),
    /// so use this function carefully.
    pub fn multi_scalar_mul_affine_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c: usize,
    ) -> Result<G::Projective, Error> {
        // Sanity checks
        if c == 0 {
            Err(format!("Invalid window size value: 0"))?
        }
        if c > 25 {
            Err(format!(
                "Invalid window size value: {}. It must be smaller than 25",
                c
            ))?
        }
        if scalars.len() > bases.len() {
            Err(format!(
                "Invalid MSM length. Scalars len: {}, Bases len: {}",
                scalars.len(),
                bases.len()
            ))?
        }

        let cc = 1 << c;

        let num_bits = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

        let zero = G::zero().into_projective();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts
            .into_par_iter()
            .map(|w_start| {
                // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc];
                scalars
                    .iter()
                    .zip(bases)
                    .filter(|(s, _)| !s.is_zero())
                    .for_each(|(&scalar, base)| {
                        if scalar == fr_one {
                            // We only process unit scalars once in the first window.
                            if w_start == 0 && base.is_zero() == false {
                                buckets[cc - 1].push(*base);
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
                            if scalar != 0 && base.is_zero() == false {
                                buckets[(scalar - 1) as usize].push(*base);
                            }
                        }
                    });
                G::add_points(&mut buckets);
                let mut res = if buckets[cc - 1].len() == 0 {
                    zero
                } else {
                    buckets[cc - 1][0].into_projective()
                };

                let mut running_sum = zero;
                for b in buckets[0..cc - 1].iter_mut().rev() {
                    if b.len() != 0 && b[0].is_zero() == false {
                        running_sum.add_assign_mixed(&b[0])
                    }
                    res += &running_sum;
                }
                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        let result = window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, sum_i| {
                total += sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            })
            + lowest;

        Ok(result)
    }

    /// WARNING: This function allows scalars and bases to have different length
    /// (as long as scalars.len() <= bases.len()): internally, bases are trimmed
    /// to have the same length of the scalars; this may lead to potential message
    /// malleability issue: e.g. MSM([s1, s2], [b1, b2]) == MSM([s1, s2], [b1, b2, b3]),
    /// so use this function carefully.
    pub fn msm_inner_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c: usize,
    ) -> Result<G::Projective, Error> {
        // Sanity checks
        if c == 0 {
            Err(format!("Invalid window size value: 0"))?
        }
        if c > 25 {
            Err(format!(
                "Invalid window size value: {}. It must be smaller than 25",
                c
            ))?
        }
        if scalars.len() > bases.len() {
            Err(format!(
                "Invalid MSM length. Scalars len: {}, Bases len: {}",
                scalars.len(),
                bases.len()
            ))?
        }

        let num_bits = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

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
                scalars
                    .iter()
                    .zip(bases)
                    .filter(|(s, _)| !s.is_zero())
                    .for_each(|(&scalar, base)| {
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
                            }
                        }
                    });
                G::Projective::batch_normalization(&mut buckets);

                let mut running_sum = G::Projective::zero();
                for b in buckets.into_iter().map(|g| g.into_affine()).rev() {
                    running_sum.add_assign_mixed(&b);
                    res += &running_sum;
                }

                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        let result = window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, sum_i| {
                total += sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            })
            + lowest;

        Ok(result)
    }

    pub fn msm_inner<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> Result<G::Projective, Error> {
        let scal_len = scalars.len();

        let c: usize = if scal_len < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        };

        Self::msm_inner_c(bases, scalars, c)
    }

    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> Result<G::Projective, Error>
    where
        G::Projective: ProjectiveCurve<Affine = G>,
    {
        let scal_len = scalars.len();

        let c: usize = if scal_len < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        };

        Self::multi_scalar_mul_affine_c(bases, scalars, c)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::curves::bls12_381::G1Projective as BlsG1Projective;
    use crate::curves::bn_382::g::Projective as Bn382GProjective;
    use crate::curves::bn_382::G1Projective as Bn382G1Projective;
    use crate::curves::tweedle::dee::Projective as TweedleDee;
    use crate::curves::tweedle::dum::Projective as TweedleDum;

    use crate::UniformRand;
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    fn naive_var_base_msm<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let mut acc = <G::Projective as ProjectiveCurve>::zero();

        for (base, scalar) in bases.iter().zip(scalars.iter()) {
            acc += &base.mul(*scalar);
        }
        acc
    }

    fn test_all_variants<G: ProjectiveCurve, R: Rng>(samples: usize, c: usize, rng: &mut R) {
        let v = (0..samples)
            .map(|_| G::ScalarField::rand(rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..samples)
            .map(|_| G::rand(rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let fast = VariableBaseMSM::msm_inner(g.as_slice(), v.as_slice()).unwrap();

        let affine =
            VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(), c).unwrap();
        let inner = VariableBaseMSM::msm_inner_c(g.as_slice(), v.as_slice(), c).unwrap();

        assert_eq!(naive, fast);

        assert_eq!(naive, affine);
        assert_eq!(naive, inner);
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn test_all_variants_tweedle() {
        let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

        test_all_variants::<TweedleDee, _>(1 << 12, 16, rng);
        test_all_variants::<TweedleDum, _>(1 << 12, 16, rng);
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn test_all_variants_bn382() {
        let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

        test_all_variants::<Bn382G1Projective, _>(1 << 12, 16, rng);
        test_all_variants::<Bn382GProjective, _>(1 << 12, 16, rng);
    }

    #[cfg(feature = "bls12_381")]
    #[test]
    fn test_all_variants_bls() {
        let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

        test_all_variants::<BlsG1Projective, _>(1 << 12, 16, rng);
    }
}
