use algebra::{AffineCurve, BigInteger, Field, FpParameters, PrimeField, ProjectiveCurve};
use rayon::prelude::*;

#[cfg(feature = "gpu")]
use algebra_kernels::msm::{get_cpu_utilization, get_kernels, get_gpu_min_length};
#[cfg(feature = "gpu")]
use algebra_cl_gen::gpu::GPUError;
#[cfg(feature = "gpu")]
use crossbeam::thread;

pub struct VariableBaseMSM;

impl VariableBaseMSM {

    pub fn multi_scalar_mul_affine<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let c = if scalars.len() < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        };
        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
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
                let mut buckets = vec![Vec::with_capacity(bases.len()/cc * 2); cc];
                scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero()).for_each(|(&scalar, base)|  {
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
        window_sums[1..].iter().rev().fold(zero, |mut total, sum_i| {
            total += sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    pub fn multi_scalar_mul_mixed<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let c = if scalars.len() < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() + 2.0).ceil() as usize
        };

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
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
                scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero()).for_each(|(&scalar, base)|  {
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
        window_sums[1..].iter().rev().fold(zero, |mut total, sum_i| {
            total += sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    fn msm_inner_cpu<G>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt]
    ) -> G::Projective
    where
        G: AffineCurve,
        G::Projective: ProjectiveCurve<Affine = G>
    {
        #[cfg(not(feature = "msm_affine"))]
        return Self::multi_scalar_mul_mixed(bases, scalars);

        #[cfg(feature = "msm_affine")]
        return Self::multi_scalar_mul_affine(bases, scalars);   
    }

    #[cfg(feature = "gpu")]
    fn msm_inner_gpu<G>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt]
    ) -> G::Projective
    where
        G: AffineCurve,
        G::Projective: ProjectiveCurve<Affine = G>
    {
        let zero = G::Projective::zero();

        let mut n = bases.len();
        let cpu_n;

        let kernels = get_kernels().unwrap();
        let num_devices = kernels.len();
        let gpu_min_length = get_gpu_min_length();
        
        if gpu_min_length > n {
            cpu_n = n;
            n = 0;
        } else {
            cpu_n = ((n as f64) * get_cpu_utilization()) as usize;
            n = n - cpu_n;                
        }

        let (cpu_bases, bases) = bases.split_at(cpu_n);
        let (cpu_scalars, scalars) = scalars.split_at(cpu_n);

        let chunk_size = ((n as f64) / (num_devices as f64)).ceil() as usize;

        match thread::scope(|s| -> Result<G::Projective, GPUError> {
            let mut acc = G::Projective::zero();
            let mut threads = Vec::new();

            if n > 0 {
                for ((bases, scalars), kern) in bases
                    .chunks(chunk_size)
                    .zip(scalars.chunks(chunk_size))
                    .zip(kernels.iter())
                {
                    threads.push(s.spawn(
                        move |_| -> Result<G::Projective, GPUError> {
                            let mut acc = G::Projective::zero();
                            for (bases, scalars) in bases.chunks(kern.n).zip(scalars.chunks(kern.n)) {
                                let result = kern.msm(bases, scalars, bases.len())?;
                                acc.add_assign_mixed(&result.into_affine());
                            }
                            Ok(acc)
                        },
                    ));
                }
            }

            if cpu_n > 0 {
                threads.push(s.spawn(
                    move |_| -> Result<G::Projective, GPUError> {
                        let acc = Self::msm_inner_cpu(cpu_bases, cpu_scalars);
                        Ok(acc)
                    }
                ))
            }

            let mut results = vec![];
            for t in threads {
                results.push(t.join());
            }
            for r in results {
                acc.add_assign_mixed(&r??.into_affine());
            }

            Ok(acc)
        }) {
            Ok(res) => res.unwrap(),
            Err(_) => zero
        }
    }

    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {

        #[cfg(not(feature = "gpu"))]
        return Self::msm_inner_cpu(bases, scalars);

        #[cfg(feature = "gpu")]
        return Self::msm_inner_gpu(bases, scalars);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use algebra::curves::bls12_381::G1Projective;
    use algebra::fields::bls12_381::Fr;
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use algebra::UniformRand;

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

    #[test]
    fn test_all_variants() {
        const SAMPLES: usize = 1 << 12;

        let mut rng = XorShiftRng::seed_from_u64(234872845u64);

        let v = (0..SAMPLES)
            .map(|_| Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let mixed = VariableBaseMSM::multi_scalar_mul_mixed(g.as_slice(), v.as_slice());
        let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());

        #[cfg(feature = "gpu")]
        let gpu = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

        assert_eq!(naive, mixed);
        assert_eq!(naive, affine);

        #[cfg(feature = "gpu")]
        assert_eq!(naive, gpu);
    }

    #[test]
    fn test_with_unequal_numbers() {
        const SAMPLES: usize = 1 << 12;

        let mut rng = XorShiftRng::seed_from_u64(234872845u64);

        let v = (0..SAMPLES-1)
            .map(|_| Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let mixed = VariableBaseMSM::multi_scalar_mul_mixed(g.as_slice(), v.as_slice());
        let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());

        #[cfg(feature = "gpu")]
        let gpu = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

        assert_eq!(naive, mixed);
        assert_eq!(naive, affine);

        #[cfg(feature = "gpu")]
        assert_eq!(naive, gpu);
    }
}
