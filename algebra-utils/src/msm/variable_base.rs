use algebra::{AffineCurve, BigInteger, Field, FpParameters, PrimeField, ProjectiveCurve};
use rayon::prelude::*;

use std::any::TypeId;

#[cfg(feature = "gpu")]
use algebra_kernels::msm::{get_cpu_utilization, get_kernels, get_gpu_min_length};
#[cfg(feature = "gpu")]
use std::sync::{Arc, Mutex};

pub struct VariableBaseMSM;

impl VariableBaseMSM {

    // Function that recodes the scalars into SD numbers
    // The output is a vector
    pub fn recode_sd<G: AffineCurve>(
        scalar: &<G::ScalarField as PrimeField>::BigInt,
        c: usize
    ) -> Vec<i64> {

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;

        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        let mut vec_coeff = Vec::new();

        window_starts.iter().rev().for_each(|x| {
            let mut scal = (*scalar).clone();
            scal.divn(*x as u32);
            // We mod the remaining bits by the window size.
            let a = scal.as_ref()[0] % (1 << c);
            vec_coeff.push(a as i64);
        });

        for idx in (0..vec_coeff.len()).rev() {
            if vec_coeff[idx] >= (1 << (c-1)) {
                vec_coeff[idx] -= 1 << c;
                if idx!= 0 {
                    vec_coeff[idx-1] += 1;
                }
            }
        }

        // last element is the least significant digit
        return vec_coeff;

    }

    // Recodes in SD the portions of the scalar divided in c bits
    // Takes as input the carry_in and generates the vector of SD digits corresponding to
    // the processing chunk and the carry_out
    // For example: for 4 cores, it takes the carry_in and
    // generates carry_out, a_(i+3), a_(i+2), a_(i+1), a_i
    #[inline]
    pub fn recode_sd_chunk<G:AffineCurve>(
        scalar: &<G::ScalarField as PrimeField>::BigInt,
        c: usize,
        chunk_pos: usize,
        num_cores: usize,
        vec_coeff: &mut Vec<i64>,
        carry_in: &mut i64
    )   {

        // starting bit position of the chunk
        let start_w_bit = chunk_pos * (c * num_cores);

        for i in 0..num_cores {
            let windows_starts = start_w_bit + i * c;
            let mut scal = (*scalar).clone();
            scal.divn(windows_starts as u32);
            let mut a = (scal.as_ref()[0] % (1 << c)) as i64 + *carry_in;
            if a >= (1 << (c-1)) {
                a -= 1 << c;
                *carry_in = 1;
            } else {
                *carry_in = 0;
            }
            (*vec_coeff)[i] = a;
        }
    }

    //TODO: Seems that this function has problems not captured by the UTs.
    //      Thus for the moment its invokation is disabled.
    #[allow(unused)]
    pub fn multi_scalar_mul_affine_sd_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c:usize
    ) -> G::Projective {

        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;

        let zero = G::zero().into_projective();

        // The number of windows of the scalars in groups of c bits
        let num_w = (num_bits as f64 / c as f64).ceil() as usize;
        // The number of cores present
        let cpus = rayon::current_num_threads();
        // The number of chunks of cpus digits
        let num_chunks = (num_w as f64/cpus as f64).floor() as usize;
        // The remaining digits to process sequentially
        let remaining_digits = num_w - (num_chunks * cpus);

        let mut window_sums =  Vec::new();
        let mut small_window_sums:Vec<_>;

        let mut carry_vec = vec![0; scalars.len()];
        let mut vec_coeff = vec![vec![0; cpus]; scalars.len()];

        for i in 0..num_chunks {
            // index of the digit within the small window of cpus number of digits
            let idx: Vec<_> = (0..cpus).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, i, cpus, v1, c1);
            });

            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }

                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

        if remaining_digits != 0 {

            let idx:Vec<_> = (0..remaining_digits).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, num_chunks, cpus, v1, c1);
            });

            // Each window is of size `c`.
            // We divide up the bits 0..num_bits into windows of size `c`, and
            // in parallel process each such window.
            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }
                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

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


    pub fn multi_scalar_mul_affine_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c: usize
    ) -> G::Projective {

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

    pub fn msm_inner_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c:usize
    ) -> G::Projective {

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
        let scal_len = scalars.len();

        #[cfg(feature = "bn_382")]
        if TypeId::of::<G>() == TypeId::of::<algebra::curves::bn_382::G1Affine>()
        {
            let c: usize = if scal_len < 32 {
                3
            } else if (scal_len < 1 << 17) || (scal_len > 1 << 23) {
                (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
            } else if scal_len < 1 << 21 {
                12
            } else {
                16
            };
            return Self::multi_scalar_mul_affine_c(bases, scalars, c);
        }

        #[cfg(feature = "tweedle")]
        if TypeId::of::<G>() == TypeId::of::<algebra::curves::tweedle::dee::Affine>() {
            if scal_len < 1 << 17 {
                let c: usize = if scal_len < 32 {
                    3
                } else {
                    (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
                };
                return Self::multi_scalar_mul_affine_c(bases, scalars, c);
            } else if scal_len < 1 << 19 {
                let c: usize = 11;
                return Self::multi_scalar_mul_affine_c(bases, scalars, c);
            } else if scal_len < 1 << 23 {
                let c: usize = 11;
                return Self::multi_scalar_mul_affine_c(bases, scalars, c);
            } else if scal_len < 1 << 24 {
                let c: usize = 16;
                return Self::msm_inner_c(bases, scalars, c);
            } else {
                let c: usize = (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize;
                return Self::msm_inner_c(bases, scalars, c);
            }
        }

        let c: usize = if scal_len < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        };
        return Self::multi_scalar_mul_affine_c(bases, scalars, c);

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

        let results = Arc::new(Mutex::new(vec![]));

        rayon::scope(|s| {

            if n > 0 {
                for ((bases, scalars), kern) in bases
                    .chunks(chunk_size)
                    .zip(scalars.chunks(chunk_size))
                    .zip(kernels.iter())
                {
                    let results = Arc::clone(&results);
                    s.spawn(
                        move |_| {
                            let mut acc = G::Projective::zero();
                            for (bases, scalars) in bases.chunks(kern.n).zip(scalars.chunks(kern.n)) {
                                let result = kern.msm(bases, scalars, bases.len());
                                acc.add_assign_mixed(&result.unwrap().into_affine());
                            }
                            results.lock().unwrap().push(acc);
                        },
                    );
                }
            }

            if cpu_n > 0 {
                let results = Arc::clone(&results);
                s.spawn(
                    move |_| {
                        let acc = Self::msm_inner_cpu(cpu_bases, cpu_scalars);
                        results.lock().unwrap().push(acc);
                    }
                );
            }
        });

        let mut acc = G::Projective::zero();
        let results = Arc::try_unwrap(results).unwrap().lock().unwrap().clone();

        for r in results.into_iter() {
            acc.add_assign_mixed(&r.into_affine());
        }

        acc
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

    use algebra::curves::bn_382::G1Projective as Bn382G1Projective;
    use algebra::curves::bn_382::g::Projective as Bn382GProjective;
    use algebra::curves::tweedle::dee::Projective as TweedleDee;
    use algebra::curves::tweedle::dum::Projective as TweedleDum;
    use algebra::curves::bls12_381::G1Projective as BlsG1Projective;

    use rand::{SeedableRng, Rng};
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

    fn test_all_variants<G: ProjectiveCurve, R: Rng>(
        samples: usize,
        c: usize,
        rng: &mut R,
    ) {
        let v = (0..samples)
            .map(|_| G::ScalarField::rand(rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..samples)
            .map(|_| G::rand(rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let fast = VariableBaseMSM::msm_inner_cpu(g.as_slice(), v.as_slice());

        let affine = VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(), c);
        //let affine_sd = VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(), c);
        let inner = VariableBaseMSM::msm_inner_c(g.as_slice(), v.as_slice(), c);

        #[cfg(feature = "gpu")]
        let gpu = VariableBaseMSM::msm_inner_gpu(g.as_slice(), v.as_slice());

        assert_eq!(naive, fast);

        assert_eq!(naive, affine);
        //assert_eq!(naive, affine_sd);
        assert_eq!(naive, inner);

        #[cfg(feature = "gpu")]
        assert_eq!(naive, gpu);
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