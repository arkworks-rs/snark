macro_rules! impl_run_kernel {
    () => {
        // We drop a lock only after the parallel portion has been handled
        fn par_run_kernel_sync<T>(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
            lock: T,
        ) -> DeviceMemory<G> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let mut tables_h = vec![G::zero(); n * TABLE_SIZE];
            let mut exps_recode_h = vec![0u8; n * NUM_U8];

            let now = std::time::Instant::now();
            generate_tables_and_recoding(
                bases_h,
                &mut tables_h[..],
                exps_h,
                &mut exps_recode_h[..],
                true,
            );
            drop(lock);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );

            let now = std::time::Instant::now();
            let mut out = DeviceMemory::<G>::zeros(&ctx, n);
            let mut tables = DeviceMemory::<G>::zeros(&ctx, n * TABLE_SIZE);
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            tables.copy_from_slice(&tables_h);
            exps.copy_from_slice(&exps_recode_h);
            println!("Copied data to device: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            scalar_mul_kernel::scalar_mul(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr(), n as isize),
            )
            .expect("Kernel call failed");

            println!("Ran kernel: {}us", now.elapsed().as_micros());
            out
        }

        pub fn par_run_kernel(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
        ) -> DeviceMemory<G> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let now = std::time::Instant::now();
            let mut tables = DeviceMemory::<G>::zeros(&ctx, n * TABLE_SIZE);
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
            let mut out = DeviceMemory::<G>::zeros(&ctx, n);
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            generate_tables_and_recoding(bases_h, &mut tables[..], exps_h, &mut exps[..], true);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );
            // Accessible from CPU as usual Rust slice (though this will be slow)
            // Can this be changed to a memcpy?
            scalar_mul_kernel::scalar_mul(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr(), n as isize),
            )
            .expect("Kernel call failed");
            out
        }
    };
}

#[macro_export]
macro_rules! impl_scalar_mul_kernel {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            pub mod [<$curve _ $type _scalar_mul_kernel>] {
                use accel::*;
                use rayon::prelude::*;
                use std::sync::Mutex;
                use lazy_static::lazy_static;

                use algebra_core::{
                    biginteger::BigInteger, FpParameters, Zero,
                    curves::{ProjectiveCurve, AffineCurve, BatchGroupArithmeticSlice},
                    fields::PrimeField,
                };

                use algebra::$curve::$ProjCurve;

                pub type G = $ProjCurve;
                type PrimeF = <G as ProjectiveCurve>::ScalarField;
                pub type BigInt = <PrimeF as PrimeField>::BigInt;

                const NUM_BITS: usize = <<PrimeF as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
                const LOG2_W: usize = 5;
                const TABLE_SIZE: usize = 1 << LOG2_W;
                const NUM_U8: usize = (NUM_BITS - 1) / LOG2_W + 1;

                impl_run_kernel!();
                impl_gpu_cpu_run_kernel!([<$curve _ $type>]);

                fn scalar_recode(k: &mut BigInt) -> [u8; NUM_U8] {
                    let mut out = [0; NUM_U8];
                    for i in (0..NUM_U8).rev() {
                        out[i] = (k.as_ref()[0] % TABLE_SIZE as u64) as u8;
                        k.divn(LOG2_W as u32);
                    }
                    assert!(k.is_zero());
                    out
                }

                fn generate_tables_and_recoding(
                    bases_h: &[<G as ProjectiveCurve>::Affine],
                    tables_h: &mut [G],
                    exps_h: &[BigInt],
                    exps_recode_h: &mut [u8],
                    run_parallel: bool,
                ) {
                    let closure = |
                        ((k, exps_chunk), (table, base)):
                        ((&BigInt, &mut [u8]), (&mut [G], &<G as ProjectiveCurve>::Affine))
                    | {
                        let base = base.into_projective();
                        exps_chunk.clone_from_slice(&scalar_recode(&mut k.clone()));

                        table[0] = G::zero();
                        for i in 1..TABLE_SIZE {
                            table[i] = table[i - 1] + base;
                        }
                    };
                    if run_parallel {
                        exps_h
                            .par_iter()
                            .zip(exps_recode_h.par_chunks_mut(NUM_U8))
                            .zip(tables_h.par_chunks_mut(TABLE_SIZE).zip(bases_h.par_iter()))
                            .for_each(|x| closure(x));
                    } else {
                        exps_h
                            .iter()
                            .zip(exps_recode_h.chunks_mut(NUM_U8))
                            .zip(tables_h.chunks_mut(TABLE_SIZE).zip(bases_h.iter()))
                            .for_each(|x| closure(x));
                    }
                }

                #[kernel_mod(to_mod)]
                #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
                #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
                #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
                pub mod scalar_mul {
                    use algebra::{$curve::$ProjCurve, FpParameters, Zero};
                    use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};

                    const NUM_BITS: isize =
                        <<<$ProjCurve as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
                    const LOG2_W: isize = 5;
                    const TABLE_SIZE: isize = 1 << LOG2_W;
                    const HALF_TABLE_SIZE: isize = 1 << (LOG2_W - 1);
                    const NUM_U8: isize = (NUM_BITS - 1) / LOG2_W + 1;

                    #[kernel_func]
                    pub unsafe fn scalar_mul(
                        #[type_substitute(*const $crate::[<$curve _ $type _scalar_mul_kernel>]::G)]
                        table: *const $ProjCurve,
                        exps: *const u8,
                        #[type_substitute(*mut $crate::[<$curve _ $type _scalar_mul_kernel>]::G)]
                        out: *mut $ProjCurve,
                        n: isize,
                    ) {
                        let i = accel_core::index();
                        if i < n {
                            let mut res = $ProjCurve::zero();
                            res += &(*table.offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8) as isize));

                            for j in 1..NUM_U8 as isize {
                                for _ in 0..LOG2_W {
                                    res.double_in_place();
                                }
                                res += &(*table
                                    .offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8 + j) as isize));
                            }
                            *out.offset(i) = res;
                        }
                    }
                }
            }
        }
    }
}

#[macro_export]
macro_rules! impl_scalar_mul_kernel_glv {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            pub mod [<$curve _ $type _scalar_mul_kernel>] {
                use accel::*;
                use rayon::prelude::*;
                use std::sync::Mutex;
                use lazy_static::lazy_static;

                use algebra_core::{
                    biginteger::BigInteger, FpParameters, Zero,
                    curves::{ProjectiveCurve, AffineCurve, BatchGroupArithmeticSlice},
                    fields::PrimeField,
                };
                use std::ops::Neg;

                use algebra::$curve::$ProjCurve;

                pub type G = $ProjCurve;
                type PrimeF = <G as ProjectiveCurve>::ScalarField;
                pub type BigInt = <PrimeF as PrimeField>::BigInt;

                const NUM_BITS: usize = <<PrimeF as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
                const LOG2_W: usize = 5;
                const TABLE_SIZE: usize = 1 << LOG2_W;
                const NUM_U8: usize = 2 * ((NUM_BITS - 1) / (2 * (LOG2_W - 1)) + 2);

                impl_run_kernel!();
                impl_gpu_cpu_run_kernel!([<$curve _ $type>]);

                fn scalar_recode_glv(k1: &mut BigInt, k2: &mut BigInt) -> [u8; NUM_U8] {
                    const TABLE_SIZE_GLV: u64 = 1u64 << (LOG2_W - 1);
                    let mut out = [0; NUM_U8];
                    for i in (0..NUM_U8 / 2).rev() {
                        out[2 * i] = (k1.as_ref()[0] % TABLE_SIZE_GLV) as u8;
                        out[2 * i + 1] = (k2.as_ref()[0] % TABLE_SIZE_GLV) as u8;
                        k1.divn(LOG2_W as u32 - 1);
                        k2.divn(LOG2_W as u32 - 1);
                    }
                    assert!(k1.is_zero());
                    assert!(k2.is_zero());
                    out
                }

                fn generate_tables_and_recoding(
                    bases_h: &[<G as ProjectiveCurve>::Affine],
                    tables_h: &mut [G],
                    exps_h: &[BigInt],
                    exps_recode_h: &mut [u8],
                    run_parallel: bool,
                ) {
                    let closure = |
                        ((k, exps_chunk), (table, base)):
                        ((&BigInt, &mut [u8]), (&mut [G], &<G as ProjectiveCurve>::Affine))
                    | {
                        let ((k1_neg, mut k1), (k2_neg, mut k2)) = G::glv_scalar_decomposition(*k);
                        let base = base.into_projective();
                        exps_chunk.clone_from_slice(&scalar_recode_glv(&mut k1, &mut k2));

                        table[0] = G::zero();
                        table[TABLE_SIZE / 2] = G::zero();

                        for i in 1..TABLE_SIZE / 2 {
                            let mut res = if k1_neg {
                                table[i - 1] - base
                            } else {
                                table[i - 1] + base
                            };
                            table[i] = res;

                            G::glv_endomorphism_in_place(&mut res.x);
                            table[TABLE_SIZE / 2 + i] =
                                if k2_neg != k1_neg { res.neg() } else { res };
                        }
                    };
                    if run_parallel {
                        exps_h
                            .par_iter()
                            .zip(exps_recode_h.par_chunks_mut(NUM_U8))
                            .zip(tables_h.par_chunks_mut(TABLE_SIZE).zip(bases_h.par_iter()))
                            .for_each(|x| closure(x));
                    } else {
                        exps_h
                            .iter()
                            .zip(exps_recode_h.chunks_mut(NUM_U8))
                            .zip(tables_h.chunks_mut(TABLE_SIZE).zip(bases_h.iter()))
                            .for_each(|x| closure(x));
                    }
                }

                #[kernel_mod(to_mod)]
                #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
                #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
                #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
                pub mod scalar_mul {
                    use algebra::{$curve::$ProjCurve, FpParameters, Zero};
                    use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};

                    const NUM_BITS: isize =
                        <<<$ProjCurve as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
                    const LOG2_W: isize = 5;
                    const TABLE_SIZE: isize = 1 << LOG2_W;
                    const HALF_TABLE_SIZE: isize = 1 << (LOG2_W - 1);
                    const NUM_U8: isize = 2 * ((NUM_BITS - 1) / (2 * (LOG2_W - 1)) + 2);

                    #[kernel_func]
                    pub unsafe fn scalar_mul(
                        #[type_substitute(*const $crate::[<$curve _ $type _scalar_mul_kernel>]::G)]
                        table: *const $ProjCurve,
                        exps: *const u8,
                        #[type_substitute(*mut $crate::[<$curve _ $type _scalar_mul_kernel>]::G)]
                        out: *mut $ProjCurve,
                        n: isize,
                    ) {
                        let i = accel_core::index();
                        if i < n {
                            let mut res = $ProjCurve::zero();

                            res += &(*table.offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8) as isize));
                            res += &(*table.offset(
                                i * TABLE_SIZE + HALF_TABLE_SIZE + *exps.offset(i * NUM_U8 + 1) as isize,
                            ));

                            for j in 1..NUM_U8 as isize / 2 {
                                for _ in 0..(LOG2_W - 1) {
                                    res.double_in_place();
                                }
                                res += &(*table
                                    .offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8 + 2 * j) as isize));
                                res += &(*table.offset(
                                    i * TABLE_SIZE
                                        + HALF_TABLE_SIZE
                                        + *exps.offset(i * NUM_U8 + 2 * j + 1) as isize,
                                ));
                            }
                            *out.offset(i) = res;
                        }
                    }
                }
            }
        }
    }
}
