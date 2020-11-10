#[macro_export]
macro_rules! impl_scalar_mul_kernel {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            #[cfg(feature = "cuda")]
            use {accel::*, std::sync::{Arc, Mutex}};

            #[cfg(not(feature = "cuda"))]
            use algebra_core::accel_dummy::*;

            use algebra_core::curves::cuda::scalar_mul::ScalarMulProfiler;

            #[cfg(feature = "cuda")]
            lazy_static::lazy_static! {
                pub static ref MICROBENCH_CPU_GPU_AVG_RATIO:
                    Arc<Mutex<(Vec<f64>, usize)>> = Arc::new(Mutex::new((vec![], 0)));
            }

            #[cfg(not(feature = "cuda"))]
            static MICROBENCH_CPU_GPU_AVG_RATIO: () = ();

            const NAMESPACE: &'static str = stringify!([<$curve _ $type _cuda_namespace>]);

            #[cfg(feature = "cuda")]
            #[kernel_mod(transparent)]
            #[name([<$curve _ $type _cuda_namespace>])]
            #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
            #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "master", package = "algebra-core", default_features = false})]
            #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "master", package = "algebra", default_features = false, features = [$curve_string]})]
            pub mod scalar_mul {
                use algebra::{$curve::$ProjCurve};
                use algebra_core::{curves::ProjectiveCurve, fields::PrimeField, FpParameters, Zero};

                const NUM_BITS: isize =
                    <<<$ProjCurve as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
                const LOG2_W: isize = 5;
                const TABLE_SIZE: isize = 1 << LOG2_W;
                const NUM_U8: isize = (NUM_BITS - 1) / LOG2_W + 1;

                #[kernel_func]
                pub unsafe fn scalar_mul(
                    #[type_substitute(*const super::$ProjCurve)]
                    table: *const $ProjCurve,
                    exps: *const u8,
                    #[type_substitute(*mut super::$ProjCurve)]
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

#[macro_export]
macro_rules! impl_scalar_mul_kernel_glv {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            #[cfg(feature = "cuda")]
            use {accel::*, std::sync::{Arc, Mutex}};

            #[cfg(not(feature = "cuda"))]
            use algebra_core::accel_dummy::*;

            use algebra_core::curves::cuda::scalar_mul::ScalarMulProfiler;

            #[cfg(feature = "cuda")]
            lazy_static::lazy_static! {
                pub static ref MICROBENCH_CPU_GPU_AVG_RATIO:
                    Arc<Mutex<(Vec<f64>, usize)>> = Arc::new(Mutex::new((vec![], 0)));
            }

            #[cfg(not(feature = "cuda"))]
            static MICROBENCH_CPU_GPU_AVG_RATIO: () = ();

            const NAMESPACE: &'static str = stringify!([<$curve _ $type _cuda_namespace>]);

            #[cfg(feature = "cuda")]
            #[kernel_mod(transparent)]
            #[name([<$curve _ $type _cuda_namespace>])]
            #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
            #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "master", package = "algebra-core", default_features = false})]
            #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "master", package = "algebra", default_features = false, features = [$curve_string]})]
            pub mod scalar_mul {
                use algebra::{$curve::$ProjCurve};
                use algebra_core::{curves::ProjectiveCurve, fields::PrimeField, FpParameters, Zero};

                const NUM_BITS: isize =
                    <<<$ProjCurve as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
                const LOG2_W: isize = 5;
                const TABLE_SIZE: isize = 1 << LOG2_W;
                const HALF_TABLE_SIZE: isize = 1 << (LOG2_W - 1);
                const NUM_U8: isize = 2 * ((NUM_BITS - 1) / (2 * (LOG2_W - 1)) + 2);

                #[kernel_func]
                pub unsafe fn scalar_mul(
                    #[type_substitute(*const super::$ProjCurve)]
                    table: *const $ProjCurve,
                    exps: *const u8,
                    #[type_substitute(*mut super::$ProjCurve)]
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

#[macro_export]
macro_rules! impl_scalar_mul_parameters {
    ($ProjCurve:ident) => {
        #[allow(unused_variables)]
        fn scalar_mul_kernel(
            ctx: &Context,
            grid: usize,
            block: usize,
            table: *const $ProjCurve,
            exps: *const u8,
            out: *mut $ProjCurve,
            n: isize,
        ) -> error::Result<()> {
            #[cfg(feature = "cuda")]
            scalar_mul(ctx, grid, block, (table, exps, out, n))
        }

        fn scalar_mul_static_profiler() -> ScalarMulProfiler {
            #[cfg(feature = "cuda")]
            return (*MICROBENCH_CPU_GPU_AVG_RATIO).clone();

            #[cfg(not(feature = "cuda"))]
            MICROBENCH_CPU_GPU_AVG_RATIO
        }

        fn namespace() -> &'static str {
            NAMESPACE
        }
    };
}
