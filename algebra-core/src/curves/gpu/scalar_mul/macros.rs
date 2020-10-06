#[macro_export]
macro_rules! impl_scalar_mul_kernel {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            use accel::*;

            #[kernel_mod(transparent)]
            #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
            #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
            #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
            pub mod scalar_mul {
                use algebra::{$curve::$ProjCurve};
                use algebra_core::{curves::ProjectiveCurve, fields::PrimeField, FpParameters, Zero};

                const NUM_BITS: isize =
                    <<<$ProjCurve as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
                const LOG2_W: isize = 5;
                const TABLE_SIZE: isize = 1 << LOG2_W;
                const HALF_TABLE_SIZE: isize = 1 << (LOG2_W - 1);
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
            use accel::*;

            #[kernel_mod(transparent)]
            #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
            #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
            #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
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
