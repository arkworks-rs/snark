// Uncomment to use. Leave commented to reduce compilation overhead
// (This is very significant as we are compiling in sequence n different
// cargo crates for n different curve impls, with very low thread util)

// impl_scalar_mul_kernel_glv!(bw6_761, "bw6_761", g1, G1Projective);
// impl_scalar_mul_kernel!(bls12_381, "bls12_381", g1, G1Projective);
impl_scalar_mul_kernel!(bls12_377, "bls12_377", g1, G1Projective);
// impl_scalar_mul_kernel!(bn254, "bn254", g1, G1Projective);
// impl_scalar_mul_kernel!(mnt4_298, "mnt4_298", g1, G1Projective);
// impl_scalar_mul_kernel!(mnt4_753, "mnt4_753", g1, G1Projective);
// impl_scalar_mul_kernel!(mnt6_298, "mnt6_298", g1, G1Projective);
// impl_scalar_mul_kernel!(mnt6_753, "mnt6_753", g1, G1Projective);
//
// impl_scalar_mul_kernel_glv!(bw6_761, "bw6_761", g2, G2Projective);
// impl_scalar_mul_kernel!(bls12_381, "bls12_381", g2, G2Projective);
// impl_scalar_mul_kernel!(bls12_377, "bls12_377", g2, G2Projective);
// impl_scalar_mul_kernel!(bn254, "bn254", g2, G2Projective);
// impl_scalar_mul_kernel!(mnt4_298, "mnt4_298", g2, G2Projective);
// impl_scalar_mul_kernel!(mnt4_753, "mnt4_753", g2, G2Projective);
// impl_scalar_mul_kernel!(mnt6_298, "mnt6_298", g2, G2Projective);
// impl_scalar_mul_kernel!(mnt6_753, "mnt6_753", g2, G2Projective);
//
// impl_scalar_mul_kernel!(ed_on_bw6_761, "ed_on_bw6_761", proj, EdwardsProjective);
// impl_scalar_mul_kernel!(ed_on_bls12_381, "ed_on_bls12_381", proj, EdwardsProjective);
// impl_scalar_mul_kernel!(ed_on_bls12_377, "ed_on_bls12_377", proj, EdwardsProjective);
// impl_scalar_mul_kernel!(ed_on_bn254, "ed_on_bn254", proj, EdwardsProjective);
// impl_scalar_mul_kernel!(ed_on_mnt4_298, "ed_on_mnt4_298", proj, EdwardsProjective);
// impl_scalar_mul_kernel!(ed_on_mnt4_753, "ed_on_mnt4_753", proj, EdwardsProjective);

#[macro_export]
macro_rules! impl_scalar_mul_kernel {
    ($curve: ident, $curve_string:expr, $type: expr, $ProjCurve: ident) => {
        paste::item! {
            pub mod [<$curve _ $type _scalar_mul_kernel>] {
                use accel::*;
                use rayon::prelude::*;

                use algebra::{BigInteger, FpParameters, Zero};
                use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};

                use algebra::$curve::$ProjCurve;

                pub type G1 = $ProjCurve;
                type PrimeF = <G1 as ProjectiveCurve>::ScalarField;
                pub type BigInt = <PrimeF as PrimeField>::BigInt;

                const NUM_BITS: usize = <<PrimeF as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
                const LOG2_W: usize = 5;
                const TABLE_SIZE: usize = 1 << LOG2_W;
                const NUM_U8: usize = (NUM_BITS - 1) / LOG2_W + 1;

                fn scalar_recode(k: &mut BigInt) -> [u8; NUM_U8] {
                    let mut out = [0; NUM_U8];
                    for i in (0..NUM_U8).rev() {
                        out[i] = (k.as_ref()[0] % TABLE_SIZE as u64) as u8;
                        k.divn(LOG2_W as u32);
                    }
                    assert!(k.is_zero());
                    out
                }

                pub fn run_kernel(
                    ctx: &Context,
                    bases_h: &[G1],
                    exps_h: &[BigInt],
                    cuda_group_size: usize,
                ) -> DeviceMemory<G1> {
                    assert_eq!(bases_h.len(), exps_h.len());
                    let n = bases_h.len();
                    let mut tables = DeviceMemory::<G1>::zeros(&ctx, n * TABLE_SIZE);
                    let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
                    let mut out = DeviceMemory::<G1>::zeros(&ctx, n);

                    let now = std::time::Instant::now();

                    exps.par_chunks_mut(NUM_U8)
                        .zip(exps_h.to_vec().par_iter_mut())
                        .for_each(|(exps_chunk, mut k)| {
                            exps_chunk.clone_from_slice(&scalar_recode(&mut k));
                        });

                    println!("Recoded scalars: {}us", now.elapsed().as_micros());
                    println!("{:?}", &exps[..NUM_U8]);

                    let now = std::time::Instant::now();
                    tables
                        .par_chunks_mut(TABLE_SIZE)
                        .zip(bases_h.par_iter())
                        .for_each(|(table, base)| {
                            table[0] = G1::zero();
                            for i in 1..TABLE_SIZE {
                                table[i] = table[i - 1] + base;
                            }
                        });
                    println!("Generated tables: {}us", now.elapsed().as_micros());
                    // Accessible from CPU as usual Rust slice (though this will be slow)
                    // Can this be changed to a memcpy?
                    scalar_mul_kernel::scalar_mul(
                        &ctx,
                        n / cuda_group_size, // grid
                        cuda_group_size,     // block
                        (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr()),
                    )
                    .expect("Kernel call failed");
                    out
                }

                #[kernel_mod]
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
                    #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
                    #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
                    #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
                    pub unsafe fn scalar_mul(
                        table: *const algebra::$curve::$ProjCurve,
                        exps: *const u8,
                        out: *mut algebra::$curve::$ProjCurve,
                    ) {
                        let mut res = $ProjCurve::zero();
                        let i = accel_core::index();

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
            pub mod [<$curve _ $type _scalar_mul_kernel>] {
                use accel::*;
                use rayon::prelude::*;

                use algebra::{BigInteger, FpParameters, Zero};
                use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};
                use std::ops::Neg;

                use algebra::$curve::$ProjCurve;

                pub type G1 = $ProjCurve;
                type PrimeF = <G1 as ProjectiveCurve>::ScalarField;
                pub type BigInt = <PrimeF as PrimeField>::BigInt;

                const NUM_BITS: usize = <<PrimeF as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
                const LOG2_W: usize = 5;
                const TABLE_SIZE: usize = 1 << LOG2_W;
                const NUM_U8: usize = 2 * ((NUM_BITS - 1) / (2 * (LOG2_W - 1)) + 2);

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

                pub fn run_kernel(
                    ctx: &Context,
                    bases_h: &[G1],
                    exps_h: &[BigInt],
                    cuda_group_size: usize,
                ) -> DeviceMemory<G1> {
                    assert_eq!(bases_h.len(), exps_h.len());
                    let n = bases_h.len();
                    let mut tables = DeviceMemory::<G1>::zeros(&ctx, n * TABLE_SIZE);
                    let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
                    let mut out = DeviceMemory::<G1>::zeros(&ctx, n);

                    let now = std::time::Instant::now();
                    let k_vec: Vec<_> = exps_h
                        .iter()
                        .map(|k| G1::glv_scalar_decomposition(*k))
                        .collect();

                    println!("GLV decomp: {}us", now.elapsed().as_micros());

                    let mut k1_scalars: Vec<_> = k_vec.iter().map(|x| (x.0).1).collect();
                    let mut k2_scalars: Vec<_> = k_vec.iter().map(|x| (x.1).1).collect();
                    exps.par_chunks_mut(NUM_U8)
                        .zip(k1_scalars.par_iter_mut().zip(k2_scalars.par_iter_mut()))
                        .for_each(|(exps_chunk, (mut k1, mut k2))| {
                            exps_chunk.clone_from_slice(&scalar_recode_glv(&mut k1, &mut k2));
                        });

                    println!("{:?}", &exps[..NUM_U8]);

                    let now = std::time::Instant::now();
                    let k1_negates: Vec<_> = k_vec.iter().map(|x| (x.0).0).collect();
                    let k2_negates: Vec<_> = k_vec.iter().map(|x| (x.1).0).collect();
                    tables
                        .par_chunks_mut(TABLE_SIZE)
                        .zip(bases_h.par_iter())
                        .zip(k1_negates.par_iter().zip(k2_negates.par_iter()))
                        .for_each(|((table, base), (k1_neg, k2_neg))| {
                            table[0] = G1::zero();
                            table[TABLE_SIZE / 2] = G1::zero();

                            for i in 1..TABLE_SIZE / 2 {
                                let mut res = if *k1_neg {
                                    table[i - 1] - base
                                } else {
                                    table[i - 1] + base
                                };
                                table[i] = res;

                                G1::glv_endomorphism_in_place(&mut res.x);
                                table[TABLE_SIZE / 2 + i] =
                                    if *k2_neg != *k1_neg { res.neg() } else { res };
                            }
                        });
                    println!("Generated tables: {}us", now.elapsed().as_micros());
                    // Accessible from CPU as usual Rust slice (though this will be slow)
                    // Can this be changed to a memcpy?
                    scalar_mul_kernel::scalar_mul(
                        &ctx,
                        n / cuda_group_size, // grid
                        cuda_group_size,     // block
                        (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr()),
                    )
                    .expect("Kernel call failed");
                    out
                }

                #[kernel_mod]
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
                    #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
                    #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
                    #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = [$curve_string]})]
                    pub unsafe fn scalar_mul(
                        table: *const algebra::$curve::$ProjCurve,
                        exps: *const u8,
                        out: *mut algebra::$curve::$ProjCurve,
                    ) {
                        let mut res = $ProjCurve::zero();
                        let i = accel_core::index();

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
