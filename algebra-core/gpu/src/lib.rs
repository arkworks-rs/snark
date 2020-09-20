use accel::*;
use rayon::prelude::*;

use algebra::{bw6_761::G1Projective, BigInteger, FpParameters, Zero};
use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};
use std::ops::Neg;

pub type G1 = G1Projective;
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

pub fn gpu_scalar_mul(
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

    if G1::has_glv() {
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
                    table[TABLE_SIZE / 2 + i] = if *k2_neg != *k1_neg { res.neg() } else { res };
                }
            });
        println!("Generated tables: {}us", now.elapsed().as_micros());
        // Accessible from CPU as usual Rust slice (though this will be slow)
        // Can this be changed to a memcpy?
        kernel::scalar_mul(
            &ctx,
            n / cuda_group_size, // grid
            cuda_group_size,     // block
            (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr()),
        )
        .expect("Kernel call failed");
    } else {
        ()
    }
    out
}

mod kernel {
    #![allow(unused)]
    use accel::*;
    #[kernel_mod]
    pub mod scalar_mul {
        use algebra::{bw6_761::G1Projective, FpParameters, Zero};
        use algebra_core::{curves::ProjectiveCurve, fields::PrimeField};

        pub type G1 = G1Projective;
        type PrimeF = <G1 as ProjectiveCurve>::ScalarField;
        pub type BigInt = <PrimeF as PrimeField>::BigInt;

        const NUM_BITS: isize = <<PrimeF as PrimeField>::Params as FpParameters>::MODULUS_BITS as isize;
        const LOG2_W: isize = 5;
        const TABLE_SIZE: isize = 1 << LOG2_W;
        const HALF_TABLE_SIZE: isize = 1 << (LOG2_W - 1);
        const NUM_U8: isize = 2 * ((NUM_BITS - 1) / (2 * (LOG2_W - 1)) + 2);

        #[kernel_func]
        #[dependencies("accel-core" = { git = "https://github.com/jon-chuang/accel", package = "accel-core" })]
        #[dependencies("algebra-core" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra-core", default_features = false})]
        #[dependencies("algebra" = { git = "https://github.com/celo-org/zexe", branch = "jonch/gpu_sc_mul", package = "algebra", default_features = false, features = ["bw6_761"]})]
        pub unsafe fn scalar_mul(table: *const crate::G1, exps: *const u8, out: *mut crate::G1) {
            if G1::has_glv() {
                let mut res = G1::zero();
                let i = accel_core::index();

                res += &(*table.offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8) as isize));
                res += &(*table
                    .offset(i * TABLE_SIZE + HALF_TABLE_SIZE + *exps.offset(i * NUM_U8 + 1) as isize));

                for j in 1..NUM_U8 as isize / 2 {
                    for _ in 0..(LOG2_W - 1) {
                        res.double_in_place();
                    }
                    res += &(*table.offset(i * TABLE_SIZE + *exps.offset(i * NUM_U8 + 2 * j) as isize));
                    res += &(*table.offset(
                        i * TABLE_SIZE + HALF_TABLE_SIZE + *exps.offset(i * NUM_U8 + 2 * j + 1) as isize,
                    ));
                }
                *out.offset(i) = res;
            } else {
                ()
            }
        }
    }
}
