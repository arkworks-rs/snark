use accel::*;
use rayon::prelude::*;
use std::sync::Mutex;
use lazy_static::lazy_static;

use algebra_core::{
    biginteger::BigInteger, FpParameters, Zero,
    curves::{ProjectiveCurve, AffineCurve, BatchGroupArithmeticSlice},
    fields::PrimeField,
};

pub trait GPUScalarMul {

}

// This ought to be instantiated concretely
pub trait GPUParameters<G: AffineCurve> {
    type AffineGroup = G;
    // This is to be instantiated with macro
    fn scalar_mul_kernel();
}

impl<G: AffineCurve, P: GPUParameters> GPUScalarMul<P> for G {
    type PrimeF = <Self::Projective as ProjectiveCurve>::ScalarField;
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
        bases_h: &[Self],
        tables_h: &mut [<Self as AffineCurve>::Projective],
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
}

pub trait GPUScalarMulSlice {

}

impl<G: AffineCurve> GPUScalarMulSlice for [G] {

}

impl<G: ProjectiveCurve> GPUScalarMulSlice for [G] {

}
