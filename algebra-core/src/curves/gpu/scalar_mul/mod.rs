#[macro_use]
mod kernel_macros;
pub use kernel_macros::*;

#[macro_use]
mod cpu_gpu_macros;

#[macro_use]
mod run_kernel_macros;

use accel::*;
use lazy_static::lazy_static;
use std::sync::Mutex;

use crate::{
    cfg_chunks_mut,
    {
        curves::{AffineCurve, BatchGroupArithmeticSlice},
        fields::PrimeField,
    },
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

lazy_static! {
    pub static ref MICROBENCH_CPU_GPU_AVG_RATIO: Mutex<(Vec<f64>, usize)> = Mutex::new((vec![], 0));
}

// We will use average of the proportions of throughput (points/s)
// Preferably, one could make this mangled and curve specific.
pub trait GPUScalarMul<G: AffineCurve>: Sized {
    const NUM_BITS: usize;
    const LOG2_W: usize;

    fn table_size() -> usize {
        1 << Self::LOG2_W
    }

    fn num_u8() -> usize;

    fn par_run_kernel(
        ctx: &Context,
        bases_h: &[G],
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
    ) -> DeviceMemory<Self>;

    fn par_run_kernel_sync<T>(
        ctx: &Context,
        bases_h: &[G],
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
        lock: T,
    ) -> DeviceMemory<Self>;

    fn generate_tables_and_recoding(
        bases_h: &[G],
        tables_h: &mut [Self],
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        exps_recode_h: &mut [u8],
    );

    fn cpu_gpu_load_balance_run_kernel(
        ctx: &Context,
        bases_h: &[G],
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
        // size of a single job in the queue e.g. 2 << 14
        job_size: usize,
        // size of the batch for cpu scalar mul
        cpu_chunk_size: usize,
    ) -> Vec<G>;

    fn cpu_gpu_static_partition_run_kernel(
        bases_h: &mut [G],
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
        // size of the batch for cpu scalar mul
        cpu_chunk_size: usize,
    );
}

#[macro_export]
macro_rules! impl_gpu_sw_projective {
    ($Parameters:ident) => {
        impl<P: $Parameters> GPUScalarMul<GroupAffine<P>> for GroupProjective<P> {
            const NUM_BITS: usize =
                <<<Self as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
            const LOG2_W: usize = 5;

            fn num_u8() -> usize {
                if P::has_glv() {
                    2 * ((Self::NUM_BITS - 1) / (2 * (Self::LOG2_W - 1)) + 2)
                } else {
                    (Self::NUM_BITS - 1) / Self::LOG2_W + 1
                }
            }

            fn generate_tables_and_recoding(
                bases_h: &[<Self as ProjectiveCurve>::Affine],
                tables_h: &mut [Self],
                exps_h: &[<<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt],
                exps_recode_h: &mut [u8],
            ) {
                if P::has_glv() {
                    let scalar_recode_glv =
                        |k1: &mut <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt, k2: &mut <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt| -> Vec<u8> {
                            let table_size_glv: u64 = 1u64 << (Self::LOG2_W - 1);
                            let mut out = vec![0; Self::num_u8()];
                            for i in (0..Self::num_u8() / 2).rev() {
                                out[2 * i] = (k1.as_ref()[0] % table_size_glv) as u8;
                                out[2 * i + 1] = (k2.as_ref()[0] % table_size_glv) as u8;
                                k1.divn(Self::LOG2_W as u32 - 1);
                                k2.divn(Self::LOG2_W as u32 - 1);
                            }
                            assert!(k1.is_zero());
                            assert!(k2.is_zero());
                            out
                        };

                    cfg_iter!(exps_h)
                        .zip(cfg_chunks_mut!(exps_recode_h, Self::num_u8()))
                        .zip(cfg_chunks_mut!(tables_h, Self::table_size()).zip(cfg_iter!(bases_h)))
                        .for_each(|((k, exps_chunk), (table, base))| {
                            let ((k1_neg, mut k1), (k2_neg, mut k2)) =
                                P::glv_scalar_decomposition(*k);
                            let base = base.into_projective();
                            exps_chunk.clone_from_slice(&scalar_recode_glv(&mut k1, &mut k2));

                            table[0] = Self::zero();
                            table[Self::table_size() / 2] = Self::zero();

                            for i in 1..Self::table_size() / 2 {
                                let mut res = if k1_neg {
                                    table[i - 1] - base
                                } else {
                                    table[i - 1] + base
                                };
                                table[i] = res;

                                P::glv_endomorphism_in_place(&mut res.x);
                                table[Self::table_size() / 2 + i] =
                                    if k2_neg != k1_neg { res.neg() } else { res };
                            }
                        });
                } else {
                    let scalar_recode = |k: &mut <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt| -> Vec<u8> {
                        let mut out = vec![0; Self::num_u8()];
                        for i in (0..Self::num_u8()).rev() {
                            out[i] = (k.as_ref()[0] % Self::table_size() as u64) as u8;
                            k.divn(Self::LOG2_W as u32);
                        }
                        assert!(k.is_zero());
                        out
                    };
                    cfg_iter!(exps_h)
                        .zip(cfg_chunks_mut!(exps_recode_h, Self::num_u8()))
                        .zip(cfg_chunks_mut!(tables_h, Self::table_size()).zip(cfg_iter!(bases_h)))
                        .for_each(|((k, exps_chunk), (table, base))| {
                            let base = base.into_projective();
                            exps_chunk.clone_from_slice(&scalar_recode(&mut k.clone())[..]);

                            table[0] = Self::zero();
                            for i in 1..Self::table_size() {
                                table[i] = table[i - 1] + base;
                            }
                        });
                }
            }

            impl_run_kernel!();
            impl_gpu_cpu_run_kernel!();
        }
    };
}

#[macro_export]
macro_rules! impl_gpu_te_projective {
    ($Parameters:ident) => {
        impl<P: $Parameters> GPUScalarMul<GroupAffine<P>> for GroupProjective<P> {
            const NUM_BITS: usize =
                <<<Self as ProjectiveCurve>::ScalarField as PrimeField>::Params as FpParameters>::MODULUS_BITS as usize;
            const LOG2_W: usize = 5;

            fn generate_tables_and_recoding(
                bases_h: &[<Self as ProjectiveCurve>::Affine],
                tables_h: &mut [Self],
                exps_h: &[<<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt],
                exps_recode_h: &mut [u8],
            ) {
                let scalar_recode = |k: &mut <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt| -> Vec<u8> {
                    let mut out = vec![0; Self::num_u8()];
                    for i in (0..Self::num_u8()).rev() {
                        out[i] = (k.as_ref()[0] % Self::table_size() as u64) as u8;
                        k.divn(Self::LOG2_W as u32);
                    }
                    assert!(k.is_zero());
                    out
                };
                cfg_iter!(exps_h)
                    .zip(cfg_chunks_mut!(exps_recode_h, Self::num_u8()))
                    .zip(cfg_chunks_mut!(tables_h, Self::table_size()).zip(cfg_iter!(bases_h)))
                    .for_each(|((k, exps_chunk), (table, base))| {
                        let base = base.into_projective();
                        exps_chunk.clone_from_slice(&scalar_recode(&mut k.clone())[..]);

                        table[0] = Self::zero();
                        for i in 1..Self::table_size() {
                            table[i] = table[i - 1] + base;
                        }
                    });
            }

            fn num_u8() -> usize {
                (Self::NUM_BITS - 1) / Self::LOG2_W + 1
            }

            impl_run_kernel!();
            impl_gpu_cpu_run_kernel!();
        }
    };
}

pub trait GPUScalarMulSlice<G: AffineCurve> {
    fn cpu_gpu_scalar_mul(
        &mut self,
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
        // size of the batch for cpu scalar mul
        cpu_chunk_size: usize,
    );
}

impl<G: AffineCurve> GPUScalarMulSlice<G> for [G] {
    fn cpu_gpu_scalar_mul(
        &mut self,
        exps_h: &[<<G as AffineCurve>::ScalarField as PrimeField>::BigInt],
        cuda_group_size: usize,
        // size of the batch for cpu scalar mul
        cpu_chunk_size: usize,
    ) {
        if accel::Device::init() && cfg!(feature = "gpu") {
            <G as AffineCurve>::Projective::cpu_gpu_static_partition_run_kernel(
                self,
                exps_h,
                cuda_group_size,
                cpu_chunk_size,
            );
        } else {
            let mut exps_mut = exps_h.to_vec();
            cfg_chunks_mut!(self, cpu_chunk_size)
                .zip(cfg_chunks_mut!(exps_mut, cpu_chunk_size))
                .for_each(|(b, s)| {
                    b[..].batch_scalar_mul_in_place(&mut s[..], 4);
                }
            );
        }
    }
}
