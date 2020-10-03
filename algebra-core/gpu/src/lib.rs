#[macro_use]
mod cpu_gpu;

#[macro_use]
// We keep this macro module private as the macros should not be used outside of this crate due to dependencies
mod scalar_mul;

// Uncomment to use. Leave commented to reduce compilation overhead
// (This is very significant as we are compiling in sequence n different
// cargo crates for the nvptx target for n different curve impls, with
// very low thread util)

impl_scalar_mul_kernel_glv!(bw6_761, "bw6_761", g1, G1Projective);
// impl_scalar_mul_kernel!(bls12_381, "bls12_381", g1, G1Projective);
// impl_scalar_mul_kernel!(bls12_377, "bls12_377", g1, G1Projective);
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

// #[macro_use]
// mod msm;
// pub use msm::*;
