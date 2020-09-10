/// This module implements the R1CS equivalent of `algebra::bls12_377`.
#[cfg(feature = "bls12_377")]
pub mod bls12_377;

/// This module implements the R1CS equivalent of `algebra::ed_on_bls12_377`.
#[cfg(feature = "ed_on_bls12_377")]
pub mod ed_on_bls12_377;

/// This module implements the R1CS equivalent of `algebra::ed_on_cp6_782`.
#[cfg(feature = "ed_on_cp6_782")]
pub mod ed_on_cp6_782;

#[cfg(all(not(feature = "ed_on_cp6_782"), feature = "ed_on_bw6_761"))]
pub(crate) mod ed_on_cp6_782;

/// This module implements the R1CS equivalent of `algebra::ed_on_bw6_761`.
#[cfg(feature = "ed_on_bw6_761")]
pub mod ed_on_bw6_761;

/// This module implements the R1CS equivalent of `algebra::ed_on_bn254`.
#[cfg(feature = "ed_on_bn254")]
pub mod ed_on_bn254;

/// This module implements the R1CS equivalent of `algebra::ed_on_bls12_381`.
#[cfg(feature = "ed_on_bls12_381")]
pub mod ed_on_bls12_381;

/// This module implements the R1CS equivalent of `algebra::ed_on_mnt4_298`.
#[cfg(feature = "ed_on_mnt4_298")]
pub mod ed_on_mnt4_298;

/// This module implements the R1CS equivalent of `algebra::ed_on_mnt4_753`.
#[cfg(feature = "ed_on_mnt4_753")]
pub mod ed_on_mnt4_753;

/// This module implements the R1CS equivalent of `algebra::mnt4_298`.
#[cfg(feature = "mnt4_298")]
pub mod mnt4_298;

/// This module implements the R1CS equivalent of `algebra::mnt4_753`.
#[cfg(feature = "mnt4_753")]
pub mod mnt4_753;

/// This module implements the R1CS equivalent of `algebra::mnt6_298`.
#[cfg(feature = "mnt6_298")]
pub mod mnt6_298;

/// This module implements the R1CS equivalent of `algebra::mnt6_753`.
#[cfg(feature = "mnt6_753")]
pub mod mnt6_753;
