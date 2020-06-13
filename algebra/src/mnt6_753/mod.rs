//! This module implements the MNT6_753 curve generated in
//! [[BCTV14]](https://eprint.iacr.org/2014/595). The name denotes that it is a
//! Miyaji--Nakabayashi--Takano curve of embedding degree 6, defined over a 753-bit (prime) field.
//! The main feature of this curve is that its scalar field and base field respectively equal the
//! base field and scalar field of MNT4_753.
//!
//! Curve information:
//! * Base field: q = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB26C5C28C859A99B3EEBCA9429212636B9DFF97634993AA4D6C381BC3F0057974EA099170FA13A4FD90776E240000001
//! * Scalar field: r = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB117E776F218059DB80F0DA5CB537E38685ACCE9767254A4638810719AC425F0E39D54522CDD119F5E9063DE245E8001
//! * valuation(q - 1, 2) = 30
//! * valuation(r - 1, 2) = 15
//! * G1 curve equation: y^2 = x^3 + ax + b, where
//!    * a = 11
//!    * b = 0x7DA285E70863C79D56446237CE2E1468D14AE9BB64B2BB01B10E60A5D5DFE0A25714B7985993F62F03B22A9A3C737A1A1E0FCF2C43D7BF847957C34CCA1E3585F9A80A95F401867C4E80F4747FDE5ABA7505BA6FCF2485540B13DFC8468A
//! * G2 curve equation: y^2 = x^3 + Ax + B, where
//!    * A = Fq3(0, 0, a)
//!    * B = Fq3(b * NON_RESIDUE, 0, 0)
//!    * NON_RESIDUE = 11 is the cubic non-residue used to construct the extension field Fq3

#[cfg(feature = "mnt6_753")]
mod curves;
#[cfg(any(feature = "mnt6_753", feature = "mnt4_753", feature = "ed_on_mnt4_753"))]
mod fields;

#[cfg(feature = "mnt6_753")]
pub use curves::*;
#[cfg(any(feature = "mnt6_753", feature = "mnt4_753", feature = "ed_on_mnt4_753"))]
pub use fields::*;
