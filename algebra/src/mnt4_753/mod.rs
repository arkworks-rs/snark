//! This module implements the MNT4_753 curve generated in
//! [[BCTV14]](https://eprint.iacr.org/2014/595). The name denotes that it is a
//! Miyaji--Nakabayashi--Takano curve of embedding degree 4, defined over a 753-bit (prime) field.
//! The main feature of this curve is that its scalar field and base field respectively equal the
//! base field and scalar field of MNT6_753.
//!
//! Curve information:
//! * Base field: q = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB117E776F218059DB80F0DA5CB537E38685ACCE9767254A4638810719AC425F0E39D54522CDD119F5E9063DE245E8001
//! * Scalar field: r = 0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB26C5C28C859A99B3EEBCA9429212636B9DFF97634993AA4D6C381BC3F0057974EA099170FA13A4FD90776E240000001
//! * valuation(q - 1, 2) = 15
//! * valuation(r - 1, 2) = 30
//! * G1 curve equation: y^2 = x^3 + ax + b, where
//!    * a = 2
//!    * b = 0x01373684A8C9DCAE7A016AC5D7748D3313CD8E39051C596560835DF0C9E50A5B59B882A92C78DC537E51A16703EC9855C77FC3D8BB21C8D68BB8CFB9DB4B8C8FBA773111C36C8B1B4E8F1ECE940EF9EAAD265458E06372009C9A0491678EF4
//! * G2 curve equation: y^2 = x^3 + Ax + B, where
//!    * A = Fq2 = (a * NON_RESIDUE, 0)
//!    * B = Fq2(0, b * NON_RESIDUE)
//!    * NON_RESIDUE = 13 is the quadratic non-residue used to construct the extension field Fq2

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
