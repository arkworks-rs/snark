//! This module implements the MNT6_298 curve generated in
//! [[BCTV14]](https://eprint.iacr.org/2014/595).  The name denotes that it is a
//! Miyaji--Nakabayashi--Takano curve of embedding degree 6, defined over a 298-bit (prime) field.
//! The main feature of this curve is that its scalar field and base field respectively equal the
//! base field and scalar field of MNT4_298.
//!
//!
//! Curve information:
//! * Scalar field: q = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
//! * Base field: r = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
//! * valuation(q - 1, 2) = 34
//! * valuation(r - 1, 2) = 17
//! * G1 curve equation: y^2 = x^3 + ax + b, where
//!    * a = 11
//!    * b = 106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074
//! * G2 curve equation: y^2 = x^3 + Ax + B, where
//!    * A = Fq2 = (0, 0, a)
//!    * B = Fq2(b * NON_RESIDUE, 0, 0)
//!    * NON_RESIDUE = 5 is the cubic non-residue used to construct the field extension Fq3

#[cfg(feature = "mnt6_298")]
mod curves;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
mod fields;

#[cfg(feature = "mnt6_298")]
pub use curves::*;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
pub use fields::*;
