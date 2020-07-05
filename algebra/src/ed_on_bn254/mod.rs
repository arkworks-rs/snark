//! This module implements a twisted Edwards curve whose base field is the scalar field of the
//! curve BN254. This allows defining cryptographic primitives that use elliptic curves over
//! the scalar field of the latter curve. This curve is also known as [Baby-Jubjub](https://github.com/barryWhiteHat/baby_jubjub).
//!
//! Curve information:
//! * Base field: q = 21888242871839275222246405745257275088548364400416034343698204186575808495617
//! * Scalar field: r = 2736030358979909402780800718157159386076813972158567259200215660948447373041
//! * Valuation(q - 1, 2) = 28
//! * Valuation(r - 1, 2) = 4
//! * Curve equation: ax^2 + y^2 =1 + dx^2y^2, where
//!    * a = 1
//!    * d = 168696/168700 mod q
//!        = 9706598848417545097372247223557719406784115219466060233080913168975159366771

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
