//! This module implements a twisted Edwards curve whose base field is the scalar field of the
//! curve BLS12-377.  This allows defining cryptographic primitives that use elliptic curves over
//! the scalar field of the latter curve. This curve was generated as part of the paper
//! [[BCGMMW20, “Zexe”]](https://eprint.iacr.org/2018/962).
//!
//! Curve information:
//! * Base field: q = 8444461749428370424248824938781546531375899335154063827935233455917409239041
//! * Scalar field: r = 2111115437357092606062206234695386632838870926408408195193685246394721360383
//! * Valuation(q - 1, 2) = 47
//! * Valuation(r - 1, 2) = 1
//! * Curve equation: ax^2 + y^2 =1 + dx^2y^2, where
//!    * a = -1
//!    * d = 3021

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
