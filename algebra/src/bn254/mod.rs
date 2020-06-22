//! This module implements the BN254 curve that was sampled as part of the [[BCTV14]][https://eprint.iacr.org/2013/879.pdf] paper .
//! The name denotes that it is a Barreto--Naehrig curve of embedding degree 12,
//! defined over a 254-bit (prime) field. The scalar field is highly 2-adic.
//!
//! This curve is also implemented in [libff](https://github.com/scipr-lab/libff/tree/master/libff/algebra/curves/alt_bn128) under the name `bn128`.
//! It is the same as the `bn256` curve used in Ethereum (eg: [go-ethereum](https://github.com/ethereum/go-ethereum/tree/master/crypto/bn254/cloudflare)).
//!
//! #CAUTION
//! **This curve does not satisfy the 128-bit security level anymore.**
//!
//!
//! Curve information:
//! * Base field: q = 21888242871839275222246405745257275088696311157297823662689037894645226208583
//! * Scalar field: r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
//! * valuation(q - 1, 2) = 1
//! * valuation(r - 1, 2) = 28
//! * G1 curve equation: y^2 = x^3 + 3
//! * G2 curve equation: y^2 = x^3 + B, where
//!    * B = 3/(u+9) where Fq2[u]=Fq/u+1
//!        = Fq2(19485874751759354771024239261021720505790618469301721065564631296452457478373, 266929791119991161246907387137283842545076965332900288569378510910307636690)

#[cfg(feature = "bn254")]
mod curves;

mod fields;

#[cfg(feature = "bn254")]
pub use curves::*;

pub use fields::*;
