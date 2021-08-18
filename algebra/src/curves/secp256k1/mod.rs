//! The SECG curve secp256k1 from [SECG-SEC2-v2](https://www.secg.org/sec2-v2.pdf),
//! a prime order curve at a security level of 128 bit.
use crate::{Field, field_new};
use crate::biginteger::BigInteger320;
use crate::fields::secp256k1::{fq::Fq, fr::Fr};
use crate::curves::{
    models::{ModelParameters, SWModelParameters},
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
};

#[cfg(test)]
mod tests;

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Secp256k1Parameters;

impl ModelParameters for Secp256k1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

pub type Affine = GroupAffine<Secp256k1Parameters>;
pub type Projective = GroupProjective<Secp256k1Parameters>;

impl SWModelParameters for Secp256k1Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger320([0x0, 0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 7
    const COEFF_B: Fq = field_new!(Fq, BigInteger320([
        0x0000000000000000, 
        0x0000000700001ab7, 
        0x0000000000000000, 
        0x0000000000000000, 
        0x0
    ]));

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = 1
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger320([
        0x0000000000000000,
        0x402da1732fc9bebf,
        0x4551231950b75fc4,
        0x0000000000000001,
        0x0,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G_GENERATOR_X, G_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

/// G_GENERATOR_X = 
/// = 55066263022277343669578718895168534326250603453777594175500187360389116729240
pub const G_GENERATOR_X: Fq = field_new!(Fq, BigInteger320([
    0xc1c8687459e7e1c8,
    0xd7362e5ae2000924,
    0x231e295329bc66db,
    0x979f48c033fd129c,
    0x0
]));

/// G_GENERATOR_Y =
/// = 32670510020758816978083085130507043184471273380659243275938904335757337482424
pub const G_GENERATOR_Y: Fq = field_new!(Fq, BigInteger320([
    0xc61091508ba852b6,
    0xb15ea6d3a31b3418,
    0x8dfc5d5d1f1dc64d,
    0x70b6b59aac19c136,
    0x0
]));