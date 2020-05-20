pub mod g1;
pub mod g2;
pub mod g;
#[cfg(test)]
mod tests;

use algebra_core::curves::bn::{Bn, BnParameters,
    g1::{G1Affine as BnG1Affine, G1Projective as BnG1Projective},
    g2::{G2Affine as BnG2Affine, G2Projective as BnG2Projective},
};
use crate::{
    biginteger::BigInteger384 as BigInteger,
    bn_382::*,
    field_new,
};

pub type Bn_382 = Bn<Bn_382Parameters>;
pub type G1Affine = BnG1Affine<Bn_382Parameters>;
pub type G1Projective = BnG1Projective<Bn_382Parameters>;
pub type G2Affine = BnG2Affine<Bn_382Parameters>;
pub type G2Projective = BnG2Projective<Bn_382Parameters>;

pub struct Bn_382Parameters;

impl BnParameters for Bn_382Parameters {
    const U: &'static [u64] = &[0, 1073873924];
    const SIX_U_PLUS_2_NAF: &'static [i8] = &[
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, -1, 0, 1,
    ];

    const CUBIC_NONRESIDUE_TO_Q_MINUS_1_OVER_2: Fq2 = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger([
                0x16b744a7d72fb912,
                0x8db76da14b98776d,
                0xd7d0fda03758326c,
                0x9a05f3af0ce04699,
                0x1c8a66ecb161efb2,
                0x13a9f1d5f1261bfe
            ])
        ),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0])),
    );
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = self::g1::Bn_382G1Parameters;
    type G2Parameters = self::g2::Bn_382G2Parameters;
}
