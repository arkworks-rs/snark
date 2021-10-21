pub mod g;
pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

use crate::{
    biginteger::BigInteger384 as BigInteger,
    curves::bn::{
        g1::{G1Affine as BnG1Affine, G1Projective as BnG1Projective},
        g2::{G2Affine as BnG2Affine, G2Projective as BnG2Projective},
        Bn, BnParameters, TwistType,
    },
    field_new,
    fields::bn_382::*,
};

pub type Bn382 = Bn<Bn382Parameters>;
pub type G1Affine = BnG1Affine<Bn382Parameters>;
pub type G1Projective = BnG1Projective<Bn382Parameters>;
pub type G2Affine = BnG2Affine<Bn382Parameters>;
pub type G2Projective = BnG2Projective<Bn382Parameters>;

pub struct Bn382Parameters;

impl BnParameters for Bn382Parameters {
    const X: &'static [u64] = &[0, 1073873924];
    const X_IS_NEGATIVE: bool = false;
    const ATE_LOOP_COUNT: &'static [i8] = &[
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, -1, 0, 1,
    ];
    const ATE_LOOP_COUNT_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
    const TWIST_MUL_BY_Q_X: Fq2 = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger([
                0x43ac10f69cd0866e,
                0xb67658d4844670fa,
                0x64500aac20e3e056,
                0xe69857d69abfc002,
                0x521ddf42ec5832c5,
                0xee09eba205fe5d8
            ])
        ),
        field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
    );
    const TWIST_MUL_BY_Q_Y: Fq2 = field_new!(
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
    type G1Parameters = self::g1::Bn382G1Parameters;
    type G2Parameters = self::g2::Bn382G2Parameters;
}
