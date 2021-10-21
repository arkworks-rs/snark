use crate::{
    biginteger::BigInteger384,
    field_new,
    fields::bn_382::{
        fq2::{Fq2, FQ2_ONE, FQ2_ZERO},
        Fq, Fq6, Fq6Parameters,
    },
    fp12_2over3over2::{Fp12, Fp12Parameters},
};

pub type Fq12 = Fp12<Fq12Parameters>;

#[derive(Clone, Copy)]
pub struct Fq12Parameters;

impl Fp12Parameters for Fq12Parameters {
    type Fp6Params = Fq6Parameters;

    const NONRESIDUE: Fq6 = field_new!(Fq6, FQ2_ZERO, FQ2_ONE, FQ2_ZERO);

    const FROBENIUS_COEFF_FP12_C1: &'static [Fq2] = &[
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xfffffffffffffff9,
                    0xfffffff57fab5757,
                    0x7f56ac056aeaf57f,
                    0x10388572e3c2c0f5,
                    0xe6ce591c2bafc343,
                    0x3e03f4104144b1a
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x51b19a7e3871df8a,
                    0xff256c8a6096ca14,
                    0x3c5ed207a2e9ac81,
                    0xee047eb105d3e89c,
                    0x59e5bf1f71597093,
                    0x2226c77500bb1b4b
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x43ac10f69cd0866e,
                    0xb67658d4844670fa,
                    0x64500aac20e3e056,
                    0xe69857d69abfc002,
                    0x521ddf42ec5832c5,
                    0xee09eba205fe5d8
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x16b744a7d72fb912,
                    0x8db76da14b98776d,
                    0xd7d0fda03758326c,
                    0x9a05f3af0ce04699,
                    0x1c8a66ecb161efb2,
                    0x13a9f1d5f1261bfe
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x43ac10f69cd08675,
                    0xb67658df049b19a2,
                    0xe4f95ea6b5f8ead6,
                    0xd65fd263b6fcff0c,
                    0x6b4f8626c0a86f82,
                    0xb005f791c4b9abd
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xc505aa299ebdd989,
                    0x8e9201186b0dc570,
                    0x1b8a5c2a1771876a,
                    0x608bab122fa766ff,
                    0x33f4e436f0a63ea7,
                    0x1587b3a0cb438841
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x8,
                    0xc0060c0c0,
                    0xc1848c18180c00,
                    0xa451b0a144d8480c,
                    0x8a81e34d84edfc45,
                    0x202449fed6c43c73
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xae4e6581c78e2077,
                    0xda93771f754e03,
                    0x43b95e89e01954fe,
                    0xc685b76322c72065,
                    0x176a7d4a3f444ef4,
                    0x1ddc1cada1d6c43
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xbc53ef09632f7993,
                    0x4989a72cfbc5a71d,
                    0x1bc825e5621f2129,
                    0xcdf1de3d8ddb48ff,
                    0x1f325d26c4458cc2,
                    0x1523ea85ba78a1b6
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xe948bb5828d046ef,
                    0x724892603473a0aa,
                    0xa84732f14baacf13,
                    0x1a8442651bbac267,
                    0x54c5d57cff3bcfd6,
                    0x105a9769e9b26b90
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0xbc53ef09632f798c,
                    0x4989a7227b70fe75,
                    0x9b1ed1eacd0a16a9,
                    0xde2a63b0719e09f4,
                    0x600b642eff55005,
                    0x190429c6be8cecd1
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger384([
                    0x3afa55d661422678,
                    0x716dfee914fe52a7,
                    0x648dd4676b917a15,
                    0x53fe8b01f8f3a202,
                    0x3d5b5832bff780e1,
                    0xe7cd59f0f94ff4d
                ])
            ),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
    ];
}
