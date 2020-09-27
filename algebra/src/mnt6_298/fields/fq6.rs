use crate::mnt6_298::{Fq, Fq3, Fq3Parameters, FQ_ONE, FQ_ZERO};
use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    field_new,
    fields::fp6_2over3::{Fp6, Fp6Parameters},
};

pub type Fq6 = Fp6<Fq6Parameters>;

pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp3Params = Fq3Parameters;

    #[rustfmt::skip]
    const NONRESIDUE: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C1: &'static [Fq] = &[
        field_new!(Fq, BigInteger([
            0xc3177aefffbb845c,
            0x9b80c702f9961788,
            0xc5df8dcdac70a85a,
            0x29184098647b5197,
            0x1c1223d33c3,
        ])),
        field_new!(Fq, BigInteger([
            0xdf2f366476c3dfc6,
            0xc1a2299f1c7e5543,
            0xe79fefde1a054632,
            0x32edfa196a9cb651,
            0x245cfad65ca,
        ])),
        field_new!(Fq, BigInteger([
            0x1c17bb7477085b6a,
            0x2621629c22e83dbb,
            0x21c062106d949dd8,
            0x9d5b981062164ba,
            0x84ad703207,
        ])),
        field_new!(Fq, BigInteger([
            0xf82bb9b400447ba5,
            0x5fc8850498c7534a,
            0x50f3b95b083993a,
            0x794de405433502f7,
            0x1fbd57fa0b0,
        ])),
        field_new!(Fq, BigInteger([
            0xdc13fe3f893c203b,
            0x39a7226875df158f,
            0xe34ed98542eefb62,
            0x6f782a843d139e3c,
            0x177280f6ea9,
        ])),
        field_new!(Fq, BigInteger([
            0x9f2b792f88f7a497,
            0xd527e96b6f752d18,
            0xa92e6752ef5fa3bc,
            0x98906b1ca18eefd4,
            0x3384a4ca26c,
        ])),
    ];
}
