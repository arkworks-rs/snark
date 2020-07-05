use crate::mnt6_298::fq::Fq;
use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    field_new,
    fields::fp3::{Fp3, Fp3Parameters},
};

pub type Fq3 = Fp3<Fq3Parameters>;

pub struct Fq3Parameters;

impl Fp3Parameters for Fq3Parameters {
    type Fp = Fq;

    #[rustfmt::skip]
    const NONRESIDUE: Fq = field_new!(Fq, BigInteger([
        0x58eefd67fea995ca,
        0x12f14affbb33a004,
        0x4780323da44ac69b,
        0x88acf9bea707eed9,
        0x14bbbb859e8,
    ]));

    const TWO_ADICITY: u32 = 34;

    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: &'static [u64] = &[
        0x69232b75663933bd,
        0xca650efcfc00ee0,
        0x77ca3963fe36f720,
        0xe4cb46632f9bcf7e,
        0xef510453f08f9f30,
        0x9dd5b8fc72f02d83,
        0x7f8d017ed86608ab,
        0xeb2219b3697c97a4,
        0xc8663846ab96996f,
        0x833cd532053eac7d,
        0x1d5b73dfb20bd3cc,
        0x6f5f6da606b59873,
        0x62e990f43dfc42d6,
        0x6878f58,
    ];

    #[rustfmt::skip]
    const QUADRATIC_NONRESIDUE_TO_T: (Fq, Fq, Fq) = (
        field_new!(Fq, BigInteger([
            0x44a4178610a3a4e6,
            0x49321e4d00f35073,
            0xbbc01b9c400c07a1,
            0xd0127c4589095738,
            0x3730de2a45d,
        ])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0])),
    );

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[
        field_new!(Fq, BigInteger([
            0xc3177aefffbb845c,
            0x9b80c702f9961788,
            0xc5df8dcdac70a85a,
            0x29184098647b5197,
            0x1c1223d33c3,
        ])),
        field_new!(Fq, BigInteger([
            0x1c17bb7477085b6a,
            0x2621629c22e83dbb,
            0x21c062106d949dd8,
            0x9d5b981062164ba,
            0x84ad703207,
        ])),
        field_new!(Fq, BigInteger([
            0xdc13fe3f893c203b,
            0x39a7226875df158f,
            0xe34ed98542eefb62,
            0x6f782a843d139e3c,
            0x177280f6ea9,
        ])),
    ];

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[
        field_new!(Fq, BigInteger([
            0xc3177aefffbb845c,
            0x9b80c702f9961788,
            0xc5df8dcdac70a85a,
            0x29184098647b5197,
            0x1c1223d33c3,
        ])),
        field_new!(Fq, BigInteger([
            0xdc13fe3f893c203b,
            0x39a7226875df158f,
            0xe34ed98542eefb62,
            0x6f782a843d139e3c,
            0x177280f6ea9,
        ])),
        field_new!(Fq, BigInteger([
            0x1c17bb7477085b6a,
            0x2621629c22e83dbb,
            0x21c062106d949dd8,
            0x9d5b981062164ba,
            0x84ad703207,
        ])),
    ];
}
