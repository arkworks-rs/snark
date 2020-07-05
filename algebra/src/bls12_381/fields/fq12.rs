use crate::bls12_381::*;
use algebra_core::{biginteger::BigInteger384, field_new, fields::*};

pub type Fq12 = Fp12<Fq12Parameters>;

#[derive(Clone, Copy)]
pub struct Fq12Parameters;

impl Fp12Parameters for Fq12Parameters {
    type Fp6Params = Fq6Parameters;

    const NONRESIDUE: Fq6 = field_new!(Fq6, FQ2_ZERO, FQ2_ONE, FQ2_ZERO);

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP12_C1: &'static [Fq2] = &[
        // Fq2(u + 1)**(((q^0) - 1) / 6)
        FQ2_ONE,
        // Fq2(u + 1)**(((q^1) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x7089552b319d465,
                0xc6695f92b50a8313,
                0x97e83cccd117228f,
                0xa35baecab2dc29ee,
                0x1ce393ea5daace4d,
                0x8f2220fb0fb66eb,
            ])),
            field_new!(Fq, BigInteger384([
                0xb2f66aad4ce5d646,
                0x5842a06bfc497cec,
                0xcf4895d42599d394,
                0xc11b9cba40a8e8d0,
                0x2e3813cbe5a0de89,
                0x110eefda88847faf,
            ])),
        ),
        // Fq2(u + 1)**(((q^2) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xecfb361b798dba3a,
                0xc100ddb891865a2c,
                0xec08ff1232bda8e,
                0xd5c13cc6f1ca4721,
                0x47222a47bf7b5c04,
                0x110f184e51c5f59,
            ])),
            FQ_ZERO,
        ),
        // Fq2(u + 1)**(((q^3) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
                0x1d794e4fac7cf0b9,
                0xbd592fc7d825ec8,
            ])),
            field_new!(Fq, BigInteger384([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0xe2b7eedbbfd87d2,
            ])),
        ),
        // Fq2(u + 1)**(((q^4) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x30f1361b798a64e8,
                0xf3b8ddab7ece5a2a,
                0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b,
                0x3636b76660701c6e,
                0x51ba4ab241b6160,
            ])),
            FQ_ZERO,
        ),
        // Fq2(u + 1)**(((q^5) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x3726c30af242c66c,
                0x7c2ac1aad1b6fe70,
                0xa04007fbba4b14a2,
                0xef517c3266341429,
                0x95ba654ed2226b,
                0x2e370eccc86f7dd,
            ])),
            field_new!(Fq, BigInteger384([
                0x82d83cf50dbce43f,
                0xa2813e53df9d018f,
                0xc6f0caa53c65e181,
                0x7525cf528d50fe95,
                0x4a85ed50f4798a6b,
                0x171da0fd6cf8eebd,
            ])),
        ),
        // Fq2(u + 1)**(((q^6) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x43f5fffffffcaaae,
                0x32b7fff2ed47fffd,
                0x7e83a49a2e99d69,
                0xeca8f3318332bb7a,
                0xef148d1ea0f4c069,
                0x40ab3263eff0206,
            ])),
            FQ_ZERO,
        ),
        // Fq2(u + 1)**(((q^7) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xb2f66aad4ce5d646,
                0x5842a06bfc497cec,
                0xcf4895d42599d394,
                0xc11b9cba40a8e8d0,
                0x2e3813cbe5a0de89,
                0x110eefda88847faf,
            ])),
            field_new!(Fq, BigInteger384([
                0x7089552b319d465,
                0xc6695f92b50a8313,
                0x97e83cccd117228f,
                0xa35baecab2dc29ee,
                0x1ce393ea5daace4d,
                0x8f2220fb0fb66eb,
            ])),
        ),
        // Fq2(u + 1)**(((q^8) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ])),
            FQ_ZERO,
        ),
        // Fq2(u + 1)**(((q^9) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0xe2b7eedbbfd87d2,
            ])),
            field_new!(Fq, BigInteger384([
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
                0x1d794e4fac7cf0b9,
                0xbd592fc7d825ec8,
            ])),
        ),
        // Fq2(u + 1)**(((q^10) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
                0x14e4f04fe2db9068,
                0x14e56d3f1564853a,
            ])),
            FQ_ZERO,
        ),
        // Fq2(u + 1)**(((q^11) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x82d83cf50dbce43f,
                0xa2813e53df9d018f,
                0xc6f0caa53c65e181,
                0x7525cf528d50fe95,
                0x4a85ed50f4798a6b,
                0x171da0fd6cf8eebd,
            ])),
            field_new!(Fq, BigInteger384([
                0x3726c30af242c66c,
                0x7c2ac1aad1b6fe70,
                0xa04007fbba4b14a2,
                0xef517c3266341429,
                0x95ba654ed2226b,
                0x2e370eccc86f7dd,
            ])),
        ),
    ];
}
