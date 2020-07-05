use super::*;
use algebra_core::{biginteger::BigInteger256, field_new, fields::*};

pub type Fq6 = Fp6<Fq6Parameters>;

#[derive(Clone, Copy)]
pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp2Params = Fq2Parameters;

    /// NONRESIDUE = U+9
    #[rustfmt::skip]
    const NONRESIDUE: Fq2 = field_new!(Fq2,
       field_new!(Fq, BigInteger256([
            0xf60647ce410d7ff7,
            0x2f3d6f4dd31bd011,
            0x2943337e3940c6d1,
            0x1d9598e8a7e39857,
        ])),
        field_new!(Fq, BigInteger256([
            202099033278250856u64,
            8885205928937022213u64,
            5545221690922665192u64,
            39800542322357402u64,
        ])),
    );

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C1: &'static [Fq2] = &[
        // Fp2::NONRESIDUE^(((q^0) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xd35d438dc58f0d9d,
                0x0a78eb28f5c70b3d,
                0x666ea36f7879462c,
                0xe0a77c19a07df2f,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^(((q^1) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xb5773b104563ab30,
                0x347f91c8a9aa6454,
                0x7a007127242e0991,
                0x1956bcd8118214ec,
            ])),
            field_new!(Fq, BigInteger256([
                0x6e849f1ea0aa4757,
                0xaa1c7b6d89f89141,
                0xb6e713cdfae0ca3a,
                0x26694fbb4e82ebc3,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^2) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x3350c88e13e80b9c,
                0x7dce557cdb5e56b9,
                0x6001b4b8b615564a,
                0x2682e617020217e0,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^(((q^3) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xc9af22f716ad6bad,
                0xb311782a4aa662b2,
                0x19eeaf64e248c7f4,
                0x20273e77e3439f82,
            ])),
            field_new!(Fq, BigInteger256([
                0xacc02860f7ce93ac,
                0x3933d5817ba76b4c,
                0x69e6188b446c8467,
                0xa46036d4417cc55,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^4) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x71930c11d782e155,
                0xa6bb947cffbe3323,
                0xaa303344d4741444,
                0x2c3b3f0d26594943,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^(((q^5) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xf91aba2654e8e3b1,
                0x4771cb2fdc92ce12,
                0xdcb16ae0fc8bdf35,
                0x274aa195cd9d8be4,
            ])),
            field_new!(Fq, BigInteger256([
                0x5cfc50ae18811f8b,
                0x4bb28433cb43988c,
                0x4fd35f13c3b56219,
                0x301949bd2fc8883a,
            ])),
        ),
    ];
    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C2: &'static [Fq2] = &[
        // Fp2::NONRESIDUE^((2*(q^0) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xd35d438dc58f0d9d,
                0x0a78eb28f5c70b3d,
                0x666ea36f7879462c,
                0xe0a77c19a07df2f,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^((2*(q^1) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x7361d77f843abe92,
                0xa5bb2bd3273411fb,
                0x9c941f314b3e2399,
                0x15df9cddbb9fd3ec,
            ])),
            field_new!(Fq, BigInteger256([
                0x5dddfd154bd8c949,
                0x62cb29a5a4445b60,
                0x37bc870a0c7dd2b9,
                0x24830a9d3171f0fd,
            ])),
        ),
        // Fp2::NONRESIDUE^((2*(q^2) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x71930c11d782e155,
                0xa6bb947cffbe3323,
                0xaa303344d4741444,
                0x2c3b3f0d26594943,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^((2*(q^3) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x448a93a57b6762df,
                0xbfd62df528fdeadf,
                0xd858f5d00e9bd47a,
                0x6b03d4d3476ec58,
            ])),
            field_new!(Fq, BigInteger256([
                0x2b19daf4bcc936d1,
                0xa1a54e7a56f4299f,
                0xb533eee05adeaef1,
                0x170c812b84dda0b2,
            ])),
        ),
        // Fp2::NONRESIDUE^((2*(q^4) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x3350c88e13e80b9c,
                0x7dce557cdb5e56b9,
                0x6001b4b8b615564a,
                0x2682e617020217e0,
            ])),
            field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0])),
        ),
        // Fp2::NONRESIDUE^((2*(q^5) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x843420f1d8dadbd6,
                0x31f010c9183fcdb2,
                0x436330b527a76049,
                0x13d47447f11adfe4,
            ])),
            field_new!(Fq, BigInteger256([
                0xef494023a857fa74,
                0x2a925d02d5ab101a,
                0x83b015829ba62f10,
                0x2539111d0c13aea3,
            ])),
        ),
    ];

    #[inline(always)]
    fn mul_fp2_by_nonresidue(fe: &Fq2) -> Fq2 {
        // (c0+u*c1)*(9+u) = (9*c0-c1)+u*(9*c1+c0)
        let mut f = *fe;
        f.double_in_place().double_in_place().double_in_place();
        let c0 = f.c0 + fe.c0 + Fq2Parameters::mul_fp_by_nonresidue(&fe.c1);
        let c1 = f.c1 + fe.c1 + fe.c0;
        field_new!(Fq2, c0, c1)
    }
}
