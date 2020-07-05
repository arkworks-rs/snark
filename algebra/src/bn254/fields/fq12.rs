use super::*;
use algebra_core::{biginteger::BigInteger256, field_new, fields::*};

pub type Fq12 = Fp12<Fq12Parameters>;

#[derive(Clone, Copy)]
pub struct Fq12Parameters;

impl Fp12Parameters for Fq12Parameters {
    type Fp6Params = Fq6Parameters;

    const NONRESIDUE: Fq6 = field_new!(Fq6, FQ2_ZERO, FQ2_ONE, FQ2_ZERO);

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP12_C1: &'static [Fq2] = &[
        // Fp2::NONRESIDUE^(((q^0) - 1) / 6)
        FQ2_ONE,
        // Fp2::NONRESIDUE^(((q^1) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xaf9ba69633144907,
                0xca6b1d7387afb78a,
                0x11bded5ef08a2087,
                0x02f34d751a1f3a7c,
            ])),
            field_new!(Fq, BigInteger256([
                0xa222ae234c492d72,
                0xd00f02a4565de15b,
                0xdc2ff3a253dfc926,
                0x10a75716b3899551,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^2) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xca8d800500fa1bf2,
                0xf0c5d61468b39769,
                0x0e201271ad0d4418,
                0x04290f65bad856e6,
            ])),
            FQ_ZERO,
        ),
        // Fp2::NONRESIDUE^(((q^3) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x365316184e46d97d,
                0x0af7129ed4c96d9f,
                0x659da72fca1009b5,
                0x08116d8983a20d23,
            ])),
            field_new!(Fq, BigInteger256([
                0xb1df4af7c39c1939,
                0x3d9f02878a73bf7f,
                0x9b2220928caf0ae0,
                0x26684515eff054a6,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^4) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x3350c88e13e80b9c,
                0x7dce557cdb5e56b9,
                0x6001b4b8b615564a,
                0x2682e617020217e0,
        ])),
        FQ_ZERO,
        ),
        // Fp2::NONRESIDUE^(((q^5) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x86b76f821b329076,
                0x408bf52b4d19b614,
                0x53dfb9d0d985e92d,
                0x051e20146982d2a7,
            ])),
            field_new!(Fq, BigInteger256([
                0x0fbc9cd47752ebc7,
                0x6d8fffe33415de24,
                0xbef22cf038cf41b9,
                0x15c0edff3c66bf54,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^6) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x68c3488912edefaa,
                0x8d087f6872aabf4f,
                0x51e1a24709081231,
                0x2259d6b14729c0fa,
            ])),
            FQ_ZERO,
        ),
        // Fp2::NONRESIDUE^(((q^7) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x8c84e580a568b440,
                0xcd164d1de0c21302,
                0xa692585790f737d5,
                0x2d7100fdc71265ad,
            ])),
            field_new!(Fq, BigInteger256([
                0x99fdddf38c33cfd5,
                0xc77267ed1213e931,
                0xdc2052142da18f36,
                0x1fbcf75c2da80ad7,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^8) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x71930c11d782e155,
                0xa6bb947cffbe3323,
                0xaa303344d4741444,
                0x2c3b3f0d26594943,
            ])),
            FQ_ZERO,
        ),
        // Fp2::NONRESIDUE^(((q^9) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x05cd75fe8a3623ca,
                0x8c8a57f293a85cee,
                0x52b29e86b7714ea8,
                0x2852e0e95d8f9306,
            ])),
            field_new!(Fq, BigInteger256([
                0x8a41411f14e0e40e,
                0x59e26809ddfe0b0d,
                0x1d2e2523f4d24d7d,
                0x09fc095cf1414b83,
            ])),
        ),
        // Fp2::NONRESIDUE^(((q^10) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0x08cfc388c494f1ab,
                0x19b315148d1373d4,
                0x584e90fdcb6c0213,
                0x09e1685bdf2f8849,
            ])),
            FQ_ZERO,
        ),
        // Fp2::NONRESIDUE^(((q^11) - 1) / 6)
        field_new!(Fq2,
            field_new!(Fq, BigInteger256([
                0xb5691c94bd4a6cd1,
                0x56f575661b581478,
                0x64708be5a7fb6f30,
                0x2b462e5e77aecd82,
            ])),
            field_new!(Fq, BigInteger256([
                0x2c63ef42612a1180,
                0x29f16aae345bec69,
                0xf95e18c648b216a4,
                0x1aa36073a4cae0d4,
            ])),
        ),
    ];
}
