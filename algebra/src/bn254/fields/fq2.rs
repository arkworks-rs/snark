use super::*;
use algebra_core::{biginteger::BigInteger256 as BigInteger, field_new, fields::*};

pub type Fq2 = Fp2<Fq2Parameters>;

pub struct Fq2Parameters;

impl Fp2Parameters for Fq2Parameters {
    type Fp = Fq;

    /// NONRESIDUE = -1
    #[rustfmt::skip]
    const NONRESIDUE: Fq = field_new!(Fq, BigInteger([
        0x68c3488912edefaa,
        0x8d087f6872aabf4f,
        0x51e1a24709081231,
        0x2259d6b14729c0fa,
    ]));

    /// QUADRATIC_NONRESIDUE = U+2
    #[rustfmt::skip]
    const QUADRATIC_NONRESIDUE: (Fq, Fq) = (
        field_new!(Fq, BigInteger([
            12014063508332092218u64,
            1509222997478479483u64,
            14762033076929465432u64,
            2023505479389396574u64,
        ])),
        field_new!(Fq, BigInteger([
            202099033278250856u64,
            8885205928937022213u64,
            5545221690922665192u64,
            39800542322357402u64,
        ])),
    );

    /// Coefficients for the Frobenius automorphism.
    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // NONRESIDUE**(((q^0) - 1) / 2)
        field_new!(Fq, BigInteger([
            0xd35d438dc58f0d9d,
            0x0a78eb28f5c70b3d,
            0x666ea36f7879462c,
            0xe0a77c19a07df2f,
        ])),
        // NONRESIDUE**(((q^1) - 1) / 2)
        field_new!(Fq, BigInteger([
            0x68c3488912edefaa,
            0x8d087f6872aabf4f,
            0x51e1a24709081231,
            0x2259d6b14729c0fa,
        ])),
    ];

    #[inline(always)]
    fn mul_fp_by_nonresidue(fe: &Self::Fp) -> Self::Fp {
        -(*fe)
    }
}

pub const FQ2_ZERO: Fq2 = field_new!(Fq2, FQ_ZERO, FQ_ZERO);
pub const FQ2_ONE: Fq2 = field_new!(Fq2, FQ_ONE, FQ_ZERO);
