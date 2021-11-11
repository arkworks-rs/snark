use crate::{
    biginteger::BigInteger384 as BigInteger,
    field_new,
    fields::bn_382::fq::{Fq, FQ_ONE, FQ_ZERO},
    fp2::{Fp2, Fp2Parameters},
    Field,
};

pub type Fq2 = Fp2<Fq2Parameters>;

pub struct Fq2Parameters;

impl Fp2Parameters for Fq2Parameters {
    type Fp = Fq;

    /// NONRESIDUE = 7
    const NONRESIDUE: Fq = field_new!(
        Fq,
        BigInteger([
            0xffffffffffffffcf,
            0xffffffb67daf6367,
            0x7b5eb425ec6cb67f,
            0x718ba6243a5346b6,
            0x4fa46fc531ce56d5,
            0x1b21bac71c8e0dbc
        ])
    );

    // U = sqrt(7)
    /// QUADRATIC_NONRESIDUE = (0 + 2 * U)
    const QUADRATIC_NONRESIDUE: (Fq, Fq) = (
        // 0
        field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        // 2
        field_new!(
            Fq,
            BigInteger([
                0xfffffffffffffff2,
                0xffffffeaff56aeaf,
                0xfead580ad5d5eaff,
                0x20710ae5c78581ea,
                0xcd9cb238575f8686,
                0x7c07e8208289635
            ])
        ),
    );

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        field_new!(
            Fq,
            BigInteger([
                0xfffffffffffffff9,
                0xfffffff57fab5757,
                0x7f56ac056aeaf57f,
                0x10388572e3c2c0f5,
                0xe6ce591c2bafc343,
                0x3e03f4104144b1a
            ])
        ),
        // Fq(7)**(((q^1) - 1) / 2)
        field_new!(
            Fq,
            BigInteger([
                0x8,
                0xc0060c0c0,
                0xc1848c18180c00,
                0xa451b0a144d8480c,
                0x8a81e34d84edfc45,
                0x202449fed6c43c73
            ])
        ),
    ];

    #[inline(always)]
    fn mul_fp_by_nonresidue(fe: &Self::Fp) -> Self::Fp {
        // times 7
        let mut result = fe.clone();
        result.double_in_place(); // 2x
        result += fe; // 3x
        result.double_in_place(); // 6x
        result += fe; // 7x
        result
    }
}

pub const FQ2_ZERO: Fq2 = field_new!(Fq2, FQ_ZERO, FQ_ZERO);
pub const FQ2_ONE: Fq2 = field_new!(Fq2, FQ_ONE, FQ_ZERO);

#[cfg(test)]
mod test {
    #![allow(unused_imports)]
    use super::*;
    use crate::{Field, UniformRand};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_fq_mul_nonresidue() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let seven: u32 = 7;
        let non_residue = Fq::from(seven);

        for _ in 0..1000 {
            let mut a = Fq::rand(&mut rng);
            let mut b = a;
            a = Fq2Parameters::mul_fp_by_nonresidue(&a);
            b *= &non_residue;

            assert_eq!(a, b);
        }
    }
}
