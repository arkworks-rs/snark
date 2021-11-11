use crate::{
    biginteger::BigInteger384 as BigInteger,
    field_new,
    fields::bn_382::{Fq, Fq2, Fq2Parameters},
    fp2::Fp2Parameters,
    fp6_3over2::{Fp6, Fp6Parameters},
    Field,
};

pub type Fq6 = Fp6<Fq6Parameters>;

#[derive(Clone, Copy)]
pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp2Params = Fq2Parameters;

    // u = sqrt(7)
    // 3 * u has no cube nor square nor sixth root
    /// NONRESIDUE = (2 * U)
    const NONRESIDUE: Fq2 = field_new!(
        Fq2,
        // 0
        field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        // 3
        field_new!(
            Fq,
            BigInteger([
                0xffffffffffffffeb,
                0xffffffe07f020607,
                0x7e04041040c0e07f,
                0x30a99058ab4842e0,
                0xb46b0b54830f49c9,
                0xba0bdc30c3ce150
            ])
        ),
    );

    const FROBENIUS_COEFF_FP6_C1: &'static [Fq2] = &[
        // Fq2(nqr)**(((q^0) - 1) / 3)
        field_new!(
            Fq2,
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
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((q^1) - 1) / 3)
        field_new!(
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
        ),
        // Fq2(nqr)**(((q^2) - 1) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0x43ac10f69cd08675,
                    0xb67658df049b19a2,
                    0xe4f95ea6b5f8ead6,
                    0xd65fd263b6fcff0c,
                    0x6b4f8626c0a86f82,
                    0xb005f791c4b9abd
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((q^3) - 1) / 3)
        field_new!(
            Fq2,
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
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((q^4) - 1) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0xbc53ef09632f7993,
                    0x4989a72cfbc5a71d,
                    0x1bc825e5621f2129,
                    0xcdf1de3d8ddb48ff,
                    0x1f325d26c4458cc2,
                    0x1523ea85ba78a1b6
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((q^5) - 1) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0xbc53ef09632f798c,
                    0x4989a7227b70fe75,
                    0x9b1ed1eacd0a16a9,
                    0xde2a63b0719e09f4,
                    0x600b642eff55005,
                    0x190429c6be8cecd1
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
    ];

    const FROBENIUS_COEFF_FP6_C2: &'static [Fq2] = &[
        // (Fq2(nqr) ** (2/3)) ** (q^0) / (Fq2(nqr) ** (2/3))
        // Fq2(nqr)**(((2q^0) - 2) / 3)
        field_new!(
            Fq2,
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
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((2q^1) - 2) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0x43ac10f69cd08675,
                    0xb67658df049b19a2,
                    0xe4f95ea6b5f8ead6,
                    0xd65fd263b6fcff0c,
                    0x6b4f8626c0a86f82,
                    0xb005f791c4b9abd
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((2q^2) - 2) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0xbc53ef09632f7993,
                    0x4989a72cfbc5a71d,
                    0x1bc825e5621f2129,
                    0xcdf1de3d8ddb48ff,
                    0x1f325d26c4458cc2,
                    0x1523ea85ba78a1b6
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((2q^3) - 2) / 3)
        field_new!(
            Fq2,
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
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((2q^4) - 2) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0x43ac10f69cd08675,
                    0xb67658df049b19a2,
                    0xe4f95ea6b5f8ead6,
                    0xd65fd263b6fcff0c,
                    0x6b4f8626c0a86f82,
                    0xb005f791c4b9abd
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(nqr)**(((2q^5) - 2) / 3)
        field_new!(
            Fq2,
            field_new!(
                Fq,
                BigInteger([
                    0xbc53ef09632f7993,
                    0x4989a72cfbc5a71d,
                    0x1bc825e5621f2129,
                    0xcdf1de3d8ddb48ff,
                    0x1f325d26c4458cc2,
                    0x1523ea85ba78a1b6
                ])
            ),
            field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
    ];

    /// Multiply this element by the quadratic nonresidue 0 + 3 * u.
    fn mul_fp2_by_nonresidue(fe: &Fq2) -> Fq2 {
        // 3 U (c0 + U * c1)
        // == 3*7*c1 + U (3 c0)

        let seven_c1 = Fq2Parameters::mul_fp_by_nonresidue(&fe.c1); // 7*c1
        let mut c0 = seven_c1.clone();
        c0.double_in_place(); // 2*7*c1
        c0 += &seven_c1; // 3*7*c1

        let mut c1 = fe.c0;
        c1.double_in_place();
        c1 += &fe.c0;

        Fq2::new(c0, c1)
    }
}
