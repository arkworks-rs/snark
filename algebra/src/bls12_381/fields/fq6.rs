use crate::bls12_381::*;
use algebra_core::{biginteger::BigInteger384, field_new, fields::*};

pub type Fq6 = Fp6<Fq6Parameters>;

#[derive(Clone, Copy)]
pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp2Params = Fq2Parameters;

    /// NONRESIDUE = (U + 1)
    #[rustfmt::skip]
    const NONRESIDUE: Fq2 = field_new!(Fq2,
        field_new!(Fq, BigInteger384([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493,
        ])),
        field_new!(Fq, BigInteger384([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493,
        ])),
    );

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C1: &'static [Fq2] = &[
        // Fq2(u + 1)**(((q^0) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x760900000002fffd,
                0xebf4000bc40c0002,
                0x5f48985753c758ba,
                0x77ce585370525745,
                0x5c071a97a256ec6d,
                0x15f65ec3fa80e493,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((q^1) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
            field_new!(Fq, BigInteger384([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ])),
        ),
        // Fq2(u + 1)**(((q^2) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x30f1361b798a64e8,
                0xf3b8ddab7ece5a2a,
                0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b,
                0x3636b76660701c6e,
                0x51ba4ab241b6160,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((q^3) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
            field_new!(Fq, BigInteger384([
                0x760900000002fffd,
                0xebf4000bc40c0002,
                0x5f48985753c758ba,
                0x77ce585370525745,
                0x5c071a97a256ec6d,
                0x15f65ec3fa80e493,
            ])),
        ),
        // Fq2(u + 1)**(((q^4) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((q^5) - 1) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
            field_new!(Fq, BigInteger384([
                0x30f1361b798a64e8,
                0xf3b8ddab7ece5a2a,
                0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b,
                0x3636b76660701c6e,
                0x51ba4ab241b6160,
            ])),
        ),
    ];

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C2: &'static [Fq2] = &[
        // Fq2(u + 1)**(((2q^0) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x760900000002fffd,
                0xebf4000bc40c0002,
                0x5f48985753c758ba,
                0x77ce585370525745,
                0x5c071a97a256ec6d,
                0x15f65ec3fa80e493,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((2q^1) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
                0x14e4f04fe2db9068,
                0x14e56d3f1564853a,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((2q^2) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((2q^3) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x43f5fffffffcaaae,
                0x32b7fff2ed47fffd,
                0x7e83a49a2e99d69,
                0xeca8f3318332bb7a,
                0xef148d1ea0f4c069,
                0x40ab3263eff0206,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((2q^4) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0x30f1361b798a64e8,
                0xf3b8ddab7ece5a2a,
                0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b,
                0x3636b76660701c6e,
                0x51ba4ab241b6160,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
        // Fq2(u + 1)**(((2q^5) - 2) / 3)
        field_new!(Fq2,
            field_new!(Fq, BigInteger384([
                0xecfb361b798dba3a,
                0xc100ddb891865a2c,
                0xec08ff1232bda8e,
                0xd5c13cc6f1ca4721,
                0x47222a47bf7b5c04,
                0x110f184e51c5f59,
            ])),
            field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        ),
    ];

    /// Multiply this element by the quadratic nonresidue 1 + u.
    /// Make this generic.
    fn mul_fp2_by_nonresidue(fe: &Fq2) -> Fq2 {
        let mut copy = *fe;
        let t0 = copy.c0;
        copy.c0 -= &fe.c1;
        copy.c1 += &t0;
        copy
    }
}
