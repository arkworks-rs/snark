use crate::{
    biginteger::BigInteger768 as BigInteger,
    bw6_761::{Fq, Fq3, Fq3Parameters, FQ_ONE, FQ_ZERO},
    field_new,
    fields::fp6_2over3::{Fp6, Fp6Parameters},
};

pub type Fq6 = Fp6<Fq6Parameters>;

pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp3Params = Fq3Parameters;

    /// NONRESIDUE = (0, 1, 0)
    #[rustfmt::skip]
    const NONRESIDUE: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP6_C1: &'static [Fq] = &[
        field_new!(Fq, BigInteger([
            0x0202ffffffff85d5,
            0x5a5826358fff8ce7,
            0x9e996e43827faade,
            0xda6aff320ee47df4,
            0xece9cb3e1d94b80b,
            0xc0e667a25248240b,
            0xa74da5bfdcad3905,
            0x2352e7fe462f2103,
            0x7b56588008b1c87c,
            0x45848a63e711022f,
            0xd7a81ebb9f65a9df,
            0x51f77ef127e87d,
        ])),
        field_new!(Fq, BigInteger([
            0x8cfcb51bd8404a93,
            0x495e69d68495a383,
            0xd23cbc9234705263,
            0x8d2b4c2b5fcf4f52,
            0x6a798a5d20c612ce,
            0x3e825d90eb6c2443,
            0x772b249f2c9525fe,
            0x521b2ed366e4b9bb,
            0x84abb49bd7c4471d,
            0x907062359c0f17e3,
            0x3385e55030cc6f12,
            0x3f11a3a41a2606,
        ])),
        field_new!(Fq, BigInteger([
            0x7f96b51bd840c549,
            0xd59782096496171f,
            0x49b046fd9ce14bbc,
            0x4b6163bba7527a56,
            0xef6c92fb771d59f1,
            0x0425bedbac1dfdc7,
            0xd3ac39de759c0ffd,
            0x9f43ed0e063a81d0,
            0x5bd7d20b4f9a3ce2,
            0x0411f03c36cf5c3c,
            0x2d658fd49661c472,
            0x1100249ae760b93,
        ])),
        field_new!(Fq, BigInteger([
            0xf29a000000007ab6,
            0x8c391832e000739b,
            0x77738a6b6870f959,
            0xbe36179047832b03,
            0x84f3089e56574722,
            0xc5a3614ac0b1d984,
            0x5c81153f4906e9fe,
            0x4d28be3a9f55c815,
            0xd72c1d6f77d5f5c5,
            0x73a18e069ac04458,
            0xf9dfaa846595555f,
            0xd0f0a60a5be58c,
        ])),
        field_new!(Fq, BigInteger([
            0x67a04ae427bfb5f8,
            0x9d32d491eb6a5cff,
            0x43d03c1cb68051d4,
            0x0b75ca96f69859a5,
            0x0763497f5325ec60,
            0x48076b5c278dd94d,
            0x8ca3965ff91efd06,
            0x1e6077657ea02f5d,
            0xcdd6c153a8c37724,
            0x28b5b634e5c22ea4,
            0x9e01e3efd42e902c,
            0xe3d6815769a804,
        ])),
        field_new!(Fq, BigInteger([
            0x75064ae427bf3b42,
            0x10f9bc5f0b69e963,
            0xcc5cb1b14e0f587b,
            0x4d3fb306af152ea1,
            0x827040e0fccea53d,
            0x82640a1166dbffc8,
            0x30228120b0181307,
            0xd137b92adf4a6748,
            0xf6aaa3e430ed815e,
            0xb514282e4b01ea4b,
            0xa422396b6e993acc,
            0x12e5db4d0dc277,
        ])),
    ];
}
