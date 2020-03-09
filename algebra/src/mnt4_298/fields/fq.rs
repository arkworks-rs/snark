use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    fields::{Fp320, Fp320Parameters, FpParameters},
};

pub type Fq = Fp320<FqParameters>;

pub struct FqParameters;

impl Fp320Parameters for FqParameters {}
impl FpParameters for FqParameters {
    type BigInt = BigInteger;

    // MODULUS = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
    const MODULUS: BigInteger = BigInteger([
        0xc90cd65a71660001,
        0x41a9e35e51200e12,
        0xcaeec9635d1330ea,
        0xa266249da7b0548e,
        0x3bcf7bcd473,
    ]);

    const MODULUS_BITS: u32 = 298;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 22;

    const R: BigInteger = BigInteger([
        0x18c31a7b5863845c,
        0xe9de7a15e3b68df5,
        0xc5df858728faab40,
        0x29184098647b5197,
        0x1c1223d33c3,
    ]);

    const R2: BigInteger = BigInteger([
        0x65acec5613d220,
        0xa266a1adbf2bc893,
        0x66bd7673318850e1,
        0x1f32e014ad38d47b,
        0x224f0918a34,
    ]);

    const INV: u64 = 0xb071a1b67165ffff;

    const GENERATOR: BigInteger = BigInteger([
        0x259ae5b7c4d1ca15,
        0xbc20e3dfe73f0ac3,
        0x97505c422d1f08e7,
        0x49d149cf165e1b2c,
        0x3a87fe6a0cc,
    ]);

    const TWO_ADICITY: u32 = 17;

    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0x884ce85c8d89f2b9,
        0x8366528dcef9a167,
        0x8a465859c7d431ff,
        0xce4c49d76adbcbc5,
        0x39fc98494e1,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x64866b2d38b30000,
        0x20d4f1af28900709,
        0x657764b1ae899875,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    // T = (MODULUS - 1) / 2^S =
    // 3630998887399759870554727551674258816109656366292531779446068791017229177993437198515
    const T: BigInteger = BigInteger([
        0x70964866b2d38b3,
        0x987520d4f1af2890,
        0x2a47657764b1ae89,
        0x6a39d133124ed3d8,
        0x1de7bde,
    ]);

    // (T - 1) / 2 =
    // 1815499443699879935277363775837129408054828183146265889723034395508614588996718599257
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x384b24335969c59,
        0xcc3a906a78d79448,
        0x1523b2bbb258d744,
        0x351ce899892769ec,
        0xef3def,
    ]);
}
