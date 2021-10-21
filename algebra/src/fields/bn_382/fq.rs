use crate::{
    biginteger::BigInteger384 as BigInteger,
    field_new,
    fields::{Fp384, Fp384Parameters, FpParameters},
};

pub type Fq = Fp384<FqParameters>;

pub struct FqParameters;

// const U : [u64; a] = [0, 1073873924];
// const SIX_U_PLUS_2_NAF : [i8;98] = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
// 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, -1, 0, 1];
// const CUBIC_NONRESIDUE_TO_Q_MINUS_1_OVER_2 : Fq2 = Fq2(
// field_new!(Fq, BigInteger([ 0x16b744a7d72fb912, 0x8db76da14b98776d,
// 0xd7d0fda03758326c, 0x9a05f3af0ce04699, 0x1c8a66ecb161efb2,
// 0x13a9f1d5f1261bfe ])), Fq::zero()
// );
impl Fp384Parameters for FqParameters {}
impl FpParameters for FqParameters {
    type BigInt = BigInteger;

    // MODULUS = 5543634365110765627805495722742127385843376434033820803592568747918351978899288491582778380528407187068941959692289
    const MODULUS: BigInteger = BigInteger([
        0x1,
        0x1800c1818,
        0x8018309183030180,
        0xb48a3614289b0901,
        0x71503c69b09dbf88,
        0x2404893fdad8878e,
    ]);

    const MODULUS_BITS: u32 = 382;

    // Check this
    const REPR_SHAVE_BITS: u32 = 2;

    const R: BigInteger = BigInteger([
        0xfffffffffffffff9,
        0xfffffff57fab5757,
        0x7f56ac056aeaf57f,
        0x10388572e3c2c0f5,
        0xe6ce591c2bafc343,
        0x3e03f4104144b1a,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    const R2: BigInteger = BigInteger([
        0xc79c121e98884701,
        0xfd75271b6a2e235d,
        0x1530439e68fe657,
        0xf6b7a72ebfbdbfb,
        0x50c6c2ce8f44951b,
        0x17fe189b54066561,
    ]);

    const INV: u64 = 18446744073709551615;

    // GENERATOR = 14
    const GENERATOR: BigInteger = BigInteger([
        0xffffffffffffff9d,
        0xffffff6b7b52aeb7,
        0x76a537ba55d66b7f,
        0x2e8d16344c0b846b,
        0x2df8a320b2feee22,
        0x123eec4e5e4393ea,
    ]);

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const TWO_ADICITY: u32 = 67;

    const ROOT_OF_UNITY: Self::BigInt = BigInteger([
        0xe38be9090411d7d0,
        0x579d9745d8f8468b,
        0x4a5514233c9850c5,
        0xa7c5be912557804a,
        0xc69e67da380310d4,
        0x136e8eef9cf4445b,
    ]);

    const T: BigInteger = BigInteger([
        0x30018303,
        0x3003061230606030,
        0x169146c285136120,
        0xce2a078d3613b7f1,
        0x4809127fb5b10f1,
        0x0,
    ]);
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x1800c181,
        0x1801830918303018,
        0x8b48a3614289b090,
        0xe71503c69b09dbf8,
        0x2404893fdad8878,
        0x0,
    ]);
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x0,
        0xc0060c0c,
        0xc00c1848c18180c0,
        0x5a451b0a144d8480,
        0x38a81e34d84edfc4,
        0x1202449fed6c43c7,
    ]);
}

pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0]));
