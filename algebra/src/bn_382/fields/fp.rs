use algebra_core::{
    biginteger::BigInteger384 as BigInteger,
    fields::{FftParameters, Fp384, Fp384Parameters},
};

pub type Fp = Fp384<FpParameters>;

pub struct FpParameters;

impl Fp384Parameters for FpParameters {}
impl FftParameters for FpParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 67;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0xdb510d8c5d0d218f,
        0x447119a2f8d5e310,
        0x1373332ba33d5a84,
        0xb830356347b45dbb,
        0x851efb96cb691ec1,
        0x141037c57e9d0173,
    ]);

}

impl algebra_core::fields::FpParameters for FpParameters {
    const MODULUS: BigInteger = BigInteger([
        0x1,
        0x1800c1818,
        0x2012246d22424120,
        0xb48a3614289b0901,
        0x71503c69b09dbf88,
        0x2404893fdad8878e,
    ]);

    const R: BigInteger = BigInteger([
        0xfffffffffffffff9,
        0xfffffff57fab5757,
        0x1f8101041030381f,
        0x10388572e3c2c0f8,
        0xe6ce591c2bafc343,
        0x3e03f4104144b1a,
    ]);

    const R2: BigInteger = BigInteger([
        0xaa7b14a53b610887,
        0xb22034140d119ca9,
        0xe10d2796937ba75,
        0xe52454bf8b810402,
        0x1b4eec3d89fc0fd3,
        0xbc857aea27171f7,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x0,
        0xc0060c0c,
        0x9009123691212090,
        0x5a451b0a144d8480,
        0x38a81e34d84edfc4,
        0x1202449fed6c43c7,
    ]);

    const T: BigInteger = BigInteger([
        0x30018303,
        0x2402448da4484824,
        0x169146c285136120,
        0xce2a078d3613b7f1,
        0x4809127fb5b10f1,
        0x0,
    ]);

    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x1800c181,
        0x12012246d2242412,
        0x8b48a3614289b090,
        0xe71503c69b09dbf8,
        0x2404893fdad8878,
        0x0,
    ]);

    // GENERATOR = 7
    const GENERATOR: BigInteger = BigInteger([
        0xffffffffffffffcf,
        0xffffffb67daf6367,
        0xdc87071c715188df,
        0x718ba6243a5346c8,
        0x4fa46fc531ce56d5,
        0x1b21bac71c8e0dbc,
    ]);

    const MODULUS_BITS: u32 = 382;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    // Check this
    const REPR_SHAVE_BITS: u32 = 2;

    const INV: u64 = 18446744073709551615;
}
