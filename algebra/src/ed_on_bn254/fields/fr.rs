use algebra_core::{
    biginteger::BigInteger256 as BigInteger,
    fields::{FftParameters, Fp256, Fp256Parameters, FpParameters},
};

pub type Fr = Fp256<FrParameters>;

pub struct FrParameters;

impl Fp256Parameters for FrParameters {}
impl FftParameters for FrParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 4;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0x1721ada8d4d27255,
        0xcda0f5264e0e35bb,
        0x961a936922086fe6,
        0x1ab00857387dd52,
    ]);
}
impl FpParameters for FrParameters {
    /// MODULUS = 2736030358979909402780800718157159386076813972158567259200215660948447373041
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        0x677297dc392126f1,
        0xab3eedb83920ee0a,
        0x370a08b6d0302b0b,
        0x60c89ce5c263405,
    ]);

    const MODULUS_BITS: u32 = 251;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 5;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        0x073315dea08f9c76,
        0xe7acffc6a098f24b,
        0xf85a9201d818f015,
        0x1f16424e1bb7724,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        0x35e44abee7ecb21e,
        0x74646cacf5f84ec4,
        0xe472df203faa158f,
        0x445b524f1ba50a8,
    ]);

    const INV: u64 = 0x532ce5aebc48f5ef;

    #[rustfmt::skip]
    /// GENERATOR = 31
    const GENERATOR: BigInteger = BigInteger([
        0x3c284f376f3993d1,
        0x08bc9d93705cf8b8,
        0x239d5fcbd9538f3e,
        0x5ca4836185b994b,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x33b94bee1c909378,
        0xd59f76dc1c907705,
        0x9b85045b68181585,
        0x30644e72e131a02,
    ]);

    const T: BigInteger = BigInteger([
        0xa677297dc392126f,
        0xbab3eedb83920ee0,
        0x5370a08b6d0302b0,
        0x60c89ce5c26340,
    ]);

    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x533b94bee1c90937,
        0x5d59f76dc1c90770,
        0x29b85045b6818158,
        0x30644e72e131a0,
    ]);
}
