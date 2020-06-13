use algebra_core::{
    biginteger::BigInteger256 as BigInteger,
    fields::{FftParameters, Fp256, Fp256Parameters, FpParameters},
};

pub type Fr = Fp256<FrParameters>;

pub struct FrParameters;

impl Fp256Parameters for FrParameters {}
impl FftParameters for FrParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 1;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0xaa9f02ab1d6124de,
        0xb3524a6466112932,
        0x7342261215ac260b,
        0x4d6b87b1da259e2,
    ]);
}
impl FpParameters for FrParameters {
    /// MODULUS = 6554484396890773809930967563523245729705921265872317281365359162392183254199.
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        0xd0970e5ed6f72cb7,
        0xa6682093ccc81082,
        0x6673b0101343b00,
        0xe7db4ea6533afa9,
    ]);

    const MODULUS_BITS: u32 = 252;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 4;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        0x25f80bb3b99607d9,
        0xf315d62f66b6e750,
        0x932514eeeb8814f4,
        0x9a6fc6f479155c6,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        0x67719aa495e57731,
        0x51b0cef09ce3fc26,
        0x69dab7fac026e9a5,
        0x4f6547b8d127688,
    ]);

    const INV: u64 = 0x1ba3a358ef788ef9;

    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        0x720b1b19d49ea8f1,
        0xbf4aa36101f13a58,
        0x5fa8cc968193ccbb,
        0xe70cbdc7dccf3ac,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        7515249040934278747,
        5995434913520945217,
        9454073218019761536,
        522094803716528084,
    ]);

    const T: BigInteger = Self::MODULUS_MINUS_ONE_DIV_TWO;

    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        12980996557321915181,
        2997717456760472608,
        4727036609009880768,
        261047401858264042,
    ]);
}
