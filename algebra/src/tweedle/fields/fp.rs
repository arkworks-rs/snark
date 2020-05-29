use algebra_core::{
    biginteger::BigInteger256 as BigInteger,
    fields::{FftParameters, Fp256, Fp256Parameters},
};

pub type Fp = Fp256<FpParameters>;

pub struct FpParameters;

impl Fp256Parameters for FpParameters {}
impl FftParameters for FpParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 33;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0xa189b4c6deb5f0b4,
        0x84b1839059d394b6,
        0x62394c58292596e9,
        0x3017cc8a62e742c4,
    ]);

}

impl algebra_core::fields::FpParameters for FpParameters {
    // 28948022309329048855892746252171976963322203655955319056773317069363642105857
    const MODULUS: BigInteger = BigInteger([
        11619397960441266177,
        255193519591741881,
        0,
        4611686018427387904,
    ]);

    const R: BigInteger = BigInteger([
        2035294266095304701,
        17681163514934325971,
        18446744073709551615,
        4611686018427387903,
    ]);

    const R2: BigInteger = BigInteger([
        2885853259929485328,
        10494584067553537908,
        15959394653775906393,
        56485833754855950,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xd0a0327100000000,
        0x1c55093b61facdc,
        0x0,
        0x2000000000000000,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T
    const T: BigInteger = BigInteger([0xb61facdcd0a03271, 0x1c55093, 0x0, 0x20000000]);

    const T_MINUS_ONE_DIV_TWO: BigInteger =
        BigInteger([0xdb0fd66e68501938, 0xe2a849, 0x0, 0x10000000]);

    // GENERATOR = 5
    const GENERATOR: BigInteger = BigInteger([
        592367636130562029,
        13598067201466455865,
        18446744073709551615,
        4611686018427387903,
    ]);

    const MODULUS_BITS: u32 = 255;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 1;

    const INV: u64 = 11619397960441266175;
}
