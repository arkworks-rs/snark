use crate::{
    biginteger::BigInteger256 as BigInteger,
    fields::{Fp256, Fp256Parameters, FpParameters},
};

pub type Fr = Fp256<FrParameters>;

pub struct FrParameters;

impl Fp256Parameters for FrParameters {}
impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    // 28948022309329048855892746252171976963322203655955319056773317069363642105857
    const MODULUS: BigInteger = BigInteger([
        0xa14064e200000001,
        0x38aa1276c3f59b9,
        0x0,
        0x4000000000000000,
    ]);

    const R: BigInteger = BigInteger([
        0x1c3ed159fffffffd,
        0xf5601c89bb41f2d3,
        0xffffffffffffffff,
        0x3fffffffffffffff,
    ]);

    const R2: BigInteger = BigInteger([
        0x280c9c4000000010,
        0x91a4409b5400af74,
        0xdd7b28e19094c659,
        0xc8ad9107ccca0e,
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
        0x8388339ffffffed,
        0xbcb60a12f74c5739,
        0xffffffffffffffff,
        0x3fffffffffffffff,
    ]);

    const MODULUS_BITS: u32 = 255;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const TWO_ADICITY: u32 = 33;

    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0xa189b4c6deb5f0b4,
        0x84b1839059d394b6,
        0x62394c58292596e9,
        0x3017cc8a62e742c4,
    ]);

    const REPR_SHAVE_BITS: u32 = 1;

    const INV: u64 = 11619397960441266175;
}
