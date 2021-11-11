use crate::{
    biginteger::BigInteger256 as BigInteger,
    field_new,
    fields::{Fp256, Fp256Parameters, FpParameters},
};

pub struct FqParameters;

pub type Fq = Fp256<FqParameters>;

impl Fp256Parameters for FqParameters {}
impl FpParameters for FqParameters {
    type BigInt = BigInteger;

    // 28948022309329048855892746252171976963322203655954433126947083963168578338817
    const MODULUS: BigInteger = BigInteger([
        0x842cafd400000001,
        0x38aa127696286c9,
        0x0,
        0x4000000000000000,
    ]);

    const R: BigInteger = BigInteger([
        0x7379f083fffffffd,
        0xf5601c89c3d86ba3,
        0xffffffffffffffff,
        0x3fffffffffffffff,
    ]);

    const R2: BigInteger = BigInteger([
        0x8595fa8000000010,
        0x7e16a565c6895230,
        0xf4c0e6fcb03aa0a2,
        0xc8ad9106886013,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xc21657ea00000000,
        0x1c55093b4b14364,
        0x0,
        0x2000000000000000,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    const T: BigInteger = BigInteger([0xda58a1b2610b2bf5, 0xe2a849, 0x0, 0x10000000]);

    const T_MINUS_ONE_DIV_TWO: BigInteger =
        BigInteger([0xed2c50d9308595fa, 0x715424, 0x0, 0x8000000]);

    // GENERATOR = 5
    const GENERATOR: BigInteger = BigInteger([
        0x30aef343ffffffed,
        0xbcb60a132dafff0b,
        0xffffffffffffffff,
        0x3fffffffffffffff,
    ]);

    const MODULUS_BITS: u32 = 255;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const TWO_ADICITY: u32 = 34;

    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0xb5373d390b45cde8,
        0x1437982d8ca321e6,
        0x5f6fd892c6494e7e,
        0x19f5297fb35e2ae1,
    ]);

    // Check this
    const REPR_SHAVE_BITS: u32 = 1;

    const INV: u64 = 9524180637049683967;
}

pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0]));
