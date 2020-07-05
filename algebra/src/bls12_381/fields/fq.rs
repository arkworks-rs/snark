use algebra_core::{
    biginteger::BigInteger384 as BigInteger,
    field_new,
    fields::{FftParameters, Fp384, Fp384Parameters, FpParameters},
};

pub type Fq = Fp384<FqParameters>;

pub struct FqParameters;

impl Fp384Parameters for FqParameters {}
impl FftParameters for FqParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 1;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0x43f5fffffffcaaae,
        0x32b7fff2ed47fffd,
        0x7e83a49a2e99d69,
        0xeca8f3318332bb7a,
        0xef148d1ea0f4c069,
        0x40ab3263eff0206,
    ]);
}
impl FpParameters for FqParameters {
    /// MODULUS = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ]);

    const MODULUS_BITS: u32 = 381;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 3;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        0xf4df1f341c341746,
        0xa76e6a609d104f1,
        0x8de5476c4c95b6d5,
        0x67eb88a9939d83c0,
        0x9a793e85b519952d,
        0x11988fe592cae3aa,
    ]);

    const INV: u64 = 0x89f3fffcfffcfffd;

    // GENERATOR = 2
    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        0x321300000006554f,
        0xb93c0018d6c40005,
        0x57605e0db0ddbb51,
        0x8b256521ed1f9bcb,
        0x6cf28d7901622c03,
        0x11ebab9dbb81e28c,
    ]);

    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xdcff7fffffffd555,
        0xf55ffff58a9ffff,
        0xb39869507b587b12,
        0xb23ba5c279c2895f,
        0x258dd3db21a5d66b,
        0xd0088f51cbff34d,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    #[rustfmt::skip]
    const T: BigInteger = BigInteger([
        0xdcff7fffffffd555,
        0xf55ffff58a9ffff,
        0xb39869507b587b12,
        0xb23ba5c279c2895f,
        0x258dd3db21a5d66b,
        0xd0088f51cbff34d,
    ]);

    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xee7fbfffffffeaaa,
        0x7aaffffac54ffff,
        0xd9cc34a83dac3d89,
        0xd91dd2e13ce144af,
        0x92c6e9ed90d2eb35,
        0x680447a8e5ff9a6,
    ]);
}

pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0]));
