use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    fields::{FftParameters, Fp320, Fp320Parameters, FpParameters},
};

pub type Fq = Fp320<FqParameters>;

pub struct FqParameters;

impl Fp320Parameters for FqParameters {}
impl FftParameters for FqParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 34;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0x818b361df1af7be4,
        0x2ae2750d46a53957,
        0x5784a8fe792c5f8a,
        0xf9bd39c0cdcf1bb6,
        0x6a24a0f8a8,
    ]);
}
impl FpParameters for FqParameters {
    /// MODULUS = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        0xbb4334a400000001,
        0xfb494c07925d6ad3,
        0xcaeec9635cf44194,
        0xa266249da7b0548e,
        0x3bcf7bcd473,
    ]);

    const MODULUS_BITS: u32 = 298;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 22;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        0xc3177aefffbb845c,
        0x9b80c702f9961788,
        0xc5df8dcdac70a85a,
        0x29184098647b5197,
        0x1c1223d33c3,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        0x465a743c68e0596b,
        0x34f9102adb68371,
        0x4bbd6dcf1e3a8386,
        0x2ff00dced8e4b6d,
        0x149bb44a342,
    ]);

    const INV: u64 = 0xbb4334a3ffffffff;

    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        0xb1ddfacffd532b94,
        0x25e295ff76674008,
        0x8f00647b48958d36,
        0x1159f37d4e0fddb2,
        0x2977770b3d1,
    ]);

    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xdda19a5200000000,
        0x7da4a603c92eb569,
        0x657764b1ae7a20ca,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    #[rustfmt::skip]
    const T: BigInteger = BigInteger([
        0xe4975ab4eed0cd29,
        0xd73d10653ed25301,
        0x69ec1523b2bbb258,
        0x3def351ce8998927,
        0xef,
    ]);

    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xf24bad5a77686694,
        0x6b9e88329f692980,
        0xb4f60a91d95dd92c,
        0x9ef79a8e744cc493,
        0x77,
    ]);
}
