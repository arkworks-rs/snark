use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    fields::{Fp320, Fp320Parameters, FpParameters},
};

pub type Fr = Fp320<FrParameters>;

pub struct FrParameters;

impl Fp320Parameters for FrParameters {}
impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    // MODULUS = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
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

    const R: BigInteger = BigInteger([
        0xc3177aefffbb845c,
        0x9b80c702f9961788,
        0xc5df8dcdac70a85a,
        0x29184098647b5197,
        0x1c1223d33c3,
    ]);

    const R2: BigInteger = BigInteger([
        0x465a743c68e0596b,
        0x34f9102adb68371,
        0x4bbd6dcf1e3a8386,
        0x2ff00dced8e4b6d,
        0x149bb44a342,
    ]);

    const INV: u64 = 0xbb4334a3ffffffff;

    const GENERATOR: BigInteger = BigInteger([
        0xd5b8b973fb73ca15,
        0x748c22fd9269a44a,
        0x9750e8f0e8cd62f1,
        0x49d149cf165e1b2c,
        0x3a87fe6a0cc,
    ]);

    const TWO_ADICITY: u32 = 34;

    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0xdb5ef7190ef806f6,
        0x6c35bdadf87862d8,
        0x10ac552d29ade8e0,
        0x356f85ccd04b9b0f,
        0x2f5d38f9dae,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xdda19a5200000000,
        0x7da4a603c92eb569,
        0x657764b1ae7a20ca,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    // T = (MODULUS - 1) / 2^S =
    // 27702323054502562488973446286577291993024111641153199339359284829066871159442729
    const T: BigInteger = BigInteger([
        0xe4975ab4eed0cd29,
        0xd73d10653ed25301,
        0x69ec1523b2bbb258,
        0x3def351ce8998927,
        0xef,
    ]);

    // (T - 1) / 2 =
    // 13851161527251281244486723143288645996512055820576599669679642414533435579721364
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xf24bad5a77686694,
        0x6b9e88329f692980,
        0xb4f60a91d95dd92c,
        0x9ef79a8e744cc493,
        0x77,
    ]);
}
