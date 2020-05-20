use algebra_core::{
    biginteger::BigInteger768 as BigInteger,
    fields::{FftParameters, Fp768, Fp768Parameters, FpParameters},
};

pub type Fq = Fp768<FqParameters>;

pub struct FqParameters;

impl Fp768Parameters for FqParameters {}
impl FftParameters for FqParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 30;

    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0x307f66b297671883,
        0xd72a7f2b1e645f4e,
        0x67079daa9a902283,
        0xf33f7620a86c668b,
        0x8878570d66464c12,
        0xa557af5b524f522b,
        0x5fafa3f6ef19319d,
        0x1eb9e04110a65629,
        0x3f96feb3c639a0b0,
        0x4d4fe37df3ffd732,
        0xadc831bd55bcf3e9,
        0x1b9f32a8bd6ab,
    ]);
}
impl FpParameters for FqParameters {
    /// MODULUS = 41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160001
    const MODULUS: BigInteger = BigInteger([
        0xd90776e240000001,
        0x4ea099170fa13a4f,
        0xd6c381bc3f005797,
        0xb9dff97634993aa4,
        0x3eebca9429212636,
        0xb26c5c28c859a99b,
        0x99d124d9a15af79d,
        0x7fdb925e8a0ed8d,
        0x5eb7e8f96c97d873,
        0xb7f997505b8fafed,
        0x10229022eee2cdad,
        0x1c4c62d92c411,
    ]);

    const MODULUS_BITS: u32 = 753;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 15;

    const R: BigInteger = BigInteger([
        0xb99680147fff6f42,
        0x4eb16817b589cea8,
        0xa1ebd2d90c79e179,
        0xf725caec549c0da,
        0xab0c4ee6d3e6dad4,
        0x9fbca908de0ccb62,
        0x320c3bb713338498,
        0x598b4302d2f00a62,
        0x4074c9cbfd8ca621,
        0xfa47edb3865e88c,
        0x95455fb31ff9a195,
        0x7b479ec8e242,
    ]);

    const R2: BigInteger = BigInteger([
        0x3f9c69c7b7f4c8d1,
        0x70a50fa9ee48d127,
        0xcdbe6702009569cb,
        0x6bd8c6c6c49edc38,
        0x7955876cc35ee94e,
        0xc7285529be54a3f4,
        0xded52121ecec77cf,
        0x99be80f2ee12ee8e,
        0xc8a0ff01493bdcef,
        0xacc27988f3d9a316,
        0xd9e817a8fb44b3c9,
        0x5b58037e0e4,
    ]);

    const INV: u64 = 0xc90776e23fffffff;

    const GENERATOR: BigInteger = BigInteger([
        0xeee0a5d37ff6635e,
        0xff458536cfa1cff4,
        0x659af978d8169ab0,
        0x1f1841c24780e3f1,
        0x602213036dcfef3a,
        0xd1d5c8f39d72db20,
        0xeb8b63c1c0ffefab,
        0xd2488e985f6cfa4e,
        0xcce1c2a623f7a66a,
        0x2a060f4d5085b19a,
        0xa9111a596408842f,
        0x11ca8d50bf627,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xec83bb7120000000,
        0xa7504c8b87d09d27,
        0x6b61c0de1f802bcb,
        0x5ceffcbb1a4c9d52,
        0x9f75e54a1490931b,
        0xd9362e14642cd4cd,
        0xcce8926cd0ad7bce,
        0x83fedc92f45076c6,
        0xaf5bf47cb64bec39,
        0xdbfccba82dc7d7f6,
        0x88114811777166d6,
        0xe26316c96208,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    /// T = (MODULUS - 1) / 2^S =
    /// 39021010480745652133919498688765463538626870065884617224134041854204007249857398469987226430131438115069708760723898631821547688442835449306011425196003537779414482717728302293895201885929702287178426719326440397855625
    const T: BigInteger = BigInteger([
        0x3e84e93f641ddb89,
        0xfc015e5d3a82645c,
        0xd264ea935b0e06f0,
        0xa48498dae77fe5d8,
        0x2166a66cfbaf2a50,
        0x856bde76c9b170a3,
        0xa283b63667449366,
        0xb25f61cc1ff6e497,
        0x6e3ebfb57adfa3e5,
        0xbb8b36b6dfe65d41,
        0xb64b1044408a408b,
        0x71318,
    ]);

    /// (T - 1) / 2 =
    /// 19510505240372826066959749344382731769313435032942308612067020927102003624928699234993613215065719057534854380361949315910773844221417724653005712598001768889707241358864151146947600942964851143589213359663220198927812
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x1f42749fb20eedc4,
        0x7e00af2e9d41322e,
        0x69327549ad870378,
        0x52424c6d73bff2ec,
        0x90b353367dd79528,
        0x42b5ef3b64d8b851,
        0xd141db1b33a249b3,
        0xd92fb0e60ffb724b,
        0xb71f5fdabd6fd1f2,
        0xddc59b5b6ff32ea0,
        0x5b25882220452045,
        0x3898c,
    ]);
}
