use crate::{
    biginteger::BigInteger256 as BigInteger,
    fields::{Fp256, Fp256Parameters, FpParameters},
};

pub type Fr = Fp256<FrParameters>;

pub struct FrParameters;

impl Fp256Parameters for FrParameters {}

impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    // MODULUS = 8444461749428370424248824938781546531375899335154063827935233455917409239041
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        725501752471715841u64,
        6461107452199829505u64,
        6968279316240510977u64,
        1345280370688173398u64,
    ]);

    const MODULUS_BITS: u32 = 253;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 3;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        9015221291577245683u64,
        8239323489949974514u64,
        1646089257421115374u64,
        958099254763297437u64,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        2726216793283724667u64,
        14712177743343147295u64,
        12091039717619697043u64,
        81024008013859129u64,
    ]);

    const INV: u64 = 725501752471715839u64;

    // GENERATOR = 11
    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        1855201571499933546u64,
        8511318076631809892u64,
        6222514765367795509u64,
        1122129207579058019u64,
    ]);

    const TWO_ADICITY: u32 = 47;

    #[rustfmt::skip]
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0x3c3d3ca739381fb2,
        0x9a14cda3ec99772b,
        0xd7aacc7c59724826,
        0xd1ba211c5cc349c,
    ]);

    /// (r - 1)/2 =
    /// 4222230874714185212124412469390773265687949667577031913967616727958704619520
    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x8508c00000000000,
        0xacd53b7f68000000,
        0x305a268f2e1bd800,
        0x955b2af4d1652ab,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where r - 1 = 2^s * t

    /// t = (r - 1) / 2^s =
    /// 60001509534603559531609739528203892656505753216962260608619555
    #[rustfmt::skip]
    const T: BigInteger = BigInteger([
        0xedfda00000021423,
        0x9a3cb86f6002b354,
        0xcabd34594aacc168,
        0x2556,
    ]);

    /// (t - 1) / 2 =
    /// 30000754767301779765804869764101946328252876608481130304309777
    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x76fed00000010a11,
        0x4d1e5c37b00159aa,
        0x655e9a2ca55660b4,
        0x12ab,
    ]);
}
