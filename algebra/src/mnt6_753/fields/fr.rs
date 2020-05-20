use algebra_core::{
    biginteger::BigInteger768 as BigInteger,
    fields::{FftParameters, Fp768, Fp768Parameters, FpParameters},
};

pub type Fr = Fp768<FrParameters>;

pub struct FrParameters;

impl Fp768Parameters for FrParameters {}
impl FftParameters for FrParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 15;

    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        0x3b079c7556ac378,
        0x2c8c74d04a3f00d4,
        0xd3b001061b90d4cf,
        0x946e77514891b0e6,
        0x79caec8ad6dc9ea1,
        0xbefd780edc81435d,
        0xe093d4dca630b154,
        0x43a0f673199f1c12,
        0x92276c78436253ff,
        0xe249d1cf014fcd24,
        0x96f36471fb7c3ec5,
        0x1080b8906b7c4,
    ]);

    const SMALL_SUBGROUP_BASE: Option<u32> = Some(5);
    const SMALL_SUBGROUP_BASE_ADICITY: Option<u32> = Some(2);
    /// LARGE_SUBGROUP_ROOT_OF_UNITY =
    /// 12249458902762217747626832919710926618510011455364963726393752854649914979954138109976331601455448780251166045203053508523342111624583986869301658366625356826888785691823710598470775453742133593634524619429629803955083254436531
    const LARGE_SUBGROUP_ROOT_OF_UNITY: Option<BigInteger> = Some(BigInteger([
        8926681816978929800,
        10873079436792120119,
        6519893728366769435,
        7899277225737766970,
        8416573500933450083,
        12951641800297678468,
        7093775028595490583,
        14327009285082556021,
        18228411097456927576,
        2823658094446565457,
        1708328092507553067,
        109589007594791,
    ]));
}
impl FpParameters for FrParameters {
    /// MODULUS = 41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888253786114353726529584385201591605722013126468931404347949840543007986327743462853720628051692141265303114721689601
    const MODULUS: BigInteger = BigInteger([
        0x5e9063de245e8001,
        0xe39d54522cdd119f,
        0x638810719ac425f0,
        0x685acce9767254a4,
        0xb80f0da5cb537e38,
        0xb117e776f218059d,
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
        0x98a8ecabd9dc6f42,
        0x91cd31c65a034686,
        0x97c3e4a0cd14572e,
        0x79589819c788b601,
        0xed269c942108976f,
        0x1e0f4d8acf031d68,
        0x320c3bb713338559,
        0x598b4302d2f00a62,
        0x4074c9cbfd8ca621,
        0xfa47edb3865e88c,
        0x95455fb31ff9a195,
        0x7b479ec8e242,
    ]);

    const R2: BigInteger = BigInteger([
        0x84717088cfd190c8,
        0xc7d9ff8e7df03c0a,
        0xa24bea56242b3507,
        0xa896a656a0714c7d,
        0x80a46659ff6f3ddf,
        0x2f47839ef88d7ce8,
        0xa8c86d4604a3b597,
        0xe03c79cac4f7ef07,
        0x2505daf1f4a81245,
        0x8e4605754c381723,
        0xb081f15bcbfdacaf,
        0x2a33e89cb485,
    ]);

    const INV: u64 = 0xf2044cfbe45e7fff;

    const GENERATOR: BigInteger = BigInteger([
        0xa8f627f0e629635e,
        0x202afce346c36872,
        0x85e1ece733493254,
        0x6d76e610664ac389,
        0xdf542f3f04441585,
        0x3aa4885bf6d4dd80,
        0xeb8b63c1c0fffc74,
        0xd2488e985f6cfa4e,
        0xcce1c2a623f7a66a,
        0x2a060f4d5085b19a,
        0xa9111a596408842f,
        0x11ca8d50bf627,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xaf4831ef122f4000,
        0x71ceaa29166e88cf,
        0x31c40838cd6212f8,
        0x342d6674bb392a52,
        0xdc0786d2e5a9bf1c,
        0xd88bf3bb790c02ce,
        0xcce8926cd0ad7bce,
        0x83fedc92f45076c6,
        0xaf5bf47cb64bec39,
        0xdbfccba82dc7d7f6,
        0x88114811777166d6,
        0xe26316c96208,
    ]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T

    /// T = (MODULUS - 1) / 2^S =
    /// 1278640471433073529124274133033466709233725278318907137200424283478556909563327233064541435662546964154604216671394463687571830033251476599169665701965732619291119517454523942352538645255842982596454713491581459512424155325
    const T: BigInteger = BigInteger([
        0x233ebd20c7bc48bd,
        0x4be1c73aa8a459ba,
        0xa948c71020e33588,
        0xfc70d0b599d2ece4,
        0xb3b701e1b4b96a6,
        0xef3b622fceede430,
        0xdb1b33a249b342b5,
        0xb0e60ffb724bd141,
        0x5fdabd6fd1f2d92f,
        0x9b5b6ff32ea0b71f,
        0x882220452045ddc5,
        0x3898c5b25,
    ]);

    /// (T - 1) / 2 =
    /// 639320235716536764562137066516733354616862639159453568600212141739278454781663616532270717831273482077302108335697231843785915016625738299584832850982866309645559758727261971176269322627921491298227356745790729756212077662
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x119f5e9063de245e,
        0x25f0e39d54522cdd,
        0x54a4638810719ac4,
        0x7e38685acce97672,
        0x59db80f0da5cb53,
        0xf79db117e776f218,
        0xed8d99d124d9a15a,
        0xd87307fdb925e8a0,
        0xafed5eb7e8f96c97,
        0xcdadb7f997505b8f,
        0xc41110229022eee2,
        0x1c4c62d92,
    ]);
}
