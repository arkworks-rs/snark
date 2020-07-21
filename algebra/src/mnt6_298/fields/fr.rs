use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    fields::{FftParameters, Fp320, Fp320Parameters, FpParameters},
};

pub type Fr = Fp320<FrParameters>;

pub struct FrParameters;

impl Fp320Parameters for FrParameters {}
impl FftParameters for FrParameters {
    type BigInt = BigInteger;

    const TWO_ADICITY: u32 = 17;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        9821480371597472441u64,
        9468346035609379175u64,
        9963748368231707135u64,
        14865337659602750405u64,
        3984815592673u64,
    ]);

    const SMALL_SUBGROUP_BASE: Option<u32> = Some(7);
    const SMALL_SUBGROUP_BASE_ADICITY: Option<u32> = Some(2);

    /// LARGE_SUBGROUP_ROOT_OF_UNITY = x * g
    /// where x = (n - 1) / 2^17 / 7^2
    /// and represent this value in the Montgomery residue form.
    /// I.e., write
    /// 381811485921190977554243339163030148371175054922689353173385941180422489253833691237722982
    /// * R
    /// = 260534023778902228073198316993669317435810479439368306496187170459125001342456918103569322
    const LARGE_SUBGROUP_ROOT_OF_UNITY: Option<BigInteger> = Some(BigInteger([
        7711798843682337706u64,
        16456007754393011187u64,
        7470854640069402569u64,
        10767969225751706229u64,
        2250015743691u64,
    ]));
}
impl FpParameters for FrParameters {
    /// MODULUS = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
        14487189785281953793u64,
        4731562877756902930u64,
        14622846468719063274u64,
        11702080941310629006u64,
        4110145082483u64,
    ]);

    const MODULUS_BITS: u32 = 298;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 22;

    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        1784298994435064924u64,
        16852041090100268533u64,
        14258261760832875328u64,
        2961187778261111191u64,
        1929014752195u64,
    ]);

    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        28619103704175136u64,
        11702218449377544339u64,
        7403203599591297249u64,
        2248105543421449339u64,
        2357678148148u64,
    ]);

    const INV: u64 = 12714121028002250751u64;

    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        2709730703260633621u64,
        13556085429182073539u64,
        10903316137158576359u64,
        5319113788683590444u64,
        4022235209932u64,
    ]);

    #[rustfmt::skip]
    const T: BigInteger = BigInteger([
        0x70964866b2d38b3,
        0x987520d4f1af2890,
        0x2a47657764b1ae89,
        0x6a39d133124ed3d8,
        0x1de7bde,
    ]);

    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x384b24335969c59,
        0xcc3a906a78d79448,
        0x1523b2bbb258d744,
        0x351ce899892769ec,
        0xef3def,
    ]);

    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x64866b2d38b30000,
        0x20d4f1af28900709,
        0x657764b1ae899875,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);
}
