use algebra_core::{
    biginteger::BigInteger768 as BigInteger,
    field_new,
    fields::{FftParameters, Fp768, Fp768Parameters, FpParameters},
};

pub type Fq = Fp768<FqParameters>;

pub struct FqParameters;

pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));

impl Fp768Parameters for FqParameters {}
impl FftParameters for FqParameters {
    type BigInt = BigInteger;

    // The internal representation of this type is six 64-bit unsigned
    // integers in little-endian order. Values are always in
    // Montgomery form; i.e., Scalar(a) = aR mod p, with R=2^768.

    // (MODULUS - 1) % 2^TWO_ADICITY == 0
    const TWO_ADICITY: u32 = 1;

    // least_quadratic_nonresidue(MODULUS) in Sage.
    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        17481284903592032950u64,
        10104133845767975835u64,
        8607375506753517913u64,
        13706168424391191299u64,
        9580010308493592354u64,
        14241333420363995524u64,
        6665632285037357566u64,
        5559902898979457045u64,
        15504799981718861253u64,
        8332096944629367896u64,
        18005297320867222879u64,
        58811391084848524u64,
    ]);
}
impl FpParameters for FqParameters {
    /// MODULUS = 6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068299
    #[rustfmt::skip]
    const MODULUS: BigInteger = BigInteger([
       0xf49d00000000008b,
       0xe6913e6870000082,
       0x160cf8aeeaf0a437,
       0x98a116c25667a8f8,
       0x71dcd3dc73ebff2e,
       0x8689c8ed12f9fd90,
       0x03cebaff25b42304,
       0x707ba638e584e919,
       0x528275ef8087be41,
       0xb926186a81d14688,
       0xd187c94004faff3e,
       0x122e824fb83ce0a
    ]);

    const MODULUS_BITS: u32 = 761;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    // gap to 64-bit machine word
    const REPR_SHAVE_BITS: u32 = 7;

    // 2^768 % MODULUS
    #[rustfmt::skip]
    const R: BigInteger = BigInteger([
        144959613005956565u64,
        6509995272855063783u64,
        11428286765660613342u64,
        15738672438262922740u64,
        17071399330169272331u64,
        13899911246788437003u64,
        12055474021000362245u64,
        2545351818702954755u64,
        8887388221587179644u64,
        5009280847225881135u64,
        15539704305423854047u64,
        23071597697427581u64,
    ]);

    // R^2
    #[rustfmt::skip]
    const R2: BigInteger = BigInteger([
        14305184132582319705u64,
        8868935336694416555u64,
        9196887162930508889u64,
        15486798265448570248u64,
        5402985275949444416u64,
        10893197322525159598u64,
        3204916688966998390u64,
        12417238192559061753u64,
        12426306557607898622u64,
        1305582522441154384u64,
        10311846026977660324u64,
        48736111365249031u64,
    ]);

    // (-1/MODULUS) % 2^64
    const INV: u64 = 744663313386281181u64;

    /// GENERATOR = 2
    // primitive_root(MODULUS)
    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([
        289919226011913130u64,
        13019990545710127566u64,
        4409829457611675068u64,
        13030600802816293865u64,
        15696054586628993047u64,
        9353078419867322391u64,
        5664203968291172875u64,
        5090703637405909511u64,
        17774776443174359288u64,
        10018561694451762270u64,
        12632664537138156478u64,
        46143195394855163u64,
    ]);

    // (MODULUS - 1) / 2
    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x7a4e800000000045,
        0xf3489f3438000041,
        0x0b067c577578521b,
        0x4c508b612b33d47c,
        0x38ee69ee39f5ff97,
        0x4344e476897cfec8,
        0x81e75d7f92da1182,
        0xb83dd31c72c2748c,
        0x29413af7c043df20,
        0x5c930c3540e8a344,
        0x68c3e4a0027d7f9f,
        0x9174127dc1e705,
    ]);

    // T =
    // 3445725192157866269698394841137828771239834456268075054756895080104811711121745868043841591644705843820432283876893306725580879560277123879674755849562650799475802549689254425186271815711798397975949850214984556421382456559534149
    // (MODULUS - 1) / 2 ^ TWO_ADICITY
    #[rustfmt::skip]
    const T: BigInteger = BigInteger([
        0x7a4e800000000045,
        0xf3489f3438000041,
        0x0b067c577578521b,
        0x4c508b612b33d47c,
        0x38ee69ee39f5ff97,
        0x4344e476897cfec8,
        0x81e75d7f92da1182,
        0xb83dd31c72c2748c,
        0x29413af7c043df20,
        0x5c930c3540e8a344,
        0x68c3e4a0027d7f9f,
        0x9174127dc1e705,
    ]);

    // (T - 1)/2 =
    // 1722862596078933134849197420568914385619917228134037527378447540052405855560872934021920795822352921910216141938446653362790439780138561939837377924781325399737901274844627212593135907855899198987974925107492278210691228279767074
    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xbd27400000000022,
        0xf9a44f9a1c000020,
        0x05833e2bbabc290d,
        0xa62845b09599ea3e,
        0x1c7734f71cfaffcb,
        0x21a2723b44be7f64,
        0x40f3aebfc96d08c1,
        0x5c1ee98e39613a46,
        0x14a09d7be021ef90,
        0xae49861aa07451a2,
        0xb461f250013ebfcf,
        0x48ba093ee0f382,
    ]);
}
