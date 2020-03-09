use crate::mnt4_753::{self, Fq, Fq2, Fr, FQ_ZERO, FR_ONE, G1_COEFF_A_NON_RESIDUE};
use algebra_core::{
    biginteger::BigInteger768,
    curves::{
        mnt4,
        mnt4::MNT4Parameters,
        models::{ModelParameters, SWModelParameters},
    },
    field_new,
};

pub type G2Affine = mnt4::G2Affine<mnt4_753::Parameters>;
pub type G2Projective = mnt4::G2Projective<mnt4_753::Parameters>;
pub type G2Prepared = mnt4::G2Prepared<mnt4_753::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

/// MUL_BY_A_C0 = NONRESIDUE * COEFF_A
#[rustfmt::skip]
pub const MUL_BY_A_C0: Fq = G1_COEFF_A_NON_RESIDUE;

/// MUL_BY_A_C1 = NONRESIDUE * COEFF_A
#[rustfmt::skip]
pub const MUL_BY_A_C1: Fq = G1_COEFF_A_NON_RESIDUE;

impl SWModelParameters for Parameters {
    const COEFF_A: Fq2 = mnt4_753::Parameters::TWIST_COEFF_A;
    // B coefficient of MNT4-753 G2 =
    // ```
    // mnt4753_twist_coeff_b = mnt4753_Fq2(mnt4753_Fq::zero(),
    //                                  mnt4753_G1::coeff_b * mnt4753_Fq2::non_residue);
    // non_residue = mnt4753_Fq2::non_residue = mnt4753_Fq("13");
    //  = (ZERO, G1_B_COEFF * NON_RESIDUE);
    //  =
    //  (0, 39196523001581428369576759982967177918859161321667605855515469914917622337081756705006832951954384669101573360625169461998308377011601613979275218690841934572954991361632773738259652003389826903175898479855893660378722437317212)
    // ```
    #[rustfmt::skip]
    const COEFF_B: Fq2 = field_new!(Fq2,
        FQ_ZERO,
        field_new!(Fq, BigInteger768([
            15129916544657421551,
            11332543254671606602,
            11913830318987286849,
            13905314883394440110,
            16479690325073358448,
            14869098639251228898,
            10663986895980443550,
            10768989312009479656,
            9469728929095040349,
            4512954369775881939,
            8788997129423430122,
            459763387588954,
        ])),
    );

    /// COFACTOR =
    /// 1
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[1];

    /// COFACTOR^(-1) mod r =
    /// 1
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = FR_ONE;

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(elt: &Fq2) -> Fq2 {
        field_new!(Fq2, MUL_BY_A_C0 * &elt.c0, MUL_BY_A_C1 * &elt.c1,)
    }
}

const G2_GENERATOR_X: Fq2 = field_new!(Fq2, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
const G2_GENERATOR_Y: Fq2 = field_new!(Fq2, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

// Generator of G2
// These are two Fq elements each because X and Y (and Z) are elements of Fq^2
// X = 22367666623321080720060256844679369841450849258634485122226826668687008928557241162389052587294939105987791589807198701072089850184203060629036090027206884547397819080026926412256978135536735656049173059573120822105654153939204,
// 19674349354065582663569886390557105215375764356464013910804136534831880915742161945711267871023918136941472003751075703860943205026648847064247080124670799190998395234694182621794580160576822167228187443851233972049521455293042,
// Y = 6945425020677398967988875731588951175743495235863391886533295045397037605326535330657361771765903175481062759367498970743022872494546449436815843306838794729313050998681159000579427733029709987073254733976366326071957733646574,
// 17406100775489352738678485154027036191618283163679980195193677896785273172506466216232026037788788436442188057889820014276378772936042638717710384987239430912364681046070625200474931975266875995282055499803236813013874788622488,
#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger768([
    6027998872234822269,
    18316510890646896372,
    12312422472129155053,
    5734876120361342991,
    15269557702005117850,
    2881476951079071365,
    18286169900244335081,
    9121647201144159462,
    10143201970379789827,
    8274325828999741070,
    16324620334829045351,
    206265970291999
]));

#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger768([
    8103915349469236029,
    3264578743345366807,
    11420369645465736410,
    1082524563485784168,
    9521451516267963565,
    947729167637853217,
    980516398344887244,
    9338197375178553870,
    12279316873997193841,
    7682952744956285803,
    5880726343431033740,
    393937351760210
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger768([
    9473954153981524211,
    17210941273038840246,
    7908789794149338247,
    14886051892012552213,
    15661744656776151922,
    12925231519316918877,
    8556358774098608262,
    12018580250090935208,
    12085312901555144792,
    16429889710969777184,
    9875783130742037645,
    390595340749132
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger768([
    1961701679311783916,
    14701334664186945036,
    4897259467765012901,
    9133165644380234616,
    5806216571361183036,
    15854960112739074469,
    3305583651306700066,
    3531430907605494499,
    7284796601341566706,
    3173472022185515256,
    1610147509686374719,
    280032154091959
]));
