use crate::curves::mnt4753::MNT4_753Parameters;
use crate::curves::models::mnt4::MNT4Parameters;
use crate::field_new;
use crate::{
    biginteger::BigInteger768,
    curves::models::{ModelParameters, SWModelParameters},
    fields::mnt4753::{Fq, Fq2, Fr},
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT4G2Parameters;

impl ModelParameters for MNT4G2Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

/// MUL_BY_A_C0 = NONRESIDUE * COEFF_A
pub const MUL_BY_A_C0: Fq = field_new!(
    Fq,
    BigInteger768([
        0xeb354e6121cdccad,
        0x9589bfe5ea49ae4f,
        0xb12cc53998b3d124,
        0x7883d83c06c22baa,
        0xd828782cb96edc7,
        0x35e68bd867a8d558,
        0xe0860ea489bec5bd,
        0xe034be400ffa8f19,
        0xf4d51fe5c821f43d,
        0x8ee1bf11396bcec0,
        0xb819c73cb726c963,
        0x23dae1639e4b,
    ])
);

/// MUL_BY_A_C1 = NONRESIDUE * COEFF_A
pub const MUL_BY_A_C1: Fq = field_new!(
    Fq,
    BigInteger768([
        0xeb354e6121cdccad,
        0x9589bfe5ea49ae4f,
        0xb12cc53998b3d124,
        0x7883d83c06c22baa,
        0xd828782cb96edc7,
        0x35e68bd867a8d558,
        0xe0860ea489bec5bd,
        0xe034be400ffa8f19,
        0xf4d51fe5c821f43d,
        0x8ee1bf11396bcec0,
        0xb819c73cb726c963,
        0x23dae1639e4b,
    ])
);

impl SWModelParameters for MNT4G2Parameters {
    // quadratic twist E' of the G1-curve E: y^2= x^3 + a + b
    // E': y^2 = x^3 + a*alpha*x + b*alpha*X.
    // over F2 = Fq[X]/(X^2-alpha),
    const COEFF_A: Fq2 = MNT4_753Parameters::TWIST_COEFF_A;
    const COEFF_B: Fq2 = field_new!(
        Fq2,
        field_new!(Fq, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
        field_new!(
            Fq,
            BigInteger768([
                0xd1f842ef859c74ef,
                0x9d45480c3873434a,
                0xa5566d8d8d841941,
                0xc0f99a3682ad8bae,
                0xe4b39f099a706e70,
                0xce59a66ebad048e2,
                0x93fe1794e855b79e,
                0x957322b9044da5e8,
                0x836b3c49c9f33d5d,
                0x3ea13c16b209ced3,
                0x79f8ca52b73621ea,
                0x1a2270165e15a,
            ])
        ),
    );

    /// cofactor of G2, native integer representation
    /// COFACTOR = 41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888049094905534395567574915333486969589229856772141392370549616644545554517640527237829320384324374366385444967219201
    const COFACTOR: &'static [u64] = &[
        0xe41950da08bd0001,
        0x789a0f8d4a18e8ee,
        0xf04c9f26f687f44a,
        0x16d5a05cb84b6ea3,
        0x313250b76d85d63a,
        0xafc372c51bd661a0,
        0x99d124d9a15af79d,
        0x7fdb925e8a0ed8d,
        0x5eb7e8f96c97d873,
        0xb7f997505b8fafed,
        0x10229022eee2cdad,
        0x1c4c62d92c411,
    ];
    /// cofactor^{-1} mod r, for projection of curve onto G2,
    /// MG representation
    const COFACTOR_INV: Fr = field_new!(
        Fr,
        BigInteger768([
            0x1a14ef94372dbc2a,
            0x6e01a14d0f55ad00,
            0x5955ab3920afde4d,
            0xe7982fd78cbf4332,
            0xecbf393ce1701610,
            0xd111cd07a49d61b4,
            0xe58145271adb10a9,
            0x2e22af0c3ca18713,
            0x35d277c2206aed22,
            0xfb6c4c412f6bacd0,
            0x68e1c109cfc51649,
            0x4747058a5c42,
        ])
    );

    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(elt: &Fq2) -> Fq2 {
        field_new!(Fq2, MUL_BY_A_C0 * &elt.c0, MUL_BY_A_C1 * &elt.c1,)
    }
}

const G2_GENERATOR_X: Fq2 = field_new!(Fq2, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
const G2_GENERATOR_Y: Fq2 = field_new!(Fq2, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

// generator of prime order r
// x = (c0+ c1*X)
// c0 =2948396511084314467570336474470883652464396010553860807886250839750244\
// 7349913068434941060515343254862580437318493682762113105361632548148204\
// 8060521140087313727573896453838919822112450139651752134560664525878695\
// 19098351487925167
// c1=1970601131963017239107607962479975394815850677122214748623799532192544\
// 3331396169656568431378974558350664383559981183980668976846806019030432\
// 3891691379539889908020005810789940082839677683482759739215981662748576\
// 31001635633631000
pub const G2_GENERATOR_X_C0: Fq = field_new!(
    Fq,
    BigInteger768([
        0x64cd9e87e5f14e2d,
        0x6a123355ee938785,
        0xdb197417e1887231,
        0xf5199b0d7e333053,
        0x397da434e85b78a7,
        0xd1117417b290f004,
        0x3f8ccbdf316d6964,
        0x1ea26a53c24e4162,
        0x4fa40c8be29a9276,
        0x3c355554caad2580,
        0x5b05c21a27b7acc7,
        0x13635a0d01b78,
    ])
);

pub const G2_GENERATOR_X_C1: Fq = field_new!(
    Fq,
    BigInteger768([
        0xb56fc312dc34c98b,
        0x2b8029c87df25f6,
        0x9217aa6ceb0cf808,
        0xe67355775c8eb87e,
        0x90eb471ebb74c1b1,
        0x6bebff63e88338c2,
        0xde8489295782a103,
        0xf1e11281b99054d1,
        0x71e05664c68aa32,
        0x6ec60806cb661af7,
        0x31facd7fa4614cca,
        0x16aea1ba33dc0,
    ])
);

// y = (c0+ c1*Y)
// c0 =3994015267076051965394032031482732794199314140370833866692520428208447\
// 7074754642625849927569427860786384998614863651207257467076192649385174\
// 1080858031687438034917805685033693170931911017795340353772663001850993\
// 18717465441820654,
// c1=1760863742496439573704129137375665713960730644019373180410245701172669\
// 0702169238966996114255971643893157857311132388792357391583164125870757\
// 5410090350414694633665287985939528847459876974030564887446038294374489\
// 27398468360797245
pub const G2_GENERATOR_Y_C0: Fq = field_new!(
    Fq,
    BigInteger768([
        0x438d727bf4c5002e,
        0x220b34b2b9daee2f,
        0xa567a1375a9e2a27,
        0x36739870b33ba70f,
        0xb058c55679b63f3e,
        0xdb048df87997b3b7,
        0xf64a68ade535340f,
        0xe526639d49ef3eff,
        0xd52be2d6e4bee8fd,
        0x8e46b4ca897b87bb,
        0x4f6af38904883c28,
        0x1202d1e47ccef,
    ])
);

pub const G2_GENERATOR_Y_C1: Fq = field_new!(
    Fq,
    BigInteger768([
        0xf06781d5bec3ed74,
        0xa41dbc2a99750c11,
        0x6a393d84e066ddfc,
        0xbbf8387b3a74937a,
        0xecb6da0ba28e9879,
        0x380c74d14f4e2d84,
        0x5c089d226f9c345d,
        0x92e8c3c2f8040454,
        0xdebf6ce50f1d3555,
        0xe659a934e501a154,
        0xa35de638cd06f1c4,
        0x11e7b5c581f,
    ])
);
