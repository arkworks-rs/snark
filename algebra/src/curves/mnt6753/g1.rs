use crate::field_new;
use crate::{
    biginteger::BigInteger768,
    curves::models::{ModelParameters, SWModelParameters},
    fields::mnt6753::{Fq, Fr},
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT6G1Parameters;

impl ModelParameters for MNT6G1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for MNT6G1Parameters {
    // a=11, Montgomery rep.
    const COEFF_A: Fq = field_new!(
        Fq,
        BigInteger768([
            0x4768931cfff9c7d4,
            0xc45e46d6ada96ca0,
            0x479b0bdb0b3c0107,
            0x362a089610f8d41b,
            0xdbafcec2c8a91aaf,
            0x78428b0ff9d96a06,
            0xf2e4472a9080c353,
            0xc9006ed33f0e971c,
            0x0794d9d10bdb7288,
            0x3c1e44cab5419e2c,
            0x49b5fc6c81f4560c,
            0x1c287777c30ba,
        ])
    );
    // b= 1162590899954132115202734022401037471684116770178358464833890823541085\
    // 9267060079819722747939267925389062611062156601938166010098747920378738\
    // 9278326581336254542601154090758161875550558594902533757047280279443155\
    // 01122723426879114
    // Montgom. rep.
    const COEFF_B: Fq = field_new!(
        Fq,
        BigInteger768([
            0x7a85e23c6984298a,
            0xb08f89f10deb6f43,
            0x1ff8d652bcdd2b90,
            0x6fe8b22127f7f097,
            0x57007df447700e3e,
            0x2f8aca277da9258d,
            0x14385d51ca5422fb,
            0x47d8f3de65c79d1d,
            0xfa9ac2fe4bd09711,
            0x9175a8b5ef915920,
            0xf83fa70b67d17c00,
            0x10804126ecf16,
        ])
    );

    const COFACTOR: &'static [u64] = &[1];

    // inverse of cofactor mod group order r
    const COFACTOR_INV: Fr = field_new!(
        Fr,
        BigInteger768([
            0x98A8ECABD9DC6F42,
            0x91CD31C65A034686,
            0x97C3E4A0CD14572E,
            0x79589819C788B601,
            0xED269C942108976F,
            0x1E0F4D8ACF031D68,
            0x320C3BB713338559,
            0x598B4302D2F00A62,
            0x4074C9CBFD8CA621,
            0x0FA47EDB3865E88C,
            0x95455FB31FF9A195,
            0x7B479EC8E242,
        ])
    );

    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);
}

//generator of prime order r
//x =3458420969484235708806261200128850544017070333833944116801482064540723\
// 2681492354777628704146649173606059496596309331847515262279936470308751\
// 6768749271405287219577008822518325905140308790615870178675844188974261\
// 8916006546636728
pub const G1_GENERATOR_X: Fq = field_new!(
    Fq,
    BigInteger768([
        0xe3a856605652f582,
        0xea2ad6adb232d3cc,
        0x006917a62cf94e5d,
        0xb0cf88593f1f8d9c,
        0xdf4294279d098622,
        0xd1805f5f25762cae,
        0x0ce84eed156d448a,
        0x092939a0aaa29f11,
        0x4851f2bd56e6d412,
        0xd6a3f94887cc2c08,
        0xa3870d376b51b4de,
        0x1262a0793b60,
    ])
);

//y=2746050840233196514962660022438213725450297597916837111164092472158912\
// 7725376473514838234361114855175488242007431439074223827742813911899817\
// 9307281122977634480108147641177014035402987649704695003396465633446808\
// 68495474127850569
pub const G1_GENERATOR_Y: Fq = field_new!(
    Fq,
    BigInteger768([
        0xa17be03d3de9993a,
        0xd23d47f834d6e6a7,
        0xc835b816dad2a400,
        0xb067d33661cbda12,
        0x34917ee69c71eaa3,
        0x69dcbdab27c304e6,
        0xeea1a2a6d6c76015,
        0x5e60253078c4f3e3,
        0x1eee46f45880e189,
        0xd8de606656eb5e1c,
        0xbf48f43a878dac3a,
        0x37d7e759d51c,
    ])
);
