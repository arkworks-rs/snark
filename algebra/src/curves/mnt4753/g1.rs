use crate::field_new;
use crate::{
    biginteger::BigInteger768,
    curves::{
        models::{ModelParameters, SWModelParameters},
    },
    fields::mnt4753::{Fq, Fr},
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT4G1Parameters;

impl ModelParameters for MNT4G1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for MNT4G1Parameters {
    const COEFF_A: Fq = field_new!(Fq, BigInteger768([
        3553860551672651396,
        2565472393707818253,
        3424927325234966109,
        17487811826058095619,
        15730291918544907998,
        4332070408724822737,
        7212646118208244402,
        12904649141092619460,
        9289117987390442562,
        2254330573517213976,
        3065472942259520298,
        271095073719429,
    ]));

    const COEFF_B: Fq = field_new!(Fq, BigInteger768([
        2672638521926201442,
        17587766986973859626,
        1309143029066506763,
        1756412671449422902,
        5395165286423163724,
        589638022240022974,
        7360845090332416697,
        9829497896347590557,
        9341553552113883496,
        5888515763059971584,
        10173739464651404689,
        456607542322059,
    ]));

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[1];

    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger768([
        0xb99680147fff6f42,
        0x4eb16817b589cea8,
        0xa1ebd2d90c79e179,
        0x0f725caec549c0da,
        0xab0c4ee6d3e6dad4,
        0x9fbca908de0ccb62,
        0x320c3bb713338498,
        0x598b4302d2f00a62,
        0x4074c9cbfd8ca621,
        0x0fa47edb3865e88c,
        0x95455fb31ff9a195,
        0x7b479ec8e242,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);
}

pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger768([
    8680369219962409717,
    12497683146525997170,
    15236963532390397985,
    105054743605190980,
    11580223711797947725,
    5964558218084543687,
    1974179831852844611,
    13386218610606908614,
    9905737029079781539,
    3769381095189112747,
    1226496298859043045,
    409264833279765,
]));

pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger768([
    8458069647833709466,
    16863815841372543189,
    7230518365128572001,
    17250077086581959530,
    15519583030873909149,
    3465247978511199450,
    5738818931561455055,
    12688417287395938373,
    3681991682605141223,
    10698656566578986929,
    10160396483421745615,
    127251255182962,
]));
