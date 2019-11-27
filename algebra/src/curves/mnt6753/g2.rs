use crate::field_new;
use crate::{
    biginteger::BigInteger768,
    curves::{
        models::{ModelParameters, SWModelParameters, mnt6::MNT6Parameters},
    },
    fields::mnt6753::{Fq, Fq3, Fr},
};
use crate::curves::mnt6753::MNT6_753Parameters;

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT6b7G2Parameters;

impl ModelParameters for MNT6b7G2Parameters {
    type BaseField = Fq3;
    type ScalarField = Fr;
}

/// MUL_BY_A_C0 = NONRESIDUE * COEFF_A
// mnt6753_twist_mul_by_a_c0  1ac145a9a6f568978374c40a3a9bd658d0b923634ec9ea036423d447f7762531b87989f5736b66da39e53e7fab86d32a06017e8d6ee3afb57f89503f8a764110e9fd4acb4d2bbb006710f05909f6c5dc71052d9fb63c49733ad687fbb9612
pub const MUL_BY_A_C0: Fq = field_new!(Fq, BigInteger768([
    0x9733ad687fbb9612,
    0x5dc71052d9fb63c4,
    0xb006710f05909f6c,
    0x110e9fd4acb4d2bb,
    0xfb57f89503f8a764,
    0x32a06017e8d6ee3a,
    0x6da39e53e7fab86d,
    0x531b87989f5736b6,
    0xa036423d447f7762,
    0x658d0b923634ec9e,
    0x8978374c40a3a9bd,
    0x1ac145a9a6f56,
]));

/// MUL_BY_A_C1 = NONRESIDUE * COEFF_A
// mnt6753_twist_mul_by_a_c1  1ac145a9a6f568978374c40a3a9bd658d0b923634ec9ea036423d447f7762531b87989f5736b66da39e53e7fab86d32a06017e8d6ee3afb57f89503f8a764110e9fd4acb4d2bbb006710f05909f6c5dc71052d9fb63c49733ad687fbb9612
pub const MUL_BY_A_C1: Fq = field_new!(Fq, BigInteger768([
    0x9733ad687fbb9612,
    0x5dc71052d9fb63c4,
    0xb006710f05909f6c,
    0x110e9fd4acb4d2bb,
    0xfb57f89503f8a764,
    0x32a06017e8d6ee3a,
    0x6da39e53e7fab86d,
    0x531b87989f5736b6,
    0xa036423d447f7762,
    0x658d0b923634ec9e,
    0x8978374c40a3a9bd,
    0x1ac145a9a6f56,
]));

/// MUL_BY_A_C2 = COEFF_A
//pub const MUL_BY_A_C2: Fq = MNT6b7G1Parameters::COEFF_A;
// mnt6753_twist_mul_by_a_c2  1c287777c30ba49b5fc6c81f4560c3c1e44cab5419e2c0794d9d10bdb7288c9006ed33f0e971cf2e4472a9080c35378428b0ff9d96a06dbafcec2c8a91aaf362a089610f8d41b479b0bdb0b3c0107c45e46d6ada96ca04768931cfff9c7d4
pub const MUL_BY_A_C2: Fq = field_new!(Fq, BigInteger768([
    0x4768931cfff9c7d4,
    0xc45e46d6ada96ca0,
    0x479b0bdb0b3c0107,
    0x362a089610f8d41b,
    0xdbafcec2c8a91aaf,
    0x78428b0ff9d96a06,
    0xf2e4472a9080c353,
    0xc9006ed33f0e971c,
    0x794d9d10bdb7288,
    0x3c1e44cab5419e2c,
    0x49b5fc6c81f4560c,
    0x1c287777c30ba,
]));

impl SWModelParameters for MNT6b7G2Parameters {
    const COEFF_A: Fq3 = MNT6_753Parameters::TWIST_COEFF_A;
    const COEFF_B: Fq3 = field_new!(Fq3, 
// G2b0  bb87b9524d9649ecccabdcaf81f3f034b3ef25e0b4da8c56eb14b6676b0be66023aaebce2eee43852569e97bb318dc6c86bdb42aa36f437eaa061b09b76673bfcea77c0ff6a1571c2b24a37ed1abbe6756d13b566a072d93ef4b08adc8e8
        field_new!(Fq, BigInteger768([
            0x2d93ef4b08adc8e8,
            0xbe6756d13b566a07,
            0x571c2b24a37ed1ab,
            0x73bfcea77c0ff6a1,
            0x437eaa061b09b766,
            0xdc6c86bdb42aa36f,
            0x43852569e97bb318,
            0xe66023aaebce2eee,
            0x8c56eb14b6676b0b,
            0xf034b3ef25e0b4da,
            0x49ecccabdcaf81f3,
            0xbb87b9524d96,
        ])),
        field_new!(Fq, BigInteger768([0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger768([0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0])),
    );

    /// 1755483545388786116744270475466687259186947712032004459714210070280389500116987496124098574823389466285978151140155508638765729019174599527183600372094760023144398285325863550664578643924584541949466179502227232245309952839189635010671372908411609248348904807785904229403747495114436660255866932060472369629692502198423138429922875792635236729929780298333055698257230963645509826963717287902205842627121011526048163097042046361575549171961352924692480000
    const COFACTOR: &'static [u64] = &[
        0xf791c4a6c0000000,
        0x6f2920fb3d823ec,
        0x1d491f05951364e8,
        0x14d431154f3deeb0,
        0xb22ff5f7a2d737ff,
        0x7c9a2c218777f2a9,
        0xbfee11b09da07297,
        0x69d6d25c051da042,
        0xce73086230450ba5,
        0x20263932a197b6eb,
        0x24193d3622676e8d,
        0xfdf90b0a130158dc,
        0xf06212969fa553ca,
        0xaeab6a88d8766394,
        0xb367861ecb8a3fa,
        0xb59bf579e2771844,
        0xe655d99477c1a210,
        0x4d6ae2386a6cc1eb,
        0xe4add8dcc2b66488,
        0xa80fe42e26865a9f,
        0x51e513f8f3102a9c,
        0x2f48ca819e99ddc4,
        0xcb8078713bbf018b,
        0x320cc6a58,
    ];

    /// COFACTOR^(-1) mod r =
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger768([
        0x82b3d2d2e7467187,
        0x30153d6ae73f0d65,
        0x6d3f936a7e253292,
        0xf30f5f7211e53371,
        0x14c41dc7045f7bbf,
        0xeb6421cf7999b27a,
        0xd6168bfe7bc42088,
        0x801a97b2485a7407,
        0xa2d037ce6f8ac1ba,
        0xe44d9966a694b385,
        0x5fe7a77550c0744,
        0x561283152c84,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(elt: &Fq3) -> Fq3 {
        field_new!(Fq3, 
            MUL_BY_A_C0 * &elt.c1,
            MUL_BY_A_C1 * &elt.c2,
            MUL_BY_A_C2 * &elt.c0,
        )
    }
}

const G2_GENERATOR_X: Fq3 = field_new!(Fq3, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1, G2_GENERATOR_X_C2);
const G2_GENERATOR_Y: Fq3 = field_new!(Fq3, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1, G2_GENERATOR_Y_C2);

// G2x0  177aea09a1f71839241c93cbf7d91d24255ae7c04a7a8a5fcdf7a6d8b7968bd9bd46c566d4cbe5873d597e56f0487a4f32415b8a2189baa1fd9a4bdef1668493c3c21ffed0366de083979f2c359b227dfb31106118698b14220caee543428
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger768([
    0xb14220caee543428,
    0x27dfb31106118698,
    0xde083979f2c359b2,
    0x493c3c21ffed0366,
    0xaa1fd9a4bdef1668,
    0xa4f32415b8a2189b,
    0x5873d597e56f0487,
    0xbd9bd46c566d4cbe,
    0xa5fcdf7a6d8b7968,
    0xd24255ae7c04a7a8,
    0x839241c93cbf7d91,
    0x177aea09a1f71,
]));

// G2x1  12d3bb942a644ff8c958e1472e4d43d785c23326c726e33024814010fd4250c74832a3a2b1a756321b4d7f65121332c3805c4b887fa30e8838252779fbebcc8de70c9e75e12f3f463969df64afb05d59ca6e9dffed6c8916a4eb2ff7ac11a
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger768([
    0x916a4eb2ff7ac11a,
    0xd59ca6e9dffed6c8,
    0xf463969df64afb05,
    0xc8de70c9e75e12f3,
    0xe8838252779fbebc,
    0x2c3805c4b887fa30,
    0x6321b4d7f6512133,
    0xc74832a3a2b1a75,
    0x33024814010fd425,
    0x3d785c23326c726e,
    0xff8c958e1472e4d4,
    0x12d3bb942a644,
]));

// G2x2  1a644a7fb0fe164905525c33a28f7828f5add25a4697c2b44333f83f000148d19280ce3a41005c744c72e17ea89f76607f68a31749df74112d188ed29cc72faf7bf5b149435ca785a588edec7b65c955e55ac1effae7b473e9cc2ebb3c488
pub const G2_GENERATOR_X_C2: Fq = field_new!(Fq, BigInteger768([
    0x473e9cc2ebb3c488,
    0x955e55ac1effae7b,
    0x785a588edec7b65c,
    0xfaf7bf5b149435ca,
    0x4112d188ed29cc72,
    0x6607f68a31749df7,
    0xc744c72e17ea89f7,
    0x8d19280ce3a41005,
    0x2b44333f83f00014,
    0x828f5add25a4697c,
    0x64905525c33a28f7,
    0x1a644a7fb0fe1,
]));

// G2y0  134292ce7239d83fb773b722c87416a9c98ccd49629dd02320552a153d9fb62b00daf959607d9c215ea37ee10430eb759e35cd6a47d7a9197d518a1a168223d2ec3bfaa8d2c0f41884f71348f2e2c8142929aa5af8f4909dd7d7268df8be4
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger768([
    0x9dd7d7268df8be4,
    0x8142929aa5af8f49,
    0x41884f71348f2e2c,
    0x3d2ec3bfaa8d2c0f,
    0x9197d518a1a16822,
    0xb759e35cd6a47d7a,
    0xc215ea37ee10430e,
    0x62b00daf959607d9,
    0x2320552a153d9fb,
    0x6a9c98ccd49629dd,
    0x83fb773b722c8741,
    0x134292ce7239d,
]));

// G2y1  b989fe6211f492d87301ea3e508704e530f5145092c2aaa316ace9544f20421c4650aaeb64296ccda659fa4011f7b02e427951c10205d25626fd52245537ce60d691c186e88ddede563eb3616cf254485e1f10822cb3971f000db290cad4
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger768([
    0x971f000db290cad4,
    0x54485e1f10822cb3,
    0xdede563eb3616cf2,
    0xce60d691c186e88d,
    0xd25626fd52245537,
    0xb02e427951c10205,
    0x6ccda659fa4011f7,
    0x421c4650aaeb6429,
    0xaaa316ace9544f20,
    0x4e530f5145092c2,
    0x92d87301ea3e5087,
    0xb989fe6211f4,
]));

// G2y2  1048083573eb8e6a2dc9d5072c959ecb892f7c87f3424c04e5ed641bbf239a066daaf13ec37f09a979eb00ec43d13f754350edd164274ad77155a7a8a919a70166e62d3f62ee72da29e127bdbbfae24d41dcb9b8a7699c4ffabe9bc7f7fb
pub const G2_GENERATOR_Y_C2: Fq = field_new!(Fq, BigInteger768([
    0x9c4ffabe9bc7f7fb,
    0xe24d41dcb9b8a769,
    0x72da29e127bdbbfa,
    0xa70166e62d3f62ee,
    0x4ad77155a7a8a919,
    0x3f754350edd16427,
    0x9a979eb00ec43d1,
    0x9a066daaf13ec37f,
    0x4c04e5ed641bbf23,
    0x9ecb892f7c87f342,
    0x8e6a2dc9d5072c95,
    0x1048083573eb,
]));