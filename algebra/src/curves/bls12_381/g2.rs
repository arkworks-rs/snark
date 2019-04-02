use crate::{
    biginteger::{BigInteger256, BigInteger384},
    curves::{
        bls12::{G2Affine as Bls12G2Affine, G2Prepared, G2Projective as Bls12G2Projective},
        bls12_381::{
            g1::{Bls12_381G1Parameters, G1Affine},
            Bls12_381, Bls12_381Parameters,
        },
        models::{ModelParameters, SWModelParameters},
        PairingCurve, PairingEngine,
    },
    fields::{
        bls12_381::{Fq, Fq12, Fq2, Fr},
        Field,
    },
};

pub type G2Affine = Bls12G2Affine<Bls12_381Parameters>;
pub type G2Projective = Bls12G2Projective<Bls12_381Parameters>;

impl PairingCurve for G2Affine {
    type Engine = Bls12_381;
    type Prepared = G2Prepared<Bls12_381Parameters>;
    type PairWith = G1Affine;
    type PairingResult = Fq12;

    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(*self)
    }

    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        Bls12_381::pairing(*other, *self)
    }
}

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bls12_381G2Parameters;

impl ModelParameters for Bls12_381G2Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

impl SWModelParameters for Bls12_381G2Parameters {
    /// COEFF_A = [0, 0]
    const COEFF_A: Fq2 = Fq2::new(
        Bls12_381G1Parameters::COEFF_A,
        Bls12_381G1Parameters::COEFF_A,
    );

    /// COEFF_B = [4, 4]
    const COEFF_B: Fq2 = Fq2::new(
        Bls12_381G1Parameters::COEFF_B,
        Bls12_381G1Parameters::COEFF_B,
    );

    /// COFACTOR = (x^8 - 4 x^7 + 5 x^6) - (4 x^4 + 6 x^3 - 4 x^2 - 4 x + 13) //
    /// 9
    /// = 305502333931268344200999753193121504214466019254188142667664032982267604182971884026507427359259977847832272839041616661285803823378372096355777062779109
    const COFACTOR: &'static [u64] = &[
        0xcf1c38e31c7238e5,
        0x1616ec6e786f0c70,
        0x21537e293a6691ae,
        0xa628f1cb4d9e82ef,
        0xa68a205b2e5a7ddf,
        0xcd91de4547085aba,
        0x91d50792876a202,
        0x5d543a95414e7f1,
    ];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// 26652489039290660355457965112010883481355318854675681319708643586776743290055
    const COFACTOR_INV: Fr = Fr::new(BigInteger256([
        6746407649509787816,
        1304054119431494378,
        2461312685643913071,
        5956596749362435284,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

pub const G2_GENERATOR_X: Fq2 = Fq2::new(G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
pub const G2_GENERATOR_Y: Fq2 = Fq2::new(G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

/// G2_GENERATOR_X_C0 =
/// 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160
pub const G2_GENERATOR_X_C0: Fq = Fq::new(BigInteger384([
    0xf5f28fa202940a10,
    0xb3f5fb2687b4961a,
    0xa1a893b53e2ae580,
    0x9894999d1a3caee9,
    0x6f67b7631863366b,
    0x58191924350bcd7,
]));

/// G2_GENERATOR_X_C1 =
/// 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758
pub const G2_GENERATOR_X_C1: Fq = Fq::new(BigInteger384([
    0xa5a9c0759e23f606,
    0xaaa0c59dbccd60c3,
    0x3bb17e18e2867806,
    0x1b1ab6cc8541b367,
    0xc2b6ed0ef2158547,
    0x11922a097360edf3,
]));

/// G2_GENERATOR_Y_C0 =
/// 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905
pub const G2_GENERATOR_Y_C0: Fq = Fq::new(BigInteger384([
    0x4c730af860494c4a,
    0x597cfa1f5e369c5a,
    0xe7e6856caa0a635a,
    0xbbefb5e96e0d495f,
    0x7d3a975f0ef25a2,
    0x83fd8e7e80dae5,
]));

/// G2_GENERATOR_Y_C1 =
/// 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582
pub const G2_GENERATOR_Y_C1: Fq = Fq::new(BigInteger384([
    0xadc0fc92df64b05d,
    0x18aa270a2b1461dc,
    0x86adac6a3be4eba0,
    0x79495c4ec93da33a,
    0xe7175850a43ccaed,
    0xb2bc2a163de1bf2,
]));
