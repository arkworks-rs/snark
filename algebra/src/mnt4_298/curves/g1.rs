use crate::mnt4_298::{self, Fq, Fr, FR_ONE};
use algebra_core::{
    biginteger::BigInteger320,
    curves::{
        mnt4,
        models::{ModelParameters, SWModelParameters},
    },
    field_new,
};

pub type G1Affine = mnt4::G1Affine<mnt4_298::Parameters>;
pub type G1Projective = mnt4::G1Projective<mnt4_298::Parameters>;
pub type G1Prepared = mnt4::G1Prepared<mnt4_298::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for Parameters {
    /// COEFF_A = 2
    /// Reference: https://github.com/scipr-lab/libff/blob/c927821ebe02e0a24b5e0f9170cec5e211a35f08/libff/algebra/curves/mnt/mnt4/mnt4_init.cpp#L116
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger320([
            3568597988870129848,
            15257338106490985450,
            10069779447956199041,
            5922375556522222383,
            3858029504390,
    ]));

    /// COEFF_B = 423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685
    /// Reference: https://github.com/scipr-lab/libff/blob/c927821ebe02e0a24b5e0f9170cec5e211a35f08/libff/algebra/curves/mnt/mnt4/mnt4_init.cpp#L117
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger320([
            7842808090366692145,
            288200302308193399,
            4162060950790347941,
            5488589108190218591,
            1553456013645,
    ]));

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[1];

    /// COFACTOR^(-1) mod r =
    /// 1
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = FR_ONE;

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);
}

// Generator of G1
// X = 60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838,
// Y = 363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306,
/// G1_GENERATOR_X
/// Reference: https://github.com/scipr-lab/libff/blob/c927821ebe02e0a24b5e0f9170cec5e211a35f08/libff/algebra/curves/mnt/mnt4/mnt4_init.cpp#L137
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger320([
    6046301378120906932,
    15105298306031900263,
    15757949605695610691,
    6113949277267426050,
    3063081829217,
]));

/// G1_GENERATOR_Y
/// Reference: https://github.com/scipr-lab/libff/blob/c927821ebe02e0a24b5e0f9170cec5e211a35f08/libff/algebra/curves/mnt/mnt4/mnt4_init.cpp#L138
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger320([
    8798367863963590781,
    9770379341721339603,
    17697354471293810920,
    15252694996423733496,
    3845520398052,
]));
