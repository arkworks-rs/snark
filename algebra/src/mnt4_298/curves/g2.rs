use crate::mnt4_298::{self, Fq, Fq2, Fr, FQ_ZERO, G1_COEFF_A_NON_RESIDUE};
use algebra_core::{
    biginteger::BigInteger320,
    curves::{
        mnt4,
        mnt4::MNT4Parameters,
        models::{ModelParameters, SWModelParameters},
    },
    field_new,
};

pub type G2Affine = mnt4::G2Affine<mnt4_298::Parameters>;
pub type G2Projective = mnt4::G2Projective<mnt4_298::Parameters>;
pub type G2Prepared = mnt4::G2Prepared<mnt4_298::Parameters>;

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
    const COEFF_A: Fq2 = mnt4_298::Parameters::TWIST_COEFF_A;
    // B coefficient of MNT4-298 G2 =
    // ```
    // mnt4298_twist_coeff_b = mnt4298_Fq2(mnt4298_Fq::zero(),
    //                                  mnt4298_G1::coeff_b * mnt4298_Fq2::non_residue);
    // non_residue = mnt4298_Fq2::non_residue = mnt4298_Fq("13");
    //  = (ZERO, G1_B_COEFF * NON_RESIDUE);
    //  =
    //  (0, 67372828414711144619833451280373307321534573815811166723479321465776723059456513877937430)
    // ```
    #[rustfmt::skip]
    const COEFF_B: Fq2 = field_new!(Fq2,
        FQ_ZERO,
        field_new!(Fq, BigInteger320([
            9511110677122940475,
            13403516020116973437,
            1464701424831086967,
            4646785117660390394,
            1747881737068,
        ])),
    );

    /// COFACTOR =
    /// 475922286169261325753349249653048451545124879932565935237842521413255878328503110407553025
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        15480692783052488705,
        9802782456999489873,
        14622846468721090623,
        11702080941310629006,
        4110145082483,
    ];

    /// COFACTOR^(-1) mod r =
    /// 475922286169261325753349249653048451545124878207887910632124039320641839552134835598065665
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger320([
        8065818351154103109,
        7537800592537321232,
        747075088561892445,
        6335802185495034136,
        1874289794052,
    ]));

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
// X = 438374926219350099854919100077809681842783509163790991847867546339851681564223481322252708,
// 37620953615500480110935514360923278605464476459712393277679280819942849043649216370485641,
// Y = 37437409008528968268352521034936931842973546441370663118543015118291998305624025037512482,
// 424621479598893882672393190337420680597584695892317197646113820787463109735345923009077489,
#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger320([
    5356671649366391794,
    2684151262065976452,
    4683110650642896126,
    10421299515941681582,
    1618695480960
]));

#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger320([
    133394645290266480,
    15395232932057272770,
    18271324022738539173,
    9095178119640120034,
    2303787573609
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger320([
    16920448081812496532,
    15580160192086626100,
    3974467672100342742,
    8216505962266760277,
    2643162835232
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger320([
    73816197493558356,
    8663991890578965996,
    11575903875707445958,
    17953546933481201011,
    2167465829200
]));
