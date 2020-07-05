use crate::mnt4_298::{Fq, Fq2, Fq2Parameters, FQ_ONE, FQ_ZERO};
use algebra_core::{
    biginteger::BigInteger320 as BigInteger,
    field_new,
    fields::fp4::{Fp4, Fp4Parameters},
};

pub type Fq4 = Fp4<Fq4Parameters>;

pub struct Fq4Parameters;

impl Fp4Parameters for Fq4Parameters {
    type Fp2Params = Fq2Parameters;

    const NONRESIDUE: Fq2 = field_new!(Fq2, FQ_ZERO, FQ_ONE);

    // Coefficients for the Frobenius automorphism.
    // c1[0] = 1,
    // c1[1] = 7684163245453501615621351552473337069301082060976805004625011694147890954040864167002308
    // c1[2] = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758080
    // c1[3] = 468238122923807824137727898100575114475823797181717920390930116882062371863914936316755773
    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP4_C1: &'static [Fq] = &[
        FQ_ONE,
        field_new!(
            Fq,
            BigInteger([
                16439849825752526567,
                14772594681319164557,
                16175669228740845684,
                4590896976404796446,
                3810243174413
            ])
        ),
        field_new!(
            Fq,
            BigInteger([
                12702890790846888869,
                6326265861366186013,
                364584707886187945,
                8740893163049517815,
                2181130330288
            ])
        ),
        field_new!(
            Fq,
            BigInteger([
                16494084033238978842,
                8405712270147289988,
                16893921313687769205,
                7111183964905832559,
                299901908070
            ])
        ),
    ];
}
