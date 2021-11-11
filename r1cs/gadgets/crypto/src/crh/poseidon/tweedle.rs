use crate::crh::{poseidon::PoseidonHashGadget, sbox::QuinticSBoxGadget};
use algebra::fields::tweedle::{Fq, Fr};
use primitives::crh::parameters::{
    TweedleFqPoseidonParameters, TweedleFqQuinticSbox, TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
};

type TweedleFqQuinticSboxGadget = QuinticSBoxGadget<Fq, TweedleFqQuinticSbox>;
pub type TweedleFqPoseidonHashGadget = PoseidonHashGadget<
    Fq,
    TweedleFqPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFqQuinticSboxGadget,
>;

type TweedleFrQuinticSboxGadget = QuinticSBoxGadget<Fr, TweedleFrQuinticSbox>;
pub type TweedleFrPoseidonHashGadget = PoseidonHashGadget<
    Fr,
    TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
    TweedleFrQuinticSboxGadget,
>;
