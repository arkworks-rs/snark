use algebra::fields::tweedle::{
    Fq as TweedleFq,
    Fr as TweedleFr,
};
use primitives::crh::parameters::{
    TweedleFqPoseidonParameters,
    TweedleFrPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFrQuinticSbox
};
use crate::crh::{
    sbox::QuinticSBoxGadget,
    poseidon::{
        PoseidonHashGadget, PoseidonSpongeGadget,
    },
};

type TweedleFqQuinticSBoxGadget = QuinticSBoxGadget<TweedleFq, TweedleFqQuinticSbox>;
pub type TweedleFqPoseidonHashGadget = PoseidonHashGadget<
    TweedleFq,
    TweedleFqPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFqQuinticSBoxGadget,
>;

pub type TweedleFqPoseidonSpongeGadget = PoseidonSpongeGadget<
    TweedleFq,
    TweedleFqPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFqQuinticSBoxGadget,
>;

type TweedleFrQuinticSBoxGadget = QuinticSBoxGadget<TweedleFr, TweedleFrQuinticSbox>;
pub type TweedleFrPoseidonHashGadget = PoseidonHashGadget<
    TweedleFr,
    TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
    TweedleFrQuinticSBoxGadget,
>;

pub type TweedleFrPoseidonSpongeGadget = PoseidonSpongeGadget<
    TweedleFr,
    TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
    TweedleFrQuinticSBoxGadget,
>;