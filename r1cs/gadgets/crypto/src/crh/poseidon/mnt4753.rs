use algebra::fields::mnt4753::Fr as MNT4753Fr;
use primitives::crh::parameters::{
    MNT4753PoseidonParameters, MNT4InversePoseidonSBox,
};
use crate::crh::{
    sbox::InverseSBoxGadget,
    poseidon::PoseidonHashGadget,
};

type MNT4InverseSBoxGadget = InverseSBoxGadget<MNT4753Fr, MNT4InversePoseidonSBox>;
pub type MNT4PoseidonHashGadget = PoseidonHashGadget<
    MNT4753Fr,
    MNT4753PoseidonParameters,
    MNT4InversePoseidonSBox,
    MNT4InverseSBoxGadget
>;