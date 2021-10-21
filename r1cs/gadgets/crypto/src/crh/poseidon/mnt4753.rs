use crate::crh::{poseidon::PoseidonHashGadget, sbox::InverseSBoxGadget};
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use primitives::crh::parameters::{MNT4753PoseidonParameters, MNT4InversePoseidonSBox};

type MNT4InverseSBoxGadget = InverseSBoxGadget<MNT4753Fr, MNT4InversePoseidonSBox>;
pub type MNT4PoseidonHashGadget = PoseidonHashGadget<
    MNT4753Fr,
    MNT4753PoseidonParameters,
    MNT4InversePoseidonSBox,
    MNT4InverseSBoxGadget,
>;
