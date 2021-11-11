use crate::crh::{poseidon::PoseidonHashGadget, sbox::InverseSBoxGadget};
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use primitives::crh::parameters::{MNT6753PoseidonParameters, MNT6InversePoseidonSBox};

type MNT6InverseSBoxGadget = InverseSBoxGadget<MNT6753Fr, MNT6InversePoseidonSBox>;
pub type MNT6PoseidonHashGadget = PoseidonHashGadget<
    MNT6753Fr,
    MNT6753PoseidonParameters,
    MNT6InversePoseidonSBox,
    MNT6InverseSBoxGadget,
>;
