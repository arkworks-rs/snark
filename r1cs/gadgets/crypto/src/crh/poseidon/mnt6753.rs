use algebra::fields::mnt6753::Fr as MNT6753Fr;
use primitives::crh::parameters::MNT6753PoseidonParameters;
use crate::crh::poseidon::PoseidonHashGadget;

pub type MNT6PoseidonHashGadget = PoseidonHashGadget<MNT6753Fr, MNT6753PoseidonParameters>;
