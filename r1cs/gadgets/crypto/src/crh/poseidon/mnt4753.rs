use algebra::fields::mnt4753::Fr as MNT4753Fr;
use primitives::crh::parameters::MNT4753PoseidonParameters;
use crate::crh::poseidon::PoseidonHashGadget;

pub type MNT4PoseidonHashGadget = PoseidonHashGadget<MNT4753Fr, MNT4753PoseidonParameters>;
