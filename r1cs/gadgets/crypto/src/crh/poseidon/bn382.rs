use algebra::fields::bn_382::{
    Fq as BN382Fq,
    Fr as BN382Fr,
};
use primitives::crh::parameters::{
    BN382FqPoseidonParameters,
    BN382FrPoseidonParameters,
    BN382FqQuinticSbox,
    BN382FrQuinticSbox
};
use crate::crh::{
    sbox::QuinticSBoxGadget,
    poseidon::PoseidonHashGadget,
};

type BN382FqQuinticSBoxGadget = QuinticSBoxGadget<BN382Fq, BN382FqQuinticSbox>;
pub type BN382FqPoseidonHashGadget = PoseidonHashGadget<
    BN382Fq,
    BN382FqPoseidonParameters,
    BN382FqQuinticSbox,
    BN382FqQuinticSBoxGadget,
>;

type BN382FrQuinticSBoxGadget = QuinticSBoxGadget<BN382Fr, BN382FrQuinticSbox>;
pub type BN382FrPoseidonHashGadget = PoseidonHashGadget<
    BN382Fr,
    BN382FrPoseidonParameters,
    BN382FrQuinticSbox,
    BN382FrQuinticSBoxGadget,
>;