use crate::{
    bytes::ToBytes,
    curves::{
        bn::BnParameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
};
use derivative::Derivative;
use std::io::{Result as IoResult, Write};

pub type G1Affine<P> = GroupAffine<<P as BnParameters>::G1Parameters>;

pub type G1Projective<P> = GroupProjective<<P as BnParameters>::G1Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BnParameters"),
    Debug(bound = "P: BnParameters"),
    PartialEq(bound = "P: BnParameters"),
    Eq(bound = "P: BnParameters")
)]
pub struct G1Prepared<P: BnParameters>(pub G1Affine<P>);

impl<P: BnParameters> G1Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    pub fn from_affine(p: G1Affine<P>) -> Self {
        G1Prepared(p)
    }
}

impl<P: BnParameters> Default for G1Prepared<P> {
    fn default() -> Self {
        G1Prepared(G1Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: BnParameters> ToBytes for G1Prepared<P> {
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}
