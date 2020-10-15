use crate::pairing::bn::PairingGadget as BnPG;
use algebra::curves::bn_382::Bn382Parameters;

pub type PairingGadget = BnPG<Bn382Parameters>;