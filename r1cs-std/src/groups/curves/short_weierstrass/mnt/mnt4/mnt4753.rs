use algebra::curves::mnt4753::MNT4_753Parameters;
use super::{G1Gadget, G2Gadget, G1PreparedGadget, G2PreparedGadget};

pub type MNT4G1Gadget = G1Gadget<MNT4_753Parameters>;
pub type MNT4G2Gadget = G2Gadget<MNT4_753Parameters>;

pub type MNT4G1PreparedGadget = G1PreparedGadget<MNT4_753Parameters>;
pub type MNT4G2PreparedGadget = G2PreparedGadget<MNT4_753Parameters>;