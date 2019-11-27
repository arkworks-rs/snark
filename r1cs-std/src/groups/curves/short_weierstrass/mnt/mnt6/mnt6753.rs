use algebra::curves::mnt6753::MNT6_753Parameters;

use super::{G1Gadget, G2Gadget, G1PreparedGadget, G2PreparedGadget};

pub type MNT6G1Gadget = G1Gadget<MNT6_753Parameters>;
pub type MNT6G2Gadget = G2Gadget<MNT6_753Parameters>;

pub type MNT6G1PreparedGadget = G1PreparedGadget<MNT6_753Parameters>;
pub type MNT6G2PreparedGadget = G2PreparedGadget<MNT6_753Parameters>;