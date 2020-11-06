extern crate algebra;
extern crate algebra_core;

use algebra::bw6_761::G1Projective as GroupProjective;
use algebra_core::{BigInteger768 as BaseFieldBigInt, BigInteger768 as FrWideBigInt};
use glv_lattice_basis::*;

fn main() {
    print_glv_params::<GroupProjective, FrWideBigInt, BaseFieldBigInt>();
}
