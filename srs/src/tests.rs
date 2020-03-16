use crate::groupmap::{GroupMap, BWParameters};
use algebra::{
    curves::{
        bn_382::g1::Bn_382G1Parameters,
        models::SWModelParameters
    },
    fields::bn_382::Fq,
};
use rand;

type BNGroupMap = dyn GroupMap::<F = algebra::fields::bn_382::Fp>;

pub fn test_group_map_on_curve<G: SWModelParameters>() {
    let params =  BWParameters::<G>::setup();
    let t : Fq = rand::random();
    let g = BWParameters::<G>::to_group(&params, t);
    assert!(g.is_on_curve());
}

pub fn test_group_map<G: SWModelParameters>() {
}
