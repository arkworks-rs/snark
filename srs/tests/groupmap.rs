use srs::groupmap::{GroupMap, BWParameters};
use algebra::{
    curves::{
        bn_382::G1Affine,
        bn_382::g1::Bn_382G1Parameters,
        models::SWModelParameters
    },
    fields::bn_382::Fq,
};
use rand;

type G = Bn_382G1Parameters;

#[test]
fn test_group_map_on_curve() {
    let params =  BWParameters::<G>::setup();
    let t : Fq = rand::random();
    let (x, y) = BWParameters::<G>::to_group(&params, t);
    let g = G1Affine::new(x, y, false);
    assert!(g.is_on_curve());
}
