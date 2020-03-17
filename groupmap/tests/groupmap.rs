use groupmap::{GroupMap, BWParameters};
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

fn first_xy(xs: & [Fq; 3]) -> (Fq, Fq) {
    for x in xs.iter() {
        match groupmap::get_y::<G>(*x) {
            Some(y) => return (*x, y),
            None => ()
        }
    }
    panic!("get_xy")
}

#[test]
fn test_batch_group_map_on_curve() {
    let n = 1000;
    let params =  BWParameters::<G>::setup();
    let ts : Vec<Fq> = (0..1000).map(|_| rand::random()).collect();
    for xs in BWParameters::<G>::batch_to_group_x(&params, ts).iter() {
        let (x, y) = first_xy(xs);
        let g = G1Affine::new(x, y, false);
        assert!(g.is_on_curve());
    };
}
