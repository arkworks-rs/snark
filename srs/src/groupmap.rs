use crate::CompressedSRS;
use algebra::{
    curves::{models::SWModelParameters, AffineCurve},
    fields::{SquareRootField, Field},
};
use core::fmt::Error;

pub trait GroupMap {
    type F : SquareRootField;
    fn setup() -> Self;
    fn batch_to_group_x(p: &Self, ts: Vec<Self::F>) -> Vec<[Self::F; 3]>;
}

struct BWParameters<G: SWModelParameters> {
    u: G::BaseField,
    fu: G::BaseField,
    sqrt_neg_three_u_squared_plus_u_over_2: G::BaseField,
    sqrt_neg_three_u_squared_minus_u_over_2: G::BaseField,
    sqrt_neg_three_u_squared: G::BaseField,
    inv_three_u_squared: G::BaseField,
}

fn curve_eqn<G : SWModelParameters>(x : G::BaseField) -> G::BaseField {
    let mut res = x;
    res *= & x; // x^2
    res += & G::COEFF_A; // x^2 + A
    res *= & x; // x^3 + A x
    res += & G::COEFF_B; // x^3 + A x + B
    res
}

fn find_first<A, F: Fn(u64) -> Option<A>>(start: u64, f : F) -> A {
    for i in start.. {
        match f(i) {
            Some(x) => return x,
            None => ()
        }
    }
    panic!("find_first")
}

fn potential_xs_helper<G: SWModelParameters>(
    params: &BWParameters<G>,
    t2: G::BaseField,
    alpha: G::BaseField,
) -> [G::BaseField; 3] {
    let x1 = {
        let temp = t2;
        temp.square_in_place();
        temp *= &alpha;
        temp *= &params.sqrt_neg_three_u_squared;
        params.sqrt_neg_three_u_squared_minus_u_over_2 - &temp
    };

    let x2 = -params.u - & x1;

    let x3 = {
        let t2_plus_fu = t2 + &params.fu;
        let t2_inv = alpha * & t2_plus_fu ;
        let mut temp = t2_plus_fu.square();
        temp *= &t2_inv;
        temp *= &params.inv_three_u_squared;
        params.u - &temp
    };

    [x1, x2, x3]
}


fn potential_xs<G: SWModelParameters>(
    params: BWParameters<G>,
    t: G::BaseField,
) -> [G::BaseField; 3] {
    let t2 = t.square();
    let mut alpha_inv = t2;
    alpha_inv += &params.fu;
    alpha_inv *= &t2;

    let alpha = match alpha_inv.inverse() {
        Some(x) => x,
        None => G::BaseField::zero()
    };

    let x1 = {
        let temp = t2;
        temp.square_in_place();
        temp *=(&alpha);
        temp *= (&params.sqrt_neg_three_u_squared);
        params.sqrt_neg_three_u_squared_minus_u_over_2 - &temp
    };

    let x2 = -params.u - & x1;

    let x3 = {
        let t2_plus_fu = t2 + &params.fu;
        let t2_inv = alpha * & t2_plus_fu ;
        let mut temp = t2_plus_fu.square();
        temp *= &t2_inv;
        temp *= &params.inv_three_u_squared;
        params.u - &temp
    };

    [x1, x2, x3]
}

impl<G: SWModelParameters> GroupMap for BWParameters<G> {
    type F = G::BaseField;

    fn setup() -> Self {
        assert!(G::COEFF_A.is_zero());

        let (u, fu) = find_first(1, |u| {
            let u : G::BaseField = u.into();
            let fu : G::BaseField = curve_eqn::<G>(u);
            if fu.is_zero() {
                return None
            } else {
                Some((u, fu))
            }
        });

        let three_u_squared = u.square() * & 3.into();
        let inv_three_u_squared = three_u_squared.inverse().unwrap();
        let sqrt_neg_three_u_squared = (-three_u_squared).sqrt().unwrap();
        let two_inv =
            G::BaseField::from(2).inverse().unwrap();
        let sqrt_neg_three_u_squared_plus_u_over_2 = 
            (sqrt_neg_three_u_squared + & u) * &two_inv;
        let sqrt_neg_three_u_squared_minus_u_over_2 =
            (sqrt_neg_three_u_squared - & u) * &two_inv;

        BWParameters::<G> {
            u,
            fu,
            sqrt_neg_three_u_squared_plus_u_over_2,
            sqrt_neg_three_u_squared_minus_u_over_2,
            sqrt_neg_three_u_squared,
            inv_three_u_squared,
        }
    }

    fn batch_to_group_x(p: &BWParameters<G>, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]> {
        let t2_alpha_invs : Vec<_> = ts.iter().map(|t| {
            let t2 = t.square();
            let mut alpha_inv = t2;
            alpha_inv += &p.fu;
            alpha_inv *= &t2;
            (t2, alpha_inv)
        }).collect();

        let mut alphas : Vec<G::BaseField> = t2_alpha_invs.iter().map(|(_, a)| a.clone()).collect();
        algebra::fields::batch_inversion::<G::BaseField>(&mut alphas);

        let potential_xs = t2_alpha_invs.iter().zip(alphas).map(|((t2,_), alpha)| {
            potential_xs_helper(p, t2.clone(), alpha.clone())
        });
    }

    /*
    fn batch_to_group_x(p: &Self::Parameters, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]> {
        fn get_y(x: G::BaseField) -> Option<G::BaseField> {
            let fx = curve_eqn::<G>(x);
            // how do we choose +/- sqrt?
            if let Some(y) = fx.sqrt() {
                Some(y)
            } else {
                None
            }
        }

        fn get_xyi(t: G::BaseField) -> Option<G::BaseField, G::BaseField, usize> {
            let xvec = &Self::potential_xs(t);
            let yvec = xvec.into_iter().map(|x| &Self::get_y(x));
            let i = yvec.into_iter().find(|y| y.is_some());
            if let Some(i) = b {
                Some(xvec[i], yvec[i], i)
            } else {
                None
            }
        }

        fn to_group(t: G::BaseField) -> Option<G> {
            if let Some(x, y, i) = get_xyi(t) {
                Some(G { x, y })
            } else {
                None
            }
        }

        fn to_group_bit(t: G::BaseField) -> Option<usize, G::BaseField> {
            if let Some(x, y, i) = get_xyi(t) {
                Some(i, y)
            } else {
                None
            }
        }
    }
    */
}
