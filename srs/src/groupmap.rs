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

impl<G: SWModelParameters> GroupMap for BWParameters<G> {
    type F = G::BaseField;

    fn setup() -> Self {
        assert!(G::COEFF_A.is_zero());

        let (u, fu) = find_first(1, |u| {
            let u : G::BaseField = u.into();
            let fu : G::BaseField = curve_eqn(u);
            if fu.is_zero() {
                return None
            } else {
                Some((u, fu))
            }
        });

        let three_u_squared = u.square() * & 3.into();
        let inv_three_u_squared = three_u_squared.inverse();
        let sqrt_neg_three_u_squared = (-three_u_squared).sqrt().unwrap();

        let three_u = &u.square();
        three_u.mul(G::BaseField::from(3));
        three_u.minus();
        let inv_three_u = three_u.copy();
        three_u.sqrt(); // sqrt(-3u^2)
        inv_three_u.inv(); // 1/(-3u^2)
        let three_u_plus_u_over_2 = three_u.copy();
        let three_u_minus_u_over_2 = three_u.copy();
        three_u_plus_u_over_2.add(&u);
        three_u_minus_u_over_2.sub(&u);
        three_u_plus_u_over_2.div(G::BaseField::from(2)); // sqrt(-3u^2) + u / 2
        three_u_minus_u_over_2.div(G::BaseField::from(2)); // sqrt(-3u^2) - u / 2

        Self::Parameters {
            u,
            fu,
            three_u_plus_u_over_2,
            three_u_minus_u_over_2,
            three_u,
            inv_three_u,
        }
    }

    /*
    fn batch_to_group_x(p: &Self::Parameters, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]> {
        fn get_y(x: G::BaseField) -> Option<G::BaseField> {
            let fx = ((x * x * x) + (G::A_COEFF * x) + G::B_COEFF);
            // how do we choose +/- sqrt?
            if let Some(y) = fx.sqrt() {
                Some(y)
            } else {
                None
            }
        }

        fn potential_xs(
            params: Self::BWParameters,
            t: G::BaseField,
        ) -> Result<Vec<[G::BaseField; 3]>, Error> {
            let alpha_inv = &t.square();
            alpha_inv.add(&params.fu);
            alpha_inv.mul(&t.square());
            // check alpha_inv non zero
            let alpha = alpha_inv.copy();
            alpha.inv(); // 1/alpha_inv = alpha

            let temp = &t.square();
            temp.square();
            temp.mul(&alpha);
            temp.add(&params.three_u);
            let x1 = params.three_u_minus_u_over_2.copy();
            x1.sub(&temp);
            let x2 = F::zero().sub(&params.u);
            x2.sub(&x1);
            alpha_inv.square();
            temp = &t.square();
            temp.square();
            temp.inv();
            temp.mul(&params.inv_three_u);
            temp.mul(&alpha_inv);
            let x3 = params.u.copy();
            x3.sub(&temp);

            vec![x1, x2, x3]
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
