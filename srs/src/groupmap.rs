use crate::CompressedSRS;
use algebra::{
    curves::{models::SWModelParameters, AffineCurve},
    fields::Field,
};
use core::fmt::Error;

pub trait GroupMap<G: SWModelParameters> {
    type Parameters;
    fn setup() -> Self::Parameters;
    fn batch_to_group_x(p: &Self::Parameters, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]>;
}

struct BWParameters<G: SWModelParameters> {
    u: G::BaseField,
    fu: G::BaseField,
    three_u_plus_u_over_2: G::BaseField,
    three_u_minus_u_over_2: G::BaseField,
    three_u: G::BaseField,
    inv_three_u: G::BaseField,
}

impl<G: SWModelParameters> GroupMap<G> for BWParameters<G> {
    type Parameters = BWParameters<G>;
    type FF = G::BaseField;

    fn setup() -> Self::Parameters {
        fn get_u<G: SWModelParameters>() -> Option<F> {
            let mut u: u8;
            for u in 0..200 {
                let x_zero = (<Self::FF as Field>::from(u)).square();
                x_zero.mul(G::BaseField::from(3));
                x_zero.minus();
                x_zero.sqrt();
                x_zero.ub(&u);
                x_zero.div(G::BaseField::from(2));
                let (x, y) = curve_eqn(x_zero);
                if y.legendre() == 1 {
                    Some(G::BaseField::from(u))
                };
            }
            None
        }

        let u = get_u();
        let fu = &u.mul(&u);
        fu.add(&G::COEFF_A);
        fu.mul(&u);
        fu.add(&G::COEFF_B);
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

    fn batch_to_group_x(p: &Self::Parameters, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]> {
        fn curve_eqn(x: G::BaseField) -> G {
            G {
                x,
                y: ((x * x * x) + (G::A_COEFF * x) + G::B_COEFF).sqrt(),
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

            BWSurface { x1, x2, x3 }
        }

        fn to_group(t: G::BaseField) -> Option<G> {
            let (x1, x2, x3) = &Self::potential_xs(t);
            if let Some(b1) = &Self::try_decode(x1) {
                Some(b1)
            } else if let Some(b2) = &Self::try_decode(x2) {
                Some(b2)
            } else if let Some(b3) = &Self::try_decode(x3) {
                Some(b3)
            } else {
                None
            }
        }
    }
}
