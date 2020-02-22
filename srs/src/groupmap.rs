use crate::CompressedSRS;
use algebra::curves::{AffineCurve, PrimeField};
use core::fmt::Error;

#[derive(Copy, Clone, Debug)]
pub struct CompressedSRS<G: AffineCurve> {
    pub y_coords: Vec<G::BaseField>,
    pub x_bits: usize, // could be [bool; 2] ?
}

pub trait GroupMap<Params, G: AffineCurve> {
    fn setup() -> Params;

    fn batch_to_group_x(p: Params, ts: Vec<G::BaseField>) -> Vec<[G::BaseField; 3]>;
}

impl GroupMap<Params, G: AffineCurve> {

    struct BWParameters<G> {
        u: G::BaseField,
        fu: G::BaseField,
        three_u_plus_u_over_2: G::BaseField,
        three_u_minus_u_over_2: G::BaseField,
        three_u: G::BaseField,
        inv_three_u: G::BaseField
    }

    fn get_u() -> G::BaseField {
        for u in 0..G::P {
            let x_zero = &u.square();
            x_zero.mul(F::from_int(-3));
            x_zero.sqrt();
            x_zero.minus(&u);
            x_zero.div(F::from_int(2));
            let (x, y) = curve_eqn(x_zero);
            if legendre(y) == 1 {
                u
            }
        }
    }

    // should this be generic? can we get rid of A_COEFF?
    fn setup() -> Result<Self::BWParameters, Error> {
        let u = get_u();
        let fu = &u.mul(&u);
        fu.add(&G::A_COEFF);
        fu.mul(&u);
        fu.add(&G::B_COEFF);
        let three_u = &u.square();
        three_u.mul(F::from_int(-3));
        let inv_three_u = three_u.copy();
        three_u.sqrt(); // sqrt(-3u^2)
        inv_three_u.inv(); // 1/(-3u^2)
        let three_u_plus_u_over_2 = three_u.copy();
        let three_u_minus_u_over_2 = three_u.copy();
        three_u_plus_u_over_2.add(&u);
        three_u_minus_u_over_2.sub(&u);
        three_u_plus_u_over_2.div(F::from_int(2)); // sqrt(-3u^2) + u / 2
        three_u_minus_u_over_2.div(F::from_int(2)); // sqrt(-3u^2) - u / 2

        Self::BWParameters {
            u,
            fu,
            three_u_plus_u_over_2,
            three_u_minus_u_over_2,
            three_u,
            inv_three_u
        }

   }

    fn curve_eqn(x: G::BaseField) -> G {
        G {
            x,
            y: ((x * x * x) + (G::A_COEFF * x) + G::B_COEFF).sqrt();
        }
    }

    fn potential_xs(params : Self::BWParameters, t: G::BaseField) -> Result<Self::BWSurface, Error> {
        let alpha_inv = &t.square();
        alpha_inv.add(&params.fu);
        alpha_inv.mul(&t.square());
        // check alpha_inv non zero
        let alpha = alpha_inv.copy();
        alpha.inv(); // 1/alpha_inv = alpha

        let temp = &t.square();
        temp.square();
        temp.mul(&alpha);
        temp.(&params.three_u);
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

        BWSurface {
            x1,
            x2,
            x3
        }

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
