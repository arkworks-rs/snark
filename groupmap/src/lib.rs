use algebra::{
    curves::{models::SWModelParameters, AffineCurve},
    fields::{SquareRootField, Field},
};

pub trait GroupMap {
    type F : SquareRootField;
    fn setup() -> Self;
    fn to_group(p: &Self, u: Self::F) -> (Self::F, Self::F);
    fn batch_to_group_x(p: &Self, ts: Vec<Self::F>) -> Vec<[Self::F; 3]>;
}

pub struct BWParameters<G: SWModelParameters> {
    u: G::BaseField,
    fu: G::BaseField,
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
        let mut temp = t2;
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
    params: &BWParameters<G>,
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

    potential_xs_helper(params, t2, alpha)
}

pub fn get_y<G:SWModelParameters>(x: G::BaseField) -> Option<G::BaseField> {
    let fx = curve_eqn::<G>(x);
    if let Some(y) = fx.sqrt() {
        Some(y)
    } else {
        None
    }
}

fn get_xy<G:SWModelParameters>(params: &BWParameters<G>, t: G::BaseField) -> (G::BaseField, G::BaseField) {
    let xvec = potential_xs(&params, t);
    for x in xvec.iter() {
        match get_y::<G>(*x) {
            Some(y) => return (*x, y),
            None => ()
        }
    }
    panic!("get_xy")
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
        let sqrt_neg_three_u_squared_minus_u_over_2 =
            (sqrt_neg_three_u_squared - & u) * &two_inv;

        BWParameters::<G> {
            u,
            fu,
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
        potential_xs.collect()
    }

    fn to_group(p: &BWParameters<G>, t: G::BaseField) -> (G::BaseField, G::BaseField) {
        get_xy(p, t)
    }
}
