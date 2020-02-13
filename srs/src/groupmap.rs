use crate::CompressedSRS;
use algebra::curves::{AffineCurve, PrimeField};
use core::fmt::Error;

#[derive(Copy, Clone, Debug)]
pub struct BN382CompressedSRS<G: AffineCurve> {
    pub y_coords: Vec<G::BaseField>,
    pub x_bits: usize, // could be [bool; 2] ?
}

#[derive(Copy, Clone, Debug)]
struct BWParameters<G: AffineCurve> {
    pub u: G::BaseField,
    pub u_over_2: G::BaseField,
    pub projp: G,
    pub conic_c: G::BaseField,
    pub a: G::BaseField,
    pub b: G::BaseField,
}

// C(z, y) : z^2 + 3u/4 y^2 = -f(u)
#[derive(Copy, Clone, Debug)]
struct BWConic<G: AffineCurve> {
    pub z: G::BaseField,
    pub y: G::BaseField,
}

// S(u, v, y) : y^2(u^2 + uv + v^2 + a) = -f(u)
#[derive(Copy, Clone, Debug)]
struct BWSurface<G: AffineCurve> {
    pub u: G::BaseField,
    pub v: G::BaseField,
    pub y: G::BaseField,
}

// do we need this
// V(x1, x2, x3, x4) : f(x1) f(x2) f(x3) = x4^2
#[derive(Copy, Clone, Debug)]
struct BWVariety<G: AffineCurve> {
    pub x1: G::BaseField,
    pub x2: G::BaseField,
    pub x3: G::BaseField,
}


impl groupmap<G: AffineCurve> {
    type Parameters = BWParameters<G>;
    type Conic = BWConic<G>;
    type Surface = BWSurface<G>;
    type Variety = BWVariety<G>;
    type F: G::BaseField + Into<G::BaseField as PrimeField>;

    fn setup() -> Result<Self::Parameters, Error> {
        let three_quarters = F::from_int(3) / F::from_int(4);
        // TODO get formula for u
        let u = three_quarters;
        let check = (three_quarters * u * u) + G::A_COEFF;
        let fu = Self::curve_eqn(u);
        assert_ne!(check, F::zero());
        assert_ne!(fu, F::zero());
        // TODO? CHECK
        assert!(!F::is_square(F::negate(fu)));
        let conic_c = (three_quarters * u * u) + G::A_COEFF;
        let conic_d = F::negate(Self::curve_eqn(u));
        // TODO get forumla for y
        let y = three_quarters;
        let z2 = conic_d - (conic_c * y * y);
        // TODO Fix? this should never return None, it should return an error
        let projp = if F::is_square(z2) {
            Self::Conic { z: F::sqrt(z2), y }
        } else {
            None
        };
        Self::Parameters {
            u,
            u_over_2: u / F::of_int(2),
            conic_c,
            projp,
            a: G::A_COEFF,
            b: G::B_COEFF,
        }
    }
    //
    // is this not found somewhere else?
    fn curve_eqn(parameters: &Self::Parameters, x: &Self::F) -> Self::G {
        Self::G {
            x,
            y: (x * x * x) + (parameters.a * x) + parameters.b,
        }
    }

    // For a curve z^2 + c y^2 = d and a point (z0, y0) on the curve, there
    // is one other point on the curve which is also on the line through (z0, y0)
    // with slope t. This function returns that point.
    fn field_to_conic(parameters: &Self::Parameters, t: &Self::F) -> Result<Self::Conic, Error> {
        let (z0, y0) = (parameters.projp.z, parameters.projp.y);
        let ct = parameters.conic_c * t;
        let s = 2 * ((ct * y0) + z0) / ((ct * t) + &Self::F::One);
        &Self::Conic {
            z: z0 - s,
            y: y0 - (s * t),
        }
    }

    // From (16) : φ(λ) : F → S : λ → ( u, α(λ)/β(λ) - u/2, β(λ) )
    // what are alpha(lambda) and beta(lambda) ?
    fn conic_to_s(parameters: &Self::Parameters, c: &Self::Conic) -> Result<Self::Surface, Error> {
        &Self::Surface {
            u: parameters.u,
            v: (c.z / c.y) - parameters.u_over_2,
            y: c.y,
        }
    }

    // We don't need to compute the final coordinate in V
    fn s_to_v_truncated(S: &Self::Surface) -> Result<Self::Variety, Error> {
        &Self::Variety {
            x1: S.v,
            x2: (S.u + S.v).negate(),
            x3: S.u + (S.y * S.y),
        }
    }

    fn potential_xs(x: &Self::F) -> Result<Self::Surface, Error> {
        &Self::s_to_v_truncated(&Self::conic_to_s(&Self::field_to_conic(x)));
    }

    fn try_decode(x: &Self::F) -> Option<Self::G> {
        let yy = &Self::curve_eqn(x);
        if &Self::is_square(yy) {
            Some(Self::G {
                x,
                y: Self::sqrt(yy),
            })
        } else {
            None
        }
    }

    fn to_group(t: &Self::F) -> Option<G> {
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
