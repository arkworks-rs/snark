use algebra::curves::AffineCurve;
use core::fmt::Error;

#[derive(Copy, Clone, Debug)]
pub struct CompressedSRS<G: AffineCurve> {
    pub y_coords: Vec<G::BaseField>,
    pub x_bits: usize, // could be [bool; 2] ?
}

#[derive(Copy, Clone, Debug)]
struct SRSParameters<G: AffineCurve> {
    pub u: G::BaseField,
    pub u_over_2: G::BaseField,
    pub projp: G,
    pub conic_c: G::BaseField,
    pub a: G::BaseField,
    pub b: G::BaseField,
}

// C(z, y) : z^2 + 3u/4 y^2 = -f(u)
#[derive(Copy, Clone, Debug)]
struct SRSConic<G: AffineCurve> {
    pub z: G::BaseField,
    pub y: G::BaseField,
}

// S(u, v, y) : y^2(u^2 + uv + v^2 + a) = -f(u)
#[derive(Copy, Clone, Debug)]
struct SRSSurface<G: AffineCurve> {
    pub u: G::BaseField,
    pub v: G::BaseField,
    pub y: G::BaseField,
}

// do we need this
// V(x1, x2, x3, x4) : f(x1) f(x2) f(x3) = x4^2
#[derive(Copy, Clone, Debug)]
struct SRSVariety<G: AffineCurve> {
    pub x1: G::BaseField,
    pub x2: G::BaseField,
    pub x3: G::BaseField,
}

pub trait CompressedSRS {
    type Parameters: Copy, Clone;
    type Conic: Copy, Clone;
    type Surface: Copy, Clone;
    type Variety: Copy, Clone;
    type Parameters: Copy, Clone;
    type F: Copy, Clone;

    fn setup() -> Result<Self::Parameters, Error>;

    fn randomize_signature(
        pp: &Self::Parameters,
        signature: &Self::Signature,
        randomness: &[u8],
    ) -> Result<Self::Signature, Error>;
}

impl<G: AffineCurve> CompressedSRS<G> {
    type Parameters = SRSParameters<G>;
    type Conic = SRSConic<G>;
    type Surface = SRSSurface<G>;
    type Variety = SRSVariety<G>;
    type F = G::BaseField;

    // let first_map f =
    // let rec go i = match f i with Some x -> x | None -> go (i + one) in
    // go zero
    // in
    // let first f = first_map (fun x -> Option.some_if (f x) x) in
    // let three_fourths = of_int 3 / of_int 4 in
    // let curve_eqn u = (u * u * u) + (a * u) + b in
    // let u =
    // first (fun u ->
    // (* from (15), A = 0, B = Params.a *)
    // let check = (three_fourths * u * u) + a in
    // let fu = curve_eqn u in
    // (not (equal check zero))
    // && (not (equal fu zero))
    // && not (is_square (negate fu))
    // )
    // in
    // let conic_c = (three_fourths * u * u) + a in
    // let conic_d = negate (curve_eqn u) in
    // let projp =
    // first_map (fun y ->
    // let z2 = conic_d - (conic_c * y * y) in
    // if F.is_square z2 then Some {Conic.z= F.sqrt z2; y} else None )
    // in
    // {u; u_over_2= u / of_int 2; conic_c; projp; a; b}
    //
    // is this not found somewhere else?
    fn curve_eqn(parameters: &Self::Parameters, x: &Self::F) -> Self::G {
        Self::G {
            x,
            y: (x * x * x) + (parameters.a * x) + parameters.b,
        }
    }

    fn setup() -> Result<Self::Parameters, Error> {}

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

    fn to_group(t: &Self::F) -> Option<Self::G> {
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

    fn create(length: usize, num_y_coords: usize) -> Self {}

    fn decompress(c: Self) -> Vec<(Self::F, Self::F)> {}
}

#[cfg(test)]
mod tests {

    use super;
    #[test]
    fn test_projection_point_well_formed() {
        assert!(on_conic(parameters.projp));
    }

    #[test]
    fn test_field_to_conic() {
        assert!(on_conic(field_to_conic(t)));
    }

    #[test]
    fn test_conic_to_s() {
        let z2 = conic_d - (parameters.conic_c * y * y);
        if is_square(z2) {
            let c = Some(Conic { z: sqrt(z2), y });
        } else {
            let c = None;
        }
        assert!(on_s(conic_to_s(p)));
    }

    #[test]
    fn test_field_to_s() {
        ct = M.field_to_conic(t);
        st = M.conic_to_s(ct);
        assert!(on_s(st));
    }

    #[test]
    fn test_field_to_V() {
        let s = M.conic_to_s(M.field_to_conic(t));
        assert!(on_v(M._s_to_v(s)));
    }

    #[test]
    fn test_full_map() {
        let (x, y) = to_group(parameters, t);
        assert_eq!(curve_eqn(x), y * y);
    }
}
