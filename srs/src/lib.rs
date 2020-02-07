use algebra;

#[derive(Copy, Clone, Debug)]
struct CompressedSRS<G : AffineCurve> {
  pub y_coords : Vec<G::BaseField>,
  pub x_bits : usize // could be [bool; 2] ?
}

#[derive(Copy, Clone, Debug)]
struct SRSParams<G : AffineCurve> {
    pub u : G::BaseField,
    pub u_over_2: G::BaseField,
    pub projection_point: G::CurvePoint,
    pub conic_c: G::BaseField,
    // do we need these?
    pub a: G::BaseField,
    pub b: G::BaseField
}

// C(z, y) : z^2 + 3u/4 y^2 = -f(u)
struct conic<G : AffineCurve> {
    pub z : G::BaseField,
    pub y : G::BaseField
}

// S(u, v, y) : y^2(u^2 + uv + v^2 + a) = -f(u)
struct surface<G : AffineCurve> {
    pub u : G::BaseField,
    pub v : G::BaseField,
    pub w : G::BaseField
}

// do we need this
// V(x1, x2, x3, x4) : f(x1) f(x2) f(x3) = x4^2
struct variety<G : AffineCurve> {
    pub x1 : G::BaseField,
    pub x2 : G::BaseField,
    pub x3 : G::BaseField,
    pub x4 : G::BaseField
}

/* A deterministic function for constructing a valid choice of parameters for a
 given field.
 We start by finding the first `u` satisfying the constraints described above,
 then find the first `y` satisyfing the condition described above. The other
 values are derived from these two choices*.
 *Actually we have one bit of freedom in choosing `z` as z = sqrt(conic_c y^2 - conic_d),
 since there are two square roots.

    let first_map f =
      let rec go i = match f i with Some x -> x | None -> go (i + one) in
      go zero
    in
    let first f = first_map (fun x -> Option.some_if (f x) x) in
    let three_fourths = of_int 3 / of_int 4 in
    let curve_eqn u = (u * u * u) + (a * u) + b in
    let u =
      first (fun u ->
          (* from (15), A = 0, B = Params.a *)
          let check = (three_fourths * u * u) + a in
          let fu = curve_eqn u in
          (not (equal check zero))
          && (not (equal fu zero))
          && not (is_square (negate fu))
      )
    in
    /* The coefficients defining the conic z^2 + c y^2 = d in (15). */
    let conic_c = (three_fourths * u * u) + a in
    let conic_d = negate (curve_eqn u) in
    let projection_point =
      first_map (fun y ->
          let z2 = conic_d - (conic_c * y * y) in
          if F.is_square z2 then Some {Conic.z= F.sqrt z2; y} else None )
    in
    {u; u_over_2= u / of_int 2; conic_c; projection_point; a; b}
*/

fn create_params<G : AffineCurve>() -> SRSParams {
}

/* For a curve z^2 + c y^2 = d and a point (z0, y0) on the curve, there
 is one other point on the curve which is also on the line through (z0, y0)
 with slope t. This function returns that point. */

fn field_to_conic<G : AffineCurve> (t : G::BaseField) -> Conic {
    let z0, y0 = ( params.projection_point.z, params.projection_point.y );
    let ct = params.conic_c * t;
    let s = 2 * ((ct * y0) + z0) / ((ct * t) + G::BaseField::One);
    Conic {
        z : z0 - s,
        y : y0 - (s * t)
    }
}

  /* From (16) : φ(λ) : F → S : λ → ( u, α(λ)/β(λ) - u/2, β(λ) ) */
fn conic_to_s (c : Conic) -> Surface {
    Surface {
    u : params.u,
    v : (z / y) - params.u_over_2,
    w : y
    }
}

fn curve_eqn<G : AffineCurve> ( x : G::BaseField ) -> G::AffinePoint {
    G::AffinePoint {
        x : x,
        y : (x * x * x) + (params.a * x) + params.b
    }
}
  /* This is here for explanatory purposes. See s_to_v_truncated. */
fn _s_to_v (s : Surface) -> Variety {
    let h = (u * u) + (u * v) + (v * v) + params.a;
    Variety {
        w : v,
        x : negate (u + v),
        y : u + (y * y),
        z : curve_eqn(u + (y * y)) * h / y)
    }
}

/* We don't actually need to compute the final coordinate in V */
fn s_to_v_truncated (S : Surface) -> Variety {
    Variety {
        w : v,
        x : negate (u + v),
        y : u + (y * y))
    }
}

fn potential_xs (x : Field) -> Surface {
    s_to_v_truncated(conic_to_s(field_to_conic(x)));
}

fn try_decode ( x : Field ) -> Group {
    let yy = curve_eqn(x);
    if is_square(yy) {
        Some(&Group{
            x : x,
            y : F.sqrt yy)
    } else {
        None
    }
}

// make less ugly :O
fn to_group (x : Field) -> Option(Group) {
  let x1, x2, x3 = potential_xs(t) in
  if let Some(b1) = try_decode(x1) {
      Some(b1)
  } else if let Some(b2) = try_decode(x2) {
      Some(b2)
  } else if let Some(b3) = try_decode(x3) {
      Some(b3)
  } else {
      None
  }
}

fn create(length : usize, num_y_coords : usize) -> CompressedCRS<F> {
}

fn decompress(c : CompressedCRS<F>) -> Vec<(F, F)> {
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_projection_point_well_formed() {
        assert!(on_conic(params.projection_point));
    }

    #[test]
    fn test_field_to_conic() {
        assert!(on_conic(M.field_to_conic(t)));
    }

    #[test]
    fn test_conic_to_s() {
        let z2 = conic_d - (params.conic_c * y * y)
        is_square z2 then Some {Conic.z= sqrt z2; y} else None )
        assert!(on_s(M.conic_to_s(p)));
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
        let (x, y) = to_group(params, t);
        assert_eq!(curve_eqn(x), y * y);
    }

}
