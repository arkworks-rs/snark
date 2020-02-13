use algebra::curves::AffineCurve;
use core::fmt::Error;

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
