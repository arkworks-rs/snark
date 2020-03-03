use algebra::curves::AffineCurve;
use core::fmt::Error;

#[cfg(test)]
mod tests {
    #[test]
    fn test_full_map() {
        let (x, y) = to_group(parameters, t);
        assert_eq!(curve_eqn(x), y * y);
    }
}
