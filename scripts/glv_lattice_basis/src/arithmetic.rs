use algebra_core::biginteger::BigInteger;

// Naive long division
pub fn div_with_remainder<BigInt: BigInteger>(
    numerator: BigInt,
    divisor: BigInt,
) -> (BigInt, BigInt) {
    assert!(divisor != BigInt::from(0), "Divisor cannot be zero");
    let mut remainder = numerator;
    let mut quotient = BigInt::from(0);

    let div_num_bits = divisor.num_bits();

    while remainder >= divisor {
        let mut current_divisor = divisor;
        let mut num_bits = 1 + remainder.num_bits() - div_num_bits;
        current_divisor.muln(num_bits);
        while current_divisor > remainder {
            current_divisor.div2();
            num_bits -= 1;
        }
        remainder.sub_noborrow(&current_divisor);

        let mut pow2_quot = BigInt::from(1);
        pow2_quot.muln(num_bits);
        quotient.add_nocarry(&pow2_quot);
    }

    let mut reconstructed_numerator =
        BigInt::mul_no_reduce_lo(&quotient.as_ref(), &divisor.as_ref());
    reconstructed_numerator.add_nocarry(&remainder);
    assert_eq!(reconstructed_numerator, numerator);
    (quotient, remainder)
}
