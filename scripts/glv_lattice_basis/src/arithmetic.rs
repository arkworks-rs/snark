use algebra_core::bigint::BigInteger;

// Naive long division
fn div_with_remainder<BigInt: BigInteger>(
    numerator: BigInt,
    divisor: BigInt
) -> (BigInt, BigInt)
{
    assert!(divisor != BigInt::from(0));
    let mut remainder = numerator;
    let mut quotient = BigInt::from(0);
    let limbs = BigIntNum::NUM_LIMBS;
    while remainder >= divisor {
        let mut current_divisor = divisor;
        let mut i = 0;
        while remainder.0[limbs - i - 1] == 0u64 && i + 1 < limbs {
            i += 1;
        }
        let biggest_non_zero = limbs - i - 1;
        let num_bits_non_zero = (biggest_non_zero * 64)
            - remainder.0[biggest_non_zero].leading_zeros();

        current_divisor.muln(num_bits_non_zero);

        let mut n_bits = num_bits_non_zero;
        while current_divisor > remainder {
            current_divisor.div2();
            n_bits -= 1;
        }
        remainder -= current_divisor;

        let mut pow2_quot = BigInt::from(1);
        pow2_quot.muln(n_bits);
        quotient += pow2_quot;
    }
    assert_eq!(quotient.mul_no_reduce(&divisor) + remainder, numerator);
    (quotient, remainder)
}
