use algebra_core::{fields::FpParameters, PrimeField};
use core::cmp::min;

pub(crate) fn mixed_radix_fft_permute(
    two_adicity: u32,
    q_adicity: u32,
    q: usize,
    n: usize,
    mut i: usize,
) -> usize {
    // This is the permutation obtained by splitting into 2 groups two_adicity times
    // and then q groups q_adicity many times. It can be efficiently described
    // as follows i = 2^0 b_0 + 2^1 b_1 + ... + 2^{two_adicity - 1}
    // b_{two_adicity - 1} + 2^two_adicity ( x_0 + q^1 x_1 + .. +
    // q^{q_adicity-1} x_{q_adicity-1}) We want to return
    // j = b_0 (n/2) + b_1 (n/ 2^2) + ... + b_{two_adicity-1} (n/ 2^two_adicity)
    // + x_0 (n / 2^two_adicity / q) + .. + x_{q_adicity-1} (n / 2^two_adicity /
    // q^q_adicity)
    let mut res = 0;
    let mut shift = n;

    for _ in 0..two_adicity {
        shift /= 2;
        res += (i % 2) * shift;
        i /= 2;
    }

    for _ in 0..q_adicity {
        shift /= q;
        res += (i % q) * shift;
        i /= q;
    }

    res
}

pub(crate) fn best_mixed_domain_size<F: PrimeField>(min_size: usize) -> usize {
    let mut best = usize::max_value();
    let small_subgroup_power = F::Params::SMALL_SUBGROUP_POWER.unwrap();
    let small_subgroup_base = F::Params::SMALL_SUBGROUP_BASE.unwrap() as usize;

    for b in 0..=small_subgroup_power {
        let mut r = small_subgroup_base.pow(b);

        let mut two_adicity = 0;
        while r < min_size {
            r *= 2;
            two_adicity += 1;
        }

        if two_adicity <= F::Params::TWO_ADICITY {
            best = min(best, r);
        }
    }

    best
}
