mod fixed_base;
mod variable_base;
pub use fixed_base::*;
pub use variable_base::*;

/// The result of this function is only approximately `ln(a)`
/// [`Explanation of usage`]
///
/// [`Explanation of usage`]: https://github.com/scipr-lab/zexe/issues/79#issue-556220473
fn ln_without_floats(a: usize) -> usize {
    fn log2(x: usize) -> u32 {
        if x <= 1 {
            return 0;
        }

        let n = x.leading_zeros();
        ::core::mem::size_of::<usize>() as u32 * 8 - n
    }

    // log2(a) * ln(2)
    (log2(a) * 69 / 100) as usize
}
