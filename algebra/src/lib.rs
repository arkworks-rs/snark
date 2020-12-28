#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public)]
#![cfg_attr(use_asm, feature(llvm_asm))]
#![cfg_attr(not(use_asm), forbid(unsafe_code))]
#![cfg_attr(use_asm, deny(unsafe_code))]

#[macro_use]
extern crate derivative;

#[cfg_attr(test, macro_use)]
pub mod bytes;
pub use self::bytes::*;

pub mod bits;
pub use self::bits::*;

pub mod biginteger;
pub use self::biginteger::*;

pub mod curves;
pub use self::curves::*;

#[macro_use]
pub mod fields;
pub use self::fields::*;

pub mod groups;
pub use self::groups::*;

pub mod validity;
pub use self::validity::*;

mod rand;
pub use self::rand::*;

mod to_field_vec;
pub use to_field_vec::ToConstraintField;

#[cfg(feature = "parallel")]
pub mod msm;
#[cfg(feature = "parallel")]
pub use self::msm::*;

#[cfg(feature = "fft")]
pub mod fft;
#[cfg(feature = "fft")]
pub use self::fft::*;

pub type Error = Box<dyn std::error::Error>;

/// Returns the base-2 logarithm of `x`.
/// ```
/// use algebra::log2;
///
/// assert_eq!(log2(16), 4);
/// assert_eq!(log2(17), 5);
/// assert_eq!(log2(1), 0);
/// assert_eq!(log2(0), 0);
/// assert_eq!(log2(usize::MAX), (std::mem::size_of::<usize>() * 8) as u32);
/// assert_eq!(log2(1 << 15), 15);
/// assert_eq!(log2(2usize.pow(18)), 18);
/// ```
#[inline]
pub fn log2(x: usize) -> u32 {
    if x == 0 {
        0
    } else if x.is_power_of_two() {
        1usize.leading_zeros() - x.leading_zeros()
    } else {
        0usize.leading_zeros() - x.leading_zeros()
    }
}

/// Returns the floor of the base-2 logarithm of `x`.
/// ```
/// use algebra::log2_floor;
///    assert_eq!(log2_floor(0), 0);
///    assert_eq!(log2_floor(1), 0);
///    assert_eq!(log2_floor(2), 1);
///    assert_eq!(log2_floor(3), 1);
///    assert_eq!(log2_floor(4), 2);
///    assert_eq!(log2_floor(5), 2);
///    assert_eq!(log2_floor(6), 2);
///    assert_eq!(log2_floor(7), 2);
///    assert_eq!(log2_floor(8), 3);
/// ```
#[inline]
pub fn log2_floor(x: usize) -> u32 {
    if x == 0 {
        0
    } else {
        (x as f64).log2() as u32
    }
}