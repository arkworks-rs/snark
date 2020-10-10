#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public)]
#![cfg_attr(use_asm, feature(llvm_asm))]
#![cfg_attr(
    not(any(use_asm, feature = "prefetch", use_bw6_asm)),
    forbid(unsafe_code)
)]

#[cfg(all(test, not(feature = "std")))]
#[macro_use]
extern crate std;

/// This crate needs to be public, because we expose the `to_bytes!` macro.
/// See the similar issue [`smallvec#198`]
///
/// [`smallvec#198`]: https://github.com/servo/rust-smallvec/pull/198
#[cfg(not(feature = "std"))]
#[macro_use]
#[doc(hidden)]
pub extern crate alloc;

#[cfg(any(feature = "timing", feature = "timing_detailed"))]
pub extern crate backtrace;

#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
#[doc(hidden)]
pub use alloc::{
    borrow::{Cow, ToOwned},
    boxed::Box,
    collections::btree_map::BTreeMap,
    format,
    string::String,
    vec,
    vec::Vec,
};

#[cfg(feature = "std")]
#[allow(unused_imports)]
#[doc(hidden)]
pub use std::{
    borrow::{Cow, ToOwned},
    boxed::Box,
    collections::btree_map::BTreeMap,
    format,
    string::String,
    vec,
    vec::Vec,
};

#[macro_use]
pub mod timing;
pub use timing::*;

#[macro_use]
extern crate derivative;

#[cfg_attr(test, macro_use)]
pub mod bytes;
pub use self::bytes::*;

#[macro_use]
pub mod serialize;
pub use self::serialize::*;

#[macro_use]
pub mod fields;
pub use self::fields::*;

pub mod biginteger;
pub use self::biginteger::*;

pub mod curves;
pub use self::curves::*;

pub mod groups;
pub use self::groups::*;

mod rand;
pub use self::rand::*;

mod error;
pub use self::error::*;

mod to_field_vec;
pub use to_field_vec::ToConstraintField;

pub mod msm;
pub use self::msm::*;

pub use num_traits::{One, Zero};

pub mod prelude {
    pub use crate::biginteger::BigInteger;

    pub use crate::fields::{Field, FpParameters, PrimeField, SquareRootField};

    pub use crate::groups::Group;

    pub use crate::curves::{AffineCurve, PairingEngine, ProjectiveCurve};

    pub use crate::rand::UniformRand;

    pub use num_traits::{One, Zero};

    pub use crate::error::*;
}

#[cfg(not(feature = "std"))]
pub mod io;

#[cfg(feature = "std")]
pub use std::io;

#[cfg(feature = "derive")]
#[allow(unused_imports)]
#[macro_use]
extern crate algebra_core_derive;

#[cfg(not(feature = "std"))]
fn error(_msg: &'static str) -> io::Error {
    io::Error
}

#[cfg(feature = "std")]
fn error(msg: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::Other, msg)
}

/// Returns the base-2 logarithm of `x`.
///
/// ```
/// use algebra_core::log2;
///
/// assert_eq!(log2(16), 4);
/// assert_eq!(log2(17), 5);
/// assert_eq!(log2(1), 0);
/// assert_eq!(log2(0), 0);
/// assert_eq!(log2(usize::MAX), (core::mem::size_of::<usize>() * 8) as u32);
/// assert_eq!(log2(1 << 15), 15);
/// assert_eq!(log2(2usize.pow(18)), 18);
/// ```
pub fn log2(x: usize) -> u32 {
    if x == 0 {
        0
    } else if x.is_power_of_two() {
        1usize.leading_zeros() - x.leading_zeros()
    } else {
        0usize.leading_zeros() - x.leading_zeros()
    }
}

/// Prefetches as many cache lines as is occupied by the type T.
/// We assume 64B cache lines
#[cfg(feature = "prefetch")]
#[inline(always)]
pub fn prefetch<T>(p: *const T) {
    unsafe {
        match n_lines::<T>() {
            1 => unroll!(1, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            2 => unroll!(2, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            3 => unroll!(3, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            4 => unroll!(4, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            5 => unroll!(5, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            6 => unroll!(6, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            7 => unroll!(7, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            8 => unroll!(8, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            9 => unroll!(9, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            10 => unroll!(10, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            11 => unroll!(11, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            12 => unroll!(12, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            13 => unroll!(13, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            14 => unroll!(14, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            15 => unroll!(15, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
            _ => unroll!(16, |i| core::arch::x86_64::_mm_prefetch(
                (p as *const i8).offset(i * 64),
                core::arch::x86_64::_MM_HINT_T0
            )),
        }
    }
}

#[cfg(feature = "prefetch")]
#[inline]
pub fn clear_cache<T>(p: *const T) {
    unsafe {
        match n_lines::<T>() {
            1 => unroll!(1, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            2 => unroll!(2, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            3 => unroll!(3, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            4 => unroll!(4, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            5 => unroll!(5, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            6 => unroll!(6, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            7 => unroll!(7, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            8 => unroll!(8, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            9 => unroll!(9, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            10 => unroll!(10, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            11 => unroll!(11, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            12 => unroll!(12, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            13 => unroll!(13, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            14 => unroll!(14, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            15 => unroll!(15, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
            _ => unroll!(16, |i| core::arch::x86_64::_mm_clflush(
                (p as *const u8).offset(i * 64)
            )),
        }
    }
}

#[cfg(feature = "prefetch")]
const fn n_lines<T>() -> isize {
    ((std::mem::size_of::<T>() - 1) / 64 + 1) as isize
}

#[macro_export]
macro_rules! unroll {
    (0, |$i:ident| $s:stmt) => {};
    (1, |$i:ident| $s:stmt) => {{
        let $i: isize = 0;
        $s
    }};
    (2, |$i:ident| $s:stmt) => {{
        unroll!(1, |$i| $s);
        let $i: isize = 1;
        $s
    }};
    (3, |$i:ident| $s:stmt) => {{
        unroll!(2, |$i| $s);
        let $i: isize = 2;
        $s
    }};
    (4, |$i:ident| $s:stmt) => {{
        unroll!(3, |$i| $s);
        let $i: isize = 3;
        $s
    }};
    (5, |$i:ident| $s:stmt) => {{
        unroll!(4, |$i| $s);
        let $i: isize = 4;
        $s
    }};
    (6, |$i:ident| $s:stmt) => {{
        unroll!(5, |$i| $s);
        let $i: isize = 5;
        $s
    }};
    (7, |$i:ident| $s:stmt) => {{
        unroll!(6, |$i| $s);
        let $i: isize = 6;
        $s
    }};
    (8, |$i:ident| $s:stmt) => {{
        unroll!(7, |$i| $s);
        let $i: isize = 7;
        $s
    }};
    (9, |$i:ident| $s:stmt) => {{
        unroll!(8, |$i| $s);
        let $i: isize = 8;
        $s
    }};
    (10, |$i:ident| $s:stmt) => {{
        unroll!(9, |$i| $s);
        let $i: isize = 9;
        $s
    }};
    (11, |$i:ident| $s:stmt) => {{
        unroll!(10, |$i| $s);
        let $i: isize = 10;
        $s
    }};
    (12, |$i:ident| $s:stmt) => {{
        unroll!(11, |$i| $s);
        let $i: isize = 11;
        $s
    }};
    (13, |$i:ident| $s:stmt) => {{
        unroll!(12, |$i| $s);
        let $i: isize = 12;
        $s
    }};
    (14, |$i:ident| $s:stmt) => {{
        unroll!(13, |$i| $s);
        let $i: isize = 13;
        $s
    }};
    (15, |$i:ident| $s:stmt) => {{
        unroll!(14, |$i| $s);
        let $i: isize = 14;
        $s
    }};
    (16, |$i:ident| $s:stmt) => {{
        unroll!(15, |$i| $s);
        let $i: isize = 15;
        $s
    }};
}

#[macro_export]
macro_rules! cfg_iter {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_iter();

        #[cfg(not(feature = "parallel"))]
        let result = $e.iter();

        result
    }};
}

#[macro_export]
macro_rules! cfg_iter_mut {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_iter_mut();

        #[cfg(not(feature = "parallel"))]
        let result = $e.iter_mut();

        result
    }};
}

#[macro_export]
macro_rules! cfg_chunks_mut {
    ($e: expr, $N: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_chunks_mut($N);

        #[cfg(not(feature = "parallel"))]
        let result = $e.chunks_mut($N);

        result
    }};
}
