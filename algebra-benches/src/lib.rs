#![cfg_attr(all(feature = "nightly", not(feature = "stable")), feature(test))]
#![allow(unused_macros)]

#[cfg(all(feature = "nightly", not(feature = "stable")))]
extern crate test;

#[macro_use]
pub mod macros;

#[cfg(all(feature = "nightly", not(feature = "stable")))]
mod curves;
