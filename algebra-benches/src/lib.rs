#![cfg_attr(feature = "nightly", feature(test))]
#![allow(unused_macros)]

#[cfg(feature = "nightly")]
extern crate test;

#[macro_use]
pub mod macros;

#[cfg(feature = "nightly")]
mod curves;
