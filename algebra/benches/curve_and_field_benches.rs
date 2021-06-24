#![cfg_attr(nightly, feature(test))]
#![allow(unused_macros)]
#![feature(test)]

extern crate test;

#[macro_use]
pub mod macros;

mod curves;

#[cfg(feature = "fft")]
mod fft;