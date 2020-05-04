#![cfg_attr(nightly, feature(test))]
#![allow(unused_macros)]

#[cfg(nightly)]
extern crate test;

#[macro_use]
pub mod macros;

#[cfg(nightly)]
mod curves;
