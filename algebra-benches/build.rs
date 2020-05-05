extern crate rustc_version;

use rustc_version::{version_meta, Channel};

fn main() {
    if version_meta().expect("nightly check failed").channel == Channel::Nightly {
        println!("cargo:rustc-cfg=nightly");
    }
}
