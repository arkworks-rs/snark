extern crate rustc_version;

use rustc_version::{version_meta, Channel};

fn main() {
    if version_meta().channel.expect("nightly check failed") == Channel::Nightly {
        println!("cargo:rustc-cfg=nightly");
    }
}
