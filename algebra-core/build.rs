use std::env;
use std::fs;
use std::path::Path;

extern crate rustc_version;
use rustc_version::{version_meta, Channel};

use field_assembly::generate_macro_string;

const NUM_LIMBS: usize = 8;

fn main() {
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("field_assembly.rs");
    let is_nightly = version_meta().channel == Channel::Nightly;

    if cfg!(feature = "llvm_asm") && is_nightly {
        fs::write(&dest_path, generate_macro_string(NUM_LIMBS)).unwrap();
    } else {
        fs::write(&dest_path, "").unwrap();
    }

    println!("cargo:rerun-if-changed=build.rs");

    if is_nightly {
        println!("cargo:rustc-cfg=nightly");
    }
}
