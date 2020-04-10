use std::env;
use std::fs;
use std::path::Path;

use field_assembly::generate_macro_string;

const NUM_LIMBS: usize = 16;

fn main() {
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("field_assembly.rs");

    fs::write(
        &dest_path,
        generate_macro_string(NUM_LIMBS)
    ).unwrap();

    println!("cargo:rerun-if-changed=build.rs");
}
