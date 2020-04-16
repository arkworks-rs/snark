use std::env;
use std::fs;
use std::path::Path;

#[cfg(feature = "asm")]
use field_assembly::generate_macro_string;

#[cfg(feature = "asm")]
const NUM_LIMBS: usize = 8;


fn main() {
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("field_assembly.rs");


    #[cfg(feature = "asm")]
    fs::write(
        &dest_path,
        generate_macro_string(NUM_LIMBS)
    ).unwrap();


    #[cfg(not(feature = "asm"))]
    fs::write(
        &dest_path,
        ""
    ).unwrap();

    println!("cargo:rerun-if-changed=build.rs");
}
