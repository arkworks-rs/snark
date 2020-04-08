use std::env;
use std::fs;
use std::path::Path;

const MAX_LIMBS: usize = 8;

fn main() {
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("assembly.rs");

    let macro_string = generate_macro_string(MAX_LIMBS);

    fs::write(
        &dest_path,
        macro_string
    ).unwrap();

    println!("cargo:rerun-if-changed=build.rs");
}

// // Different strategies for limbs <= 6, 7 <= limbs <= 12, limbs > 12?
// fn generate_mul_add_1_asm (limbs: usize) -> String {
//     let mut asm_string = String::from("");
//     for i in 0..limbs {


fn generate_macro_string (max_limbs:usize) -> std::string::String {
    let mut macro_string = String::from(
    "macro_rules! asm_mul {
        ($a:expr, $b:expr, $limbs:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {
    ");
    for i in 2..(max_limbs+1) {
        let mut rs = String::from("");
        for k in 0..i {
            rs = format!("{}{}", rs, format!("\"r{}\", ", 8+k));
        }
        let limb_specialisation = format!(
    "           {} => {{
                    unsafe {{
                        asm!({}
                            :
                            : \"r\"(&mut $a),                            // $0
                              \"r\"(&$b),                                // $1
                              \"r\"(&$modulus),                          // $2
                              \"i\"(0u64),                               // $3
                              \"i\"($inverse)                            // $4
                            : \"rcx\", \"rbx\", \"rdx\", \"rax\", {} \"cc\", \"memory\"
                        );
                    }}
                }}

    ", i, generate_asm_mul_string(i), rs);//ASM_STR, rs);//
        macro_string = format!("{}{}", macro_string, limb_specialisation);
    }
    macro_string = format!("{}{}", macro_string,
                        "x => panic!(\"asm_mul (no-carry): number of limbs supported is 2 up to 8. You had {}\", x)
        };
    }
}");
    macro_string = format!("{}{}", macro_string,
    "macro_rules! asm_square {
        ($a:expr, $limbs:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {
    ");
    for i in 2..(max_limbs+1) {
        let mut rs = String::from("");
        for k in 0..i {
            rs = format!("{}{}", rs, format!("\"r{}\", ", 8+k));
        }
        let limb_specialisation = format!("
            {} => {{
                    unsafe {{
                        asm!({}
                            :
                            : \"r\"(&mut $a),                            // $0
                              \"r\"(&$modulus),                          // $1
                              \"i\"(0u64),                               // $2
                              \"i\"($inverse)                            // $3
                            : \"rcx\", \"rbx\", \"rdx\", \"rax\", {} \"cc\", \"memory\"
                        );
                    }}
                }}

    ", i, generate_asm_square_string(i), rs);//ASM_STR, rs);//
        macro_string = format!("{}{}", macro_string, limb_specialisation);
    }
    macro_string = format!("{}{}", macro_string,
                        "x => panic!(\"asm_mul (no-carry): number of limbs supported is 2 up to 8. You had {}\", x)
        };
    }
    }");
    macro_string
}

fn generate_asm_square_string (limbs: usize) -> String {
    let mut asm_string = String::from("");
    for i in 0..limbs {
        // First inner loop
        if i == 0 {
            asm_string = format!("{}{}", asm_string,"\"
                            movq 0($0), %rdx
                            xorq %rcx, %rcx
                                mulxq 0($0), %r8, %r9");
            for j in 1..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($0), %rax, %r{}
                                adcxq %rax, %r{}",
                                j*8, 8 + ((j+1) % limbs), 8+j));
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($0), %rax, %rcx
                                mov $2, %rdx
                                adcxq %rax, %r{}
                                adcxq %rdx, %rcx               // %rcx is carry1",
                                (limbs-1)*8, 8+limbs-1));
        } else {
            asm_string = format!("{}{}", asm_string, format!("
                            movq {}($0), %rdx", i * 8));
            for j in 0..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($0), %rax, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %r{}",
                                j * 8,
                                8 + ((j+i) % limbs),
                                8 + ((j+i+1) % limbs)));
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($0), %rax, %rcx
                                mov $2, %rdx
                                adcxq %rax, %r{}
                                adoxq %rdx, %rcx
                                adcxq %rdx, %rcx",
                                (limbs-1) * 8,
                                8 + ((i+limbs-1) % limbs)));
        }
        // Second inner loop
        asm_string = format!("{}{}", asm_string, format!("
                            movq $3, %rdx
                            mulxq %r{}, %rdx, %rax            // wrapping_mul", 8+i));
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq 0($1), %rax, %rbx
                                adcxq %r{}, %rax              // put junk in rax
                                adoxq %rbx, %r{}",
                                8 + (i % limbs),
                                8 + ((i+1) % limbs)));
        for j in 1..limbs-1 {
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %r{}",
                                j * 8,
                                8 + ((j+i) % limbs),
                                8 + ((j+i+1) % limbs)));
        }
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %r{2}
                                mov $2, %rdx
                                adcxq %rax, %r{}
                                adoxq %rcx, %r{2}
                                adcxq %rdx, %r{2}",
                                (limbs-1)*8,
                                8 + ((i+limbs-1) % limbs),
                                8 + ((i) % limbs)));
    }
    for i in 0..limbs {
        asm_string = format!("{}{}", asm_string, format!("
                            movq %r{}, {}($0)", 8+(i % limbs), i*8));
    }
    format!("{}{}",asm_string, "\"")
}

// For now, generated code only works for up to  8/10 limbss
// In the future, we can try to implement data movement to and from an address
// for higher number of limbs
fn generate_asm_mul_string (limbs: usize) -> String {
    let mut asm_string = String::from("");
    for i in 0..limbs {
        // First inner loop
        if i == 0 {
            asm_string = format!("{}{}", asm_string,"\"
                            movq 0($0), %rdx
                            xorq %rcx, %rcx
                                mulxq 0($1), %r8, %r9");
            for j in 1..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %r{}
                                adcxq %rax, %r{}",
                                j*8, 8 + ((j+1) % limbs), 8+j));
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %rcx
                                mov $3, %rdx
                                adcxq %rax, %r{}
                                adcxq %rdx, %rcx               // %rcx is carry1",
                                (limbs-1)*8, 8+limbs-1));
        } else {
            asm_string = format!("{}{}", asm_string, format!("
                            movq {}($0), %rdx", i * 8));
            for j in 0..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %r{}",
                                j * 8,
                                8 + ((j+i) % limbs),
                                8 + ((j+i+1) % limbs)));
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %rcx
                                mov $3, %rdx
                                adcxq %rax, %r{}
                                adoxq %rdx, %rcx
                                adcxq %rdx, %rcx",
                                (limbs-1) * 8,
                                8 + ((i+limbs-1) % limbs)));
        }
        // Second inner loop
        asm_string = format!("{}{}", asm_string, format!("
                            movq $4, %rdx
                            mulxq %r{}, %rdx, %rax            // wrapping_mul", 8+i));
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq 0($2), %rax, %rbx
                                adcxq %r{}, %rax              // put junk in rax
                                adoxq %rbx, %r{}",
                                8 + (i % limbs),
                                8 + ((i+1) % limbs)));
        for j in 1..limbs-1 {
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($2), %rax, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %r{}",
                                j * 8,
                                8 + ((j+i) % limbs),
                                8 + ((j+i+1) % limbs)));
        }
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($2), %rax, %r{2}
                                mov $3, %rdx
                                adcxq %rax, %r{}
                                adoxq %rcx, %r{2}
                                adcxq %rdx, %r{2}",
                                (limbs-1)*8,
                                8 + ((i+limbs-1) % limbs),
                                8 + ((i) % limbs)));
    }
    for i in 0..limbs {
        asm_string = format!("{}{}", asm_string, format!("
                            movq %r{}, {}($0)", 8+(i % limbs), i*8));
    }
    format!("{}{}",asm_string, "\"")
}
