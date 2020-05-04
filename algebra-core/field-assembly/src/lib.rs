extern crate std;

#[macro_use]
pub mod utils;
use utils::*;

pub mod context;
use context::*;

use mince::assemble;

use std::cell::RefCell;

const MAX_REGS: usize = 6;

pub fn generate_macro_string(num_limbs: usize) -> std::string::String {
    if num_limbs > 3 * MAX_REGS {
        panic!(
            "Number of limbs must be <= {} and MAX_REGS >= 6",
            3 * MAX_REGS
        );
    }
    let mut macro_string = String::from(
        "
        macro_rules! llvm_asm_mul {
        ($limbs:expr, $a:expr, $b:expr, $modulus:expr, $mod_prime:expr) => {
            match $limbs {",
    );
    macro_string += &generate_matches(num_limbs, true);

    macro_string += &"
        macro_rules! llvm_asm_square {
        ($limbs:expr, $a:expr, $modulus:expr, $mod_prime:expr) => {
            match $limbs {";
    macro_string += &generate_matches(num_limbs, false);
    macro_string
}

#[assemble]
fn generate_llvm_asm_mul_string(
    a: &str,
    b: &str,
    modulus: &str,
    zero: &str,
    mod_prime: &str,
    limbs: usize,
) -> String {
    reg!(a0, a1, a, limbs);
    reg!(b0, b1, b, limbs);
    reg!(m, m1, modulus, limbs);

    xorq(RCX, RCX);
    for i in 0..limbs {
        if i == 0 {
            mul_1!(a1[0], b1, zero, limbs);
        } else {
            mul_add_1!(a1, b1, zero, i, limbs);
        }
        mul_add_shift_1!(m1, mod_prime, zero, i, limbs);
    }
    for i in 0..limbs {
        movq(R[i], a1[i]);
    }
}

fn generate_matches(num_limbs: usize, is_mul: bool) -> String {
    let mut ctx = Context::new();
    for limbs in 2..(num_limbs + 1) {
        ctx.reset();

        ctx.add_declaration("a", "r", "&mut $a");
        if is_mul {
            ctx.add_declaration("b", "r", "&$b");
        }
        ctx.add_declaration("modulus", "r", "&$modulus");
        ctx.add_declaration("0", "i", "0u64");
        ctx.add_declaration("mod_prime", "i", "$mod_prime");

        ctx.add_limb(limbs);
        if limbs > MAX_REGS {
            ctx.add_buffer(2 * limbs);
            ctx.add_declaration("buf", "r", "&mut spill_buffer");
        }

        let llvm_asm_string = generate_llvm_asm_mul_string(
            &ctx.clone().get("a"),
            &ctx.clone().try_get("b", "a"),
            &ctx.clone().get("modulus"),
            &ctx.clone().get("0"),
            &ctx.clone().get("mod_prime"),
            limbs,
        );

        ctx.add_llvm_asm(llvm_asm_string);
        ctx.add_clobber_from_vec(vec!["rcx", "rbx", "rdx", "rax"]);
        for j in 0..std::cmp::min(limbs, 8) {
            ctx.add_clobber(REG_CLOBBER[j]);
        }
        ctx.add_clobber_from_vec(vec!["cc", "memory"]);
        ctx.build();
    }
    ctx.end(num_limbs);
    ctx.get_string()
}
