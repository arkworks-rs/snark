extern crate std;

#[macro_use]
pub mod assembler;
pub mod arithmetic;
pub mod context;

use assembler::*;
use arithmetic as ar;
use context::*;

const MAX_REGS: usize = 6;

pub fn generate_macro_string (num_limbs:usize) -> std::string::String {
    if num_limbs > 3 * MAX_REGS {
        panic!("Number of limbs must be <= {} and MAX_REGS >= 6", 3*MAX_REGS);
    }
    let mut macro_string = String::from(
    "macro_rules! asm_mul {
        ($limbs:expr, $a:expr, $b:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {");
    macro_string = format!("{}{}", macro_string, generate_matches(num_limbs, true));

    macro_string = format!("{}{}", macro_string,
    "macro_rules! asm_square {
        ($limbs:expr, $a:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {");
    macro_string = format!("{}{}", macro_string, generate_matches(num_limbs, false));
    macro_string
}

fn generate_asm_mul_string (ctx: &Context, limbs: usize) -> String {
    let a = ctx.clone().get("a");
    let b = ctx.clone().try_get("b", "a");
    let modulus = ctx.clone().get("modulus");
    let zero = ctx.clone().get("0");
    let inverse = ctx.clone().get("inverse");

    generate_array!(a0, a1, a, limbs);
    generate_array!(b0, b1, b, limbs);
    generate_array!(m, m1, modulus, limbs);
    // if limbs > 8 {
    //     generate_array!(s, s1, spills, limbs * 2);
    // }

    let mut asm = Assembler::new(limbs);

    asm.begin();

    // if limbs <= 8 {
        asm.xorq(RCX, RCX);
        for i in 0..limbs {
            if i == 0 {
                ar::mul_1(&mut asm, a1[0], &b1, &zero);
            } else {
                ar::mul_add_1(&mut asm, &a1, &b1, &zero, i);
            }
            ar::mul_add_shift_1(&mut asm, &m1, &inverse, &zero, i);
        }
        for i in 0..asm.limbs {
            asm.movq(R[i], a1[i]);
        }

    // } else {
    //     asm.xorq(RCX, RCX);
    //     for i in 0..8 {
    //         if i == 0 {
    //             ar::mul_1_mov(&mut asm, a1[0], &b1, 0);
    //         } else {
    //             ar::mul_add_1(&mut asm, &a1, &b1, i);
    //         }
    //     }
    //     for i in 0..8 {
    //         ar::mul_add_1(&mut asm, &m1, 0);
    //     }
    //     for i in 0..asm.limbs {
    //         asm.movq(R[i], a1[i]);
    //     }
    //
    // }

    asm.end();

    asm.get_asm_string()
}

fn generate_matches (num_limbs: usize, is_mul: bool) -> String {
    let mut ctx = Context::new();
    for limbs in 2..(num_limbs+1) {
        ctx.reset();

        ctx.add_declaration("a", "r", "&mut $a");
        if is_mul { ctx.add_declaration("b", "r", "&$b"); }
        ctx.add_declaration("modulus", "r", "&$modulus");
        ctx.add_declaration("0", "i", "0u64");
        ctx.add_declaration("inverse", "i", "$inverse");

        ctx.add_limb(limbs);
        if limbs > 8 {
            ctx.add_buffer(2*limbs);
            ctx.add_declaration("buf", "r", "&mut spill_buffer");
        }

        let asm_string = generate_asm_mul_string(&ctx, limbs);

        ctx.add_asm(asm_string);
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
