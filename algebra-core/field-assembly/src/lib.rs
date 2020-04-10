extern crate std;

const MAX_REGS: usize = 8;

pub fn generate_macro_string (num_limbs:usize) -> std::string::String {
    let mut macro_string = String::from(
    "macro_rules! asm_mul {
        ($limbs:expr, $a:expr, $b:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {");
    macro_string = generate_matches(num_limbs, macro_string, true);

    macro_string = format!("{}{}", macro_string,
    "macro_rules! asm_square {
        ($limbs:expr, $a:expr, $modulus:expr, $inverse:expr) => {
            match $limbs {");
    macro_string = generate_matches(num_limbs, macro_string, false);
    macro_string
}

fn generate_matches (num_limbs: usize, mut macro_string: String, is_mul: bool) -> String {
    for i in 2..(num_limbs+1) {
        let mut rs_clobber = String::from("");
        let (mut b_declare, mut regs_declare, mut b, mut regs) = ("                   // $3", String::from(""), "$0", "");

        // logic to format macro based on how many limbs there are, whether it is a mul
        for k in 0..i { rs_clobber = format!("{}{}", rs_clobber, format!("\"r{}\", ", 8+k)); }
        let mut limb_specialisation = format!("
                {} => {{", i);
        if is_mul {
            b_declare = ",                  // $3
                              \"r\"(&$b)";
            b = "$4";
            regs_declare = String::from("                        // $4");
        }
        if i > MAX_REGS {
            let extra_reg = if i <= 2*MAX_REGS { 2*(i-MAX_REGS) } else { i };
            limb_specialisation = format!("{}{}", limb_specialisation, format!("
                    let mut regs = [0u64; {}];", extra_reg));
            if is_mul { regs = "$5"; } else { regs = "$4";}
            regs_declare = format!(",                       // ${}
                              \"r\"(&mut regs)                  // {}", 3+(is_mul as usize), regs);
        }

        // Actual asm declaration
        limb_specialisation = format!("{}{}", limb_specialisation, format!("
                    unsafe {{
                        asm!({asm_string}
                            :
                            : \"r\"(&mut $a),                   // $0
                              \"r\"(&$modulus),                 // $1
                              \"i\"(0u64),                      // $2
                              \"i\"($inverse){b_declare}{regs_declare}
                            : \"rcx\", \"rbx\", \"rdx\", \"rax\", {rs_clobber}\"cc\", \"memory\"
                        );
                    }}
                }}",
                asm_string = generate_asm_mul_string(i, "$0", b, regs),
                rs_clobber=rs_clobber,
                b_declare=b_declare,
                regs_declare=regs_declare));
        macro_string = format!("{}{}", macro_string, limb_specialisation);
    }
    macro_string = format!("{}{}", macro_string, format!("
            x => panic!(\"asm_mul (no-carry): number of limbs supported is 2 up to {}. You had {{}}.\", x)
        }};
    }}
}}

", num_limbs));
    macro_string
}

fn generate_asm_mul_string (limbs: usize, a: &str, b: &str, regs: &str) -> String {
    let extra_reg = if limbs <= MAX_REGS { 0 } else { limbs - MAX_REGS };
    let reg_max = std::cmp::min(limbs, MAX_REGS);
    let block_size = if limbs <= MAX_REGS { 0 } else if limbs <= 2*MAX_REGS { limbs-MAX_REGS } else { MAX_REGS };
    let n_spill_blocks = 1 + limbs / MAX_REGS;

    let mut asm_string = String::from("");
    let mut store = "";

    for i in 0..limbs {
        // First inner loop
        if i == 0 {
            asm_string = format!("{}{}", asm_string,format!("\"
                            movq 0({a}), %rdx
                            xorq %rcx, %rcx
                                mulxq 0({b}), %r8, %r9",
                                a=a, b=b));
            for j in 1..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %r{}
                                adcxq %rax, %r{}",
                                j*8, 8 + ((j+1) % limbs), 8+j, b=b));
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rcx
                                mov $2, %rbx
                                adcxq %rax, %r{}
                                adcxq %rbx, %rcx               // %rcx is carry1",
                                (limbs-1)*8, 8+limbs-1, b=b));
        } else {
            asm_string = format!("{}{}", asm_string, format!("
                            movq {}($0), %rdx", i * 8));
            for j in 0..limbs-1 {
                let index_lo = (j+i) % reg_max;
                let index_hi = (j+i+1) % reg_max;
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %r{}",
                                j * 8, 8 + index_lo, 8 + index_hi, b=b));
                // Store the lower word if register spills
                if index_lo < extra_reg {
                    let reg_index = (j+i) % limbs;
                    let block_index = reg_index / MAX_REGS;
                    asm_string = format!("{}{}", asm_string, format!("
                                mov %r{r}, {}({regs})
                                mov {}({regs}), %r{r}",
                                8 * (block_index*block_size + index_lo),
                                8 * (((block_index+block_size) % n_spill_blocks) + index_lo),
                                r = 8 + index_lo, regs=regs));
                }
            }
            let index_lo = (i+limbs-1) % reg_max;
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rcx
                                mov $2, %rbx
                                adcxq %rax, %r{}
                                adoxq %rbx, %rcx
                                adcxq %rbx, %rcx",
                                (limbs-1) * 8,
                                8 + index_lo,
                                b=b));
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
                                mov $2, %rbx
                                adcxq %rax, %r{}
                                adoxq %rcx, %r{2}
                                adcxq %rbx, %r{2}",
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
