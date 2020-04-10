extern crate std;

const MAX_REGS: usize = 5;

// Only works for up to
pub fn generate_macro_string (num_limbs:usize) -> std::string::String {
    if (num_limbs > 2 * MAX_REGS) || (MAX_REGS < 4) {
        panic!("Number of limbs must be <= {} and MAX_REGS >= 4", 2*MAX_REGS);
    }
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
        let mut limb_specialisation = format!("
                {} => {{", i);
        // logic to format macro based on how many limbs there are, whether it is a mul
        let (mut b_declare, mut regs_declare, mut b, mut regs) = ("                   // $3", String::from(""), "$0", "");
        let mut rs_clobber = String::from("");
        for k in 0..i { rs_clobber = format!("{}{}", rs_clobber, format!("\"r{}\", ", 8+k)); }
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

fn generate_asm_mul_string (limbs: usize, a: &str, b: &str, spill: &str) -> String {
    let mut asm_string = String::from("");
    let spilled = std::collections::HashMap::new();

    for i in 0..limbs {
        // First inner loop
        if i == 0 {
            asm_string = format!("{}{}", asm_string,format!("\"
                            movq 0({a}), %rdx
                            xorq %rcx, %rcx
                                mulxq 0({b}), %r8, %r9",
                                a=a, b=b));
            if is_spill(limbs, 0) {
                asm_string = spill_swap(asm_string, 0, limbs, spill);
                spilled.insert("%r8", 0);
            }

            for j in 1..limbs-1 {

                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, {}
                                adcxq %rax, {}",
                                j*8, reg_from_index(limbs, (j+1) % limbs),
                                reg_from_index(limbs, j), b=b));
                if is_spill(limbs, j) {
                    asm_string = spill_swap(asm_string, j, limbs, spill);
                    spilled.insert(&reg_from_index(limbs, (j+1) % limbs), j);
                }
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rcx
                                mov $2, %rbx
                                adcxq %rax, {}
                                adcxq %rbx, %rcx               // %rcx is carry1",
                                (limbs-1)*8, reg_from_index(limbs, limbs-1), b=b));
            if is_spill(limbs, limbs-1) {
                asm_string = spill_swap(asm_string, limbs-1, limbs, spill);
                spilled.insert(&reg_from_index(limbs, limbs-1), limbs-1);
            }
        } else {
            asm_string = format!("{}{}", asm_string, format!("
                            movq {}($0), %rdx", i * 8));
            for j in 0..limbs-1 {
                asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rbx
                                adcxq %rax, {}
                                adoxq %rbx, {}",
                                j * 8, reg_from_index(limbs, (j+i) % limbs), reg_from_index(limbs, (j+i+1) % limbs), b=b));
                if is_spill(limbs, (j+i) % limbs) {
                    asm_string = spill_swap(asm_string, (j+i) % limbs, limbs, spill);
                    spilled.insert(&reg_from_index(limbs, (j+i) % limbs), (j+i) % limbs);
                }
            }
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}({b}), %rax, %rcx
                                mov $2, %rbx
                                adcxq %rax, {}
                                adoxq %rbx, %rcx
                                adcxq %rbx, %rcx",
                                (limbs-1) * 8, reg_from_index(limbs, (i+limbs-1) % limbs), b=b));
            if is_spill(limbs, (i+limbs-1) % limbs) { asm_string = spill_swap(asm_string, (i+limbs-1) % limbs, limbs, spill); }
        }
            // Second inner loop
        asm_string = format!("{}{}", asm_string, format!("
                            movq $3, %rdx
                            mulxq %r{}, %rdx, %rax            // wrapping_mul", 8+i));
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq 0($1), %rax, %rbx
                                adcxq {}, %rax              // put junk in rax
                                adoxq %rbx, {}",
                                reg_from_index(limbs, i % limbs),
                                reg_from_index(limbs, (i+1) % limbs)));
        for j in 1..limbs-1 {
            asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, %rbx
                                adcxq %rax, {}
                                adoxq %rbx, {}",
                                j * 8,
                                reg_from_index(limbs, (j+i) % limbs),
                                reg_from_index(limbs, (j+i+1) % limbs)));
            if i == limbs-1 {
                if is_spill(limbs, (j+i) % limbs) { asm_string = final_swap(asm_string, (j+i) % limbs, limbs, spill, a);}
            } else {
                if is_spill(limbs, (j+i) % limbs) { asm_string = spill_swap(asm_string, (j+i) % limbs, limbs, spill); }
            }

        }
        asm_string = format!("{}{}", asm_string, format!("
                                mulxq {}($1), %rax, {2}
                                mov $2, %rbx
                                adcxq %rax, {}
                                adoxq %rcx, {2}
                                adcxq %rbx, {2}",
                                (limbs-1)*8,
                                reg_from_index(limbs, (i+limbs-1) % limbs),
                                reg_from_index(limbs, i % limbs)));
        if i == limbs-1 {
            if is_spill(limbs, (i+limbs-1) % limbs) { asm_string = final_swap(asm_string, (i+limbs-1) % limbs, limbs, spill, a); }
            if is_spill(limbs, i % limbs) { asm_string = final_swap(asm_string, i % limbs, limbs, spill, a); }
        } else {
            if is_spill(limbs, (i+limbs-1) % limbs) { asm_string = spill_swap(asm_string, (i+limbs-1) % limbs, limbs, spill); }
            if is_spill(limbs, i % limbs) { asm_string = spill_swap(asm_string, i % limbs, limbs, spill); }
        }
    }
    for i in 0..limbs {
        if !is_spill(limbs, i) {
            asm_string = format!("{}{}", asm_string, format!("
                                movq %r{}, {}($0)", 8+(i % limbs), i*8));
        }
    }
    format!("{}{}", asm_string, "\"")
}

fn reg_from_index (limbs: usize, index: usize) -> String {
    let index = get_index(limbs, index);
    if index < 8 {
        format!("%r{}", index+8)
    } else {
        match index {
            8 => String::from("%rsi"),
            9 => String::from("%rdi"),
            _ => panic!("More than 10 registers is not supported")
        }
    }
}

fn is_spill(limbs: usize, index: usize) -> bool {
    let half = 1 + (MAX_REGS / 2);
    if limbs <= MAX_REGS { false } else if limbs <= (MAX_REGS+3) {
        if limbs == (MAX_REGS+1) {
            index % half == 0
        } else if limbs == (MAX_REGS+2) {
            (index % half == 0) | (index % half == 1)
        } else {
            (index % (half + 1) == 0) | (index % (half + 1) == 1) | (index % (half + 1) == 2)
        }
    } else { true }
}

fn get_index(limbs: usize, index: usize) -> usize {
    let half = 1 + (MAX_REGS / 2);
    if limbs <= MAX_REGS { index } else if limbs <= (MAX_REGS+3) {
        if limbs == (MAX_REGS+1) {
            if is_spill(limbs, index) { index % half } else {
                if index > half { index - 1 } else { index }
            }
        } else if limbs == (MAX_REGS+2) {
            if is_spill(limbs, index) { index % half } else {
                if index > half { index - 2 } else { index }
            }
        } else {
            if is_spill(limbs, index) { index % (half + 1) } else {
                if index > (half + 1) { index - 3 } else { index }
            }
        }
    } else { index % MAX_REGS }
}

fn get_spill_index(limbs: usize, index: usize) -> usize {
    let half = 1 + (MAX_REGS / 2);
    if limbs <= MAX_REGS {
        panic!("no spill for {} limbs", limbs);
    } else if limbs <= (MAX_REGS+3) {
        if limbs == (MAX_REGS+1) {
            if index >= half { 1 } else { index }
        } else if limbs == (MAX_REGS+2) {
            if index >= half { index - half + 2 } else { index }
        } else {
            if index >= (half + 1) { index - (half + 1) + 3 } else { index }
        }
    } else if index >= MAX_REGS {
        (index - MAX_REGS) + (limbs % MAX_REGS)
    } else {
        index
    }
}

fn swap_spill_index(limbs: usize, index: usize) -> usize {
    let half = 1 + (MAX_REGS / 2);
    if limbs <= MAX_REGS {
        panic!("no spill for {} limbs", limbs);
    } else if limbs <= (MAX_REGS+3) {
        if limbs == (MAX_REGS+1) {
            if index >= half { 0 } else { 1 }
        } else if limbs == (MAX_REGS+2) {
            if index >= half { index - half } else { index + 2 }
        } else {
            if index >= (half + 1) { index - (half + 1) } else { index + 3 }
        }
    } else if index >= MAX_REGS {
        index - MAX_REGS
    } else {
        index  + (limbs % MAX_REGS)
    }
}

fn spill_swap (asm_string: String, index: usize, limbs: usize, spill: &str) -> String {
    format!("{}{}", asm_string, format!("
                                mov {r}, {}({spill})
                                mov {}({spill}), {r}",
                                8 * get_spill_index(limbs, index),
                                8 * swap_spill_index(limbs, index),
                                r = reg_from_index(limbs, index), spill=spill))
    }

fn final_swap (asm_string: String, index: usize, limbs: usize, spill: &str, a: &str) -> String {
    if get_spill_index(limbs, index) < swap_spill_index(limbs, index) {
        format!("{}{}", asm_string, format!("
                                mov {r}, {}({a})
                                mov {}({spill}), {r}",
                                8 * index,
                                8 * swap_spill_index(limbs, index),
                                r = reg_from_index(limbs, index), a=a, spill=spill))
    } else {
        format!("{}{}", asm_string, format!("
                                mov {r}, {}({a})",
                                8 * index,r = reg_from_index(limbs, index), a=a))
    }
}
