use crate::assembler::*;
// TODO: Replace assembler with macro to write movq instead of asm.movq

// Computes [rcx, r(8+limbs), ..., r8] =  0($a) * [r(8+limbs), ..., r8]
pub fn mul_1 (asm: &mut Assembler, a: &str, b: &Vec<&str>, zero: &str) {
    asm.movq(a, RDX);
    asm.mulxq(b[0], R[0], R[1]);
    for j in 1..asm.limbs-1 {
        asm.mulxq(b[j], RAX, R[((j + 1) % asm.limbs)]);
        asm.adcxq(RAX, R[j]);
    }
    asm.mulxq(b[asm.limbs-1], RAX, RCX);
    asm.movq(zero, RBX);
    asm.adcxq(RAX, R[asm.limbs-1]);
    asm.adcxq(RBX, RCX);
}


pub fn mul_add_1 (asm: &mut Assembler, a: &Vec<&str>, b: &Vec<&str>, zero: &str, i: usize) {
    asm.movq(a[i], RDX);
    for j in 0..asm.limbs-1 {
        asm.mulxq(b[j], RAX, RBX);
        asm.adcxq(RAX, R[(j+i) % asm.limbs]);
        asm.adoxq(RBX, R[(j+i+1) % asm.limbs]);
    }
    asm.mulxq(b[asm.limbs-1], RAX, RCX);
    asm.movq(zero, RBX);
    asm.adcxq(RAX, R[(i+asm.limbs-1) % asm.limbs]);
    asm.adoxq(RBX, RCX);
    asm.adcxq(RBX, RCX);
}

pub fn mul_add_shift_1 (asm: &mut Assembler, a: &Vec<&str>, inverse: &str, zero: &str, i: usize) {
    asm.movq(inverse, RDX);
    asm.mulxq(R[i], RDX, RAX);
    asm.mulxq(a[0], RAX, RBX);
    asm.adcxq(R[i % asm.limbs], RAX);
    asm.adoxq(RBX, R[(i+1) % asm.limbs]);
    for j in 1..asm.limbs-1 {
        asm.mulxq(a[j], RAX, RBX);
        asm.adcxq(RAX, R[(j+i) % asm.limbs]);
        asm.adoxq(RBX, R[(j+i+1) % asm.limbs]);
    }
    asm.mulxq(a[asm.limbs-1], RAX, R[i % asm.limbs]);
    asm.movq(zero, RBX);
    asm.adcxq(RAX, R[(i+asm.limbs-1) % asm.limbs]);
    asm.adoxq(RCX, R[i % asm.limbs]);
    asm.adcxq(RBX, R[i % asm.limbs]);
}

// Computes [rcx, r(8+limbs), ..., r8] =  0($a) * [r(8+limbs), ..., r8]
pub fn mul_1_mov (asm: &mut Assembler, a: &str, b: &Vec<&str>) {
    asm.movq(a, RDX);
    asm.mulxq(b[0], R[0], R[1]);
    for j in 1..asm.limbs-1 {
        asm.mulxq(b[j], RAX, R[((j + 1) % asm.limbs)]);
        asm.adcxq(RAX, R[j]);
    }
    asm.mulxq(b[asm.limbs-1], RAX, RCX);
    asm.movq("$2", RBX);
    asm.adcxq(RAX, R[asm.limbs-1]);
    asm.adcxq(RBX, RCX);
}
