use proc_macro::TokenStream;
use quote::quote;

pub fn define_arithmetic() -> TokenStream {
    (quote! {
        {
            macro_rules! mul_1 {
                ($a:expr, $b:ident, $zero:ident, $limbs:expr) => {
                    movq($a, RDX);
                    mulxq($b[0], R[0], R[1]);
                    for j in 1..$limbs-1 {
                        mulxq($b[j], RAX, R[((j + 1) % $limbs)]);
                        adcxq(RAX, R[j]);
                    }
                    mulxq($b[$limbs-1], RAX, RCX);
                    movq($zero, RBX);
                    adcxq(RAX, R[$limbs-1]);
                    adcxq(RBX, RCX);
                }
            }

            macro_rules! mul_add_1 {
                ($a:ident, $b:ident, $zero:ident, $i:ident, $limbs:expr) => {
                    movq($a[$i], RDX);
                    for j in 0..$limbs-1 {
                        mulxq($b[j], RAX, RBX);
                        adcxq(RAX, R[(j+$i) % $limbs]);
                        adoxq(RBX, R[(j+$i+1) % $limbs]);
                    }
                    mulxq($b[$limbs-1], RAX, RCX);
                    movq($zero, RBX);
                    adcxq(RAX, R[($i+$limbs-1) % $limbs]);
                    adoxq(RBX, RCX);
                    adcxq(RBX, RCX);
                }
            }

            macro_rules! mul_add_shift_1 {
                ($a:ident, $mod_prime:ident, $zero:ident, $i:ident, $limbs:expr) => {
                    movq($mod_prime, RDX);
                    mulxq(R[$i], RDX, RAX);
                    mulxq($a[0], RAX, RBX);
                    adcxq(R[$i % $limbs], RAX);
                    adoxq(RBX, R[($i+1) % $limbs]);
                    for j in 1..$limbs-1 {
                        mulxq($a[j], RAX, RBX);
                        adcxq(RAX, R[(j+$i) % $limbs]);
                        adoxq(RBX, R[(j+$i+1) % $limbs]);
                    }
                    mulxq($a[$limbs-1], RAX, R[$i % $limbs]);
                    movq($zero, RBX);
                    adcxq(RAX, R[($i+$limbs-1) % $limbs]);
                    adoxq(RCX, R[$i % $limbs]);
                    adcxq(RBX, R[$i % $limbs]);
                }
            }
        }
    })
    .into()
}
