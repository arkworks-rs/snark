use proc_macro::TokenStream;
use quote::quote;

pub fn define_intrinsics() -> TokenStream {
    (quote! {
        {
            let mut begin = || {
                asm_string.borrow_mut().push_str("\"");
            };

            let mut end = || {
                asm_string.borrow_mut().push_str("
                                        \"");
            };

            let mut comment = | comment: &str | {
                asm_string.borrow_mut().push_str(&format!("         // {}", comment));
            };

            let mut mulxq = | a: &str, b: &str, c: &str | {
                asm_string.borrow_mut().push_str(&format!("
                                        mulxq {}, {}, {}", a, b, c));
            };

            let mut adcxq = | a: &str, b: &str| {
                asm_string.borrow_mut().push_str(&format!("
                                        adcxq {}, {}", a, b));
            };

            let mut adoxq = | a: &str, b: &str | {
                asm_string.borrow_mut().push_str(&format!("
                                        adoxq {}, {}", a, b));
            };

            let mut movq = | a: &str, b: &str | {
                asm_string.borrow_mut().push_str(&format!("
                                        movq {}, {}", a, b));
            };

            let mut xorq = | a: &str, b: &str | {
                asm_string.borrow_mut().push_str(&format!("
                                        xorq {}, {}", a, b));
            };
        }
    })
    .into()
}
