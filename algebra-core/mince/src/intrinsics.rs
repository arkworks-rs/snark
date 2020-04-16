use proc_macro::TokenStream;
use quote::quote;

pub fn define_intrinsics() -> TokenStream {
    (quote! {
        {
            let mut begin = || {
                asm_string.replace_with(|_| "\"".to_string());
            };

            let mut end = || {
                asm_string.replace_with(|x| format!("{}{}", x, "
                                        \"".to_string()).clone());
            };

            let mut comment = | comment: &str | {
                asm_string.replace_with(|x| format!("{}{}", x, format!("         // {}", comment)).clone());
            };

            let mut mulxq = | a: &str, b: &str, c: &str | {
                asm_string.replace_with(|x| format!("{}{}", x, format!("
                                        mulxq {}, {}, {}", a, b, c)).clone());
            };

            let mut adcxq = | a: &str, b: &str| {
                asm_string.replace_with(|x| format!("{}{}", x, format!("
                                        adcxq {}, {}", a, b)));
            };

            let mut adoxq = | a: &str, b: &str | {
                asm_string.replace_with(|x| format!("{}{}", x, format!("
                                        adoxq {}, {}", a, b)).clone());
            };

            let mut movq = | a: &str, b: &str | {
                asm_string.replace_with(|x| format!("{}{}", x, format!("
                                        movq {}, {}", a, b)).clone());
            };

            let mut xorq = | a: &str, b: &str | {
                asm_string.replace_with(|x| format!("{}{}", x, format!("
                                        xorq {}, {}", a, b)).clone());
            };
        }
    }).into()
}
