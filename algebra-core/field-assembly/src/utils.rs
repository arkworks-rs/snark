pub const RAX: &'static str = "%rax";
pub const RBX: &'static str = "%rbx";
pub const RCX: &'static str = "%rcx";
pub const RDX: &'static str = "%rdx";
pub const RDI: &'static str = "%rdi";
pub const RSI: &'static str = "%rsi";
pub const R: [&'static str; 8] = ["%r8", "%r9", "%r10", "%r11", "%r12", "%r13", "%r14", "%r15"];

macro_rules! reg {
    ($a_reg:ident, $a:ident, $range:expr) => {
        paste::item! {
            let mut $a_reg = Vec::new();
            let mut [<$a_reg _1>] = Vec::new();
            for i in 0..$range {
                [<$a_reg _1>].push(format!("{}({})", i * 8, $a));
            }
            for i in 0..$range {
                $a_reg.push(&*[<$a_reg _1>][i]);
            }
        }
    };
}
