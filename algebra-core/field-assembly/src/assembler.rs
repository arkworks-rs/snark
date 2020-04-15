use std::rc::Rc;

pub const RAX: &'static str = "%rax";
pub const RBX: &'static str = "%rbx";
pub const RCX: &'static str = "%rcx";
pub const RDX: &'static str = "%rdx";
pub const RDI: &'static str = "%rdi";
pub const RSI: &'static str = "%rsi";
pub const R: [&'static str; 18] = ["%r8", "%r9", "%r10", "%r11", "%r12", "%r13", "%r14", "%r15",
                                    "%r16", "%r17", "%r18", "%r19", "%r20", "%r21", "%r22", "%r23", "%r24", "%r25"];

pub struct Assembler {
    pub limbs: usize,
    asm_string: Rc<String>,

}

impl<'a> Assembler {
    pub fn new (limbs: usize) -> Assembler {
        Assembler {
            limbs: limbs,
            asm_string: Rc::new(String::new()),
        }
    }

    pub fn get_asm_string (&mut self) -> String {
        Rc::make_mut(&mut self.asm_string).to_string()
    }

    pub fn begin (&mut self) {
        self.asm_string = Rc::new("\"".to_string());
    }
    pub fn end (&mut self) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), "
                                \"".to_string()));
    }

    pub fn comment (&mut self, comment: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("         // {}", comment)));
    }

    pub fn mulxq (&mut self, a: &str, b: &str, c: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("
                                    mulxq {}, {}, {}", a, b, c)));
    }

    pub fn adcxq (&mut self, a: &str, b: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("
                                    adcxq {}, {}", a, b)));
    }

    pub fn adoxq (&mut self, a: &str, b: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("
                                    adoxq {}, {}", a, b)));
    }

    pub fn movq (&mut self, a: &str, b: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("
                                    movq {}, {}", a, b)));
    }

    pub fn xorq (&mut self, a: &str, b: &str) {
        self.asm_string = Rc::new(format!("{}{}", Rc::clone(&self.asm_string), format!("
                                    xorq {}, {}", a, b)));
    }
}

macro_rules! generate_array {
    ($a_0:ident, $a_1:ident, $a:ident, $range:expr) => {
        let mut $a_0 = Vec::new();
        let mut $a_1 = Vec::new();
        for i in 0..$range {
            $a_0.push(format!("{}({})", i*8, $a));
        }
        for i in 0..$range {
            $a_1.push(&*$a_0[i]);
        }
    }
}
