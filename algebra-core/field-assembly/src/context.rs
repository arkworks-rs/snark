use std::collections::HashMap;
use std::rc::Rc;

pub const REG_CLOBBER: [&'static str; 8] = ["r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15"];

#[derive(Clone)]
pub struct Context {
    ctx_string: Rc<String>,
    declarations: HashMap<String, Declare>,
    declaration_vec: Vec<Declare>,
    clobbers: Vec<String>,
}

#[derive(Clone)]
struct Declare {
    ty: String,
    var: String,
    pos: usize,
    token: String,
}

impl Context {
    pub fn new() -> Self {
        Context {
            ctx_string: Rc::new(String::new()),
            declarations: HashMap::new(),
            declaration_vec: Vec::new(),
            clobbers: Vec::new(),
        }
    }

    pub fn get_string(&mut self) -> String {
        Rc::make_mut(&mut self.ctx_string).to_string()
    }

    pub fn reset(&mut self) {
        self.declarations.clear();
        self.declaration_vec.clear();
        self.clobbers.clear();
    }

    pub fn get(self, id: &str) -> String {
        self.declarations
            .get(&id.to_string())
            .unwrap()
            .token
            .clone()
    }

    pub fn try_get(self, id: &str, fallback_id: &str) -> String {
        match self.declarations.get(&id.to_string()) {
            Some(dec) => dec.token.clone(),
            None => self
                .declarations
                .get(&fallback_id.to_string())
                .unwrap()
                .token
                .clone(),
        }
    }

    pub fn add_declaration(&mut self, id: &str, ty: &str, var: &str) {
        self.declarations.insert(
            id.to_string(),
            Declare {
                ty: ty.to_string(),
                var: var.to_string(),
                pos: self.declarations.len(),
                token: format!("${}", self.declarations.len()),
            },
        );
        self.declaration_vec.push(Declare {
            ty: ty.to_string(),
            var: var.to_string(),
            pos: self.declaration_vec.len(),
            token: format!("${}", self.declaration_vec.len()),
        });
    }

    pub fn add_limb(&mut self, limb: usize) {
        self.ctx_string = Rc::new(format!(
            "{}{}",
            Rc::clone(&self.ctx_string),
            format!(
                "
                {} => {{",
                limb
            )
        ));
    }

    pub fn add_buffer(&mut self, extra_reg: usize) {
        self.ctx_string = Rc::new(format!(
            "{}{}",
            Rc::clone(&self.ctx_string),
            format!(
                "
                    let mut spill_buffer = MaybeUninit::<[u64; {}]>::uninit();",
                extra_reg
            )
        ));
    }

    pub fn add_asm(&mut self, ctx_string: String) {
        self.ctx_string = Rc::new(format!(
            "{}{}",
            Rc::clone(&self.ctx_string),
            format!(
                "
                    unsafe {{
                        asm!({}
                            :
                            :",
                ctx_string
            )
        ));
    }

    pub fn add_clobber_from_vec(&mut self, clobbers: Vec<&str>) {
        for clobber in clobbers {
            self.clobbers.push(format!(" \"{}\"", clobber));
        }
    }

    pub fn add_clobber(&mut self, clobber: &str) {
        self.clobbers.push(format!(" \"{}\"", clobber));
    }

    pub fn build(&mut self) {
        for i in 0..self.declarations.len() {
            let dec = &self.declaration_vec[i];
            let last = i == self.declarations.len() - 1;
            self.ctx_string = Rc::new(format!(
                "{}{}",
                Rc::clone(&self.ctx_string),
                format!(
                    "
                            \"{}\"({}){}      // {}",
                    dec.ty,
                    dec.var,
                    if last { "" } else { "," },
                    dec.pos
                )
            ));
        }
        self.ctx_string = Rc::new(format!(
            "{}{}",
            Rc::clone(&self.ctx_string),
            format!(
                "
                            : {}
                        );
                    }}
                }}",
                self.clobbers.join(",")
            )
        ));
    }

    pub fn end(&mut self, num_limbs: usize) {
        self.ctx_string = Rc::new(format!("{}{}", Rc::clone(&self.ctx_string), format!("
            x => panic!(\"asm_mul (no-carry): number of limbs supported is 2 up to {}. You had {{}}.\", x)
        }};
    }}
}}
",
        num_limbs)));
    }
}
