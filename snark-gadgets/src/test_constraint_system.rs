use algebra::Field;
use r1cs_core::{ConstraintSystem, Index, LinearCombination, SynthesisError, Variable};

use radix_trie::Trie;

#[derive(Debug)]
enum NamedObject {
    Constraint(usize),
    Var(Variable),
    Namespace,
}

/// Constraint system for testing purposes.
pub struct TestConstraintSystem<ConstraintF: Field> {
    named_objects:     Trie<String, NamedObject>,
    current_namespace: Vec<String>,
    pub constraints: Vec<(
        LinearCombination<ConstraintF>,
        LinearCombination<ConstraintF>,
        LinearCombination<ConstraintF>,
        String,
    )>,
    inputs:            Vec<(ConstraintF, String)>,
    aux:               Vec<(ConstraintF, String)>,
}

impl<ConstraintF: Field> TestConstraintSystem<ConstraintF> {
    fn eval_lc(
        terms: &[(Variable, ConstraintF)],
        inputs: &[(ConstraintF, String)],
        aux: &[(ConstraintF, String)],
    ) -> ConstraintF {
        let mut acc = ConstraintF::zero();

        for &(var, ref coeff) in terms {
            let mut tmp = match var.get_unchecked() {
                Index::Input(index) => inputs[index].0,
                Index::Aux(index) => aux[index].0,
            };

            tmp.mul_assign(&coeff);
            acc.add_assign(&tmp);
        }

        acc
    }
}

impl<ConstraintF: Field> TestConstraintSystem<ConstraintF> {
    pub fn new() -> TestConstraintSystem<ConstraintF> {
        let mut map = Trie::new();
        map.insert(
            "ONE".into(),
            NamedObject::Var(TestConstraintSystem::<ConstraintF>::one()),
        );

        TestConstraintSystem {
            named_objects:     map,
            current_namespace: vec![],
            constraints:       vec![],
            inputs:            vec![(ConstraintF::one(), "ONE".into())],
            aux:               vec![],
        }
    }

    pub fn print_named_objects(&self) {
        for &(_, _, _, ref name) in &self.constraints {
            println!("{}", name);
        }
    }

    pub fn which_is_unsatisfied(&self) -> Option<&str> {
        for &(ref a, ref b, ref c, ref path) in &self.constraints {
            let mut a = Self::eval_lc(a.as_ref(), &self.inputs, &self.aux);
            let b = Self::eval_lc(b.as_ref(), &self.inputs, &self.aux);
            let c = Self::eval_lc(c.as_ref(), &self.inputs, &self.aux);

            a.mul_assign(&b);

            if a != c {
                return Some(&*path);
            }
        }

        None
    }

    pub fn is_satisfied(&self) -> bool {
        self.which_is_unsatisfied().is_none()
    }

    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }

    pub fn set(&mut self, path: &str, to: ConstraintF) {
        match self.named_objects.get(path) {
            Some(&NamedObject::Var(ref v)) => match v.get_unchecked() {
                Index::Input(index) => self.inputs[index].0 = to,
                Index::Aux(index) => self.aux[index].0 = to,
            },
            Some(e) => panic!(
                "tried to set path `{}` to value, but `{:?}` already exists there.",
                path, e
            ),
            _ => panic!("no variable exists at path: {}", path),
        }
    }

    pub fn get(&mut self, path: &str) -> ConstraintF {
        match self.named_objects.get(path) {
            Some(&NamedObject::Var(ref v)) => match v.get_unchecked() {
                Index::Input(index) => self.inputs[index].0,
                Index::Aux(index) => self.aux[index].0,
            },
            Some(e) => panic!(
                "tried to get value of path `{}`, but `{:?}` exists there (not a variable)",
                path, e
            ),
            _ => panic!("no variable exists at path: {}", path),
        }
    }

    fn set_named_obj(&mut self, path: String, to: NamedObject) {
        if self.named_objects.get(&path).is_some() {
            panic!("tried to create object at existing path: {}", path);
        }

        self.named_objects.insert(path, to);
    }
}

fn compute_path(ns: &[String], this: String) -> String {
    if this.chars().any(|a| a == '/') {
        panic!("'/' is not allowed in names");
    }

    let mut name = String::new();

    let mut needs_separation = false;
    for ns in ns.iter().chain(Some(&this).into_iter()) {
        if needs_separation {
            name += "/";
        }

        name += ns;
        needs_separation = true;
    }

    name
}

impl<ConstraintF: Field> ConstraintSystem<ConstraintF> for TestConstraintSystem<ConstraintF> {
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let index = self.aux.len();
        let path = compute_path(&self.current_namespace, annotation().into());
        self.aux.push((f()?, path.clone()));
        let var = Variable::new_unchecked(Index::Aux(index));
        self.set_named_obj(path, NamedObject::Var(var));

        Ok(var)
    }

    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let index = self.inputs.len();
        let path = compute_path(&self.current_namespace, annotation().into());
        self.inputs.push((f()?, path.clone()));
        let var = Variable::new_unchecked(Index::Input(index));
        self.set_named_obj(path, NamedObject::Var(var));

        Ok(var)
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LB: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LC: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
    {
        let path = compute_path(&self.current_namespace, annotation().into());
        let index = self.constraints.len();
        self.set_named_obj(path.clone(), NamedObject::Constraint(index));

        let mut a = a(LinearCombination::zero());
        let mut b = b(LinearCombination::zero());
        let mut c = c(LinearCombination::zero());
        a.0.shrink_to_fit();
        b.0.shrink_to_fit();
        c.0.shrink_to_fit();

        self.constraints.push((a, b, c, path));
    }

    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        let name = name_fn().into();
        let path = compute_path(&self.current_namespace, name.clone());
        self.set_named_obj(path.clone(), NamedObject::Namespace);
        self.current_namespace.push(name);
    }

    fn pop_namespace(&mut self) {
        assert!(self.current_namespace.pop().is_some());
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }

    fn num_constraints(&self) -> usize {
        self.constraints.len()
    }
}
