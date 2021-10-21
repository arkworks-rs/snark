use rand::Rng;
use rayon::prelude::*;

use algebra::{
    groups::Group, msm::VariableBaseMSM, AffineCurve, Field, PairingEngine, PrimeField,
    ProjectiveCurve, UniformRand,
};

use crate::groth16::{r1cs_to_qap::R1CStoQAP, Parameters, Proof};

use r1cs_core::{
    ConstraintSynthesizer, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable,
};

use smallvec::SmallVec;

use std::{
    ops::{AddAssign, MulAssign, SubAssign},
    sync::Arc,
};

type CoeffVec<T> = SmallVec<[T; 2]>;

#[inline]
fn eval<E: PairingEngine>(
    lc: &LinearCombination<E::Fr>,
    constraints: &mut [CoeffVec<(E::Fr, Index)>],
    input_assignment: &[E::Fr],
    aux_assignment: &[E::Fr],
    this_constraint: usize,
) -> E::Fr {
    let mut acc = E::Fr::zero();

    for &(index, coeff) in lc.as_ref() {
        let mut tmp;

        match index.get_unchecked() {
            Index::Input(i) => {
                constraints[this_constraint].push((coeff, Index::Input(i)));
                tmp = input_assignment[i];
            },
            Index::Aux(i) => {
                constraints[this_constraint].push((coeff, Index::Aux(i)));
                tmp = aux_assignment[i];
            },
        }

        if coeff.is_one() {
            acc.add_assign(&tmp);
        } else {
            tmp.mul_assign(&coeff);
            acc.add_assign(&tmp);
        }
    }

    acc
}

pub struct ProvingAssignment<E: PairingEngine> {
    // Constraints
    pub(crate) at: Vec<CoeffVec<(E::Fr, Index)>>,
    pub(crate) bt: Vec<CoeffVec<(E::Fr, Index)>>,
    pub(crate) ct: Vec<CoeffVec<(E::Fr, Index)>>,

    // Evaluations of A and C polynomials
    pub(crate) a: Vec<E::Fr>,
    pub(crate) b: Vec<E::Fr>,
    pub(crate) c: Vec<E::Fr>,

    // Assignments of variables
    pub(crate) input_assignment: Vec<E::Fr>,
    pub(crate) aux_assignment:   Vec<E::Fr>,
    pub(crate) num_inputs:       usize,
    pub(crate) num_aux:          usize,
    pub(crate) num_constraints:  usize,
}

impl<E: PairingEngine> ProvingAssignment<E> {
    pub fn which_is_unsatisfied(&self) -> Option<usize> {
        for (i, ((a_i, b_i), c_i)) in (self.a.iter().zip(self.b.iter()))
            .zip(self.c.iter())
            .enumerate()
        {
            if *a_i * b_i != *c_i {
                return Some(i);
            }
        }
        None
    }
}

impl<E: PairingEngine> ConstraintSystem<E::Fr> for ProvingAssignment<E> {
    type Root = Self;

    #[inline]
    fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let index = self.num_aux;
        self.num_aux += 1;

        self.aux_assignment.push(f()?);
        Ok(Variable::new_unchecked(Index::Aux(index)))
    }

    #[inline]
    fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let index = self.num_inputs;
        self.num_inputs += 1;

        self.input_assignment.push(f()?);
        Ok(Variable::new_unchecked(Index::Input(index)))
    }

    #[inline]
    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E::Fr>) -> LinearCombination<E::Fr>,
        LB: FnOnce(LinearCombination<E::Fr>) -> LinearCombination<E::Fr>,
        LC: FnOnce(LinearCombination<E::Fr>) -> LinearCombination<E::Fr>,
    {
        self.at.push(CoeffVec::new());
        self.bt.push(CoeffVec::new());
        self.ct.push(CoeffVec::new());

        self.a.push(eval::<E>(
            &a(LinearCombination::zero()),
            &mut self.at,
            &self.input_assignment,
            &self.aux_assignment,
            self.num_constraints,
        ));
        self.b.push(eval::<E>(
            &b(LinearCombination::zero()),
            &mut self.bt,
            &self.input_assignment,
            &self.aux_assignment,
            self.num_constraints,
        ));
        self.c.push(eval::<E>(
            &c(LinearCombination::zero()),
            &mut self.ct,
            &self.input_assignment,
            &self.aux_assignment,
            self.num_constraints,
        ));

        self.num_constraints += 1;
    }

    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn pop_namespace(&mut self) {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }

    fn num_constraints(&self) -> usize {
        self.a.len()
    }
}

pub fn create_random_proof<E, C, R>(
    circuit: C,
    params: &Parameters<E>,
    rng: &mut R,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    let d1 = E::Fr::zero();
    let d2 = E::Fr::zero();
    let d3 = E::Fr::zero();
    let r = E::Fr::rand(rng);
    let s = E::Fr::rand(rng);

    create_proof::<E, C>(circuit, params, d1, d2, d3, r, s)
}

pub fn create_proof<E, C>(
    circuit: C,
    params: &Parameters<E>,
    d1: E::Fr,
    d2: E::Fr,
    d3: E::Fr,
    r: E::Fr,
    s: E::Fr,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    let prover_time = start_timer!(|| "Prover");
    let mut prover = ProvingAssignment {
        at:               vec![],
        bt:               vec![],
        ct:               vec![],
        a:                vec![],
        b:                vec![],
        c:                vec![],
        input_assignment: vec![],
        aux_assignment:   vec![],
        num_inputs:       0,
        num_aux:          0,
        num_constraints:  0,
    };

    // Allocate the "one" input variable
    prover.alloc_input(|| "", || Ok(E::Fr::one()))?;

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(&mut prover)?;
    end_timer!(synthesis_time);

    let witness_map_time = start_timer!(|| "R1CS to QAP witness map");
    let (full_input_assignment, h, _) = R1CStoQAP::witness_map::<E>(&prover, &d1, &d2, &d3)?;
    end_timer!(witness_map_time);

    let input_assignment = Arc::new(
        full_input_assignment[1..prover.num_inputs]
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let aux_assignment = Arc::new(
        full_input_assignment[prover.num_inputs..]
            .into_par_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );
    drop(full_input_assignment);

    let h_input_assignment = Arc::new(
        h[0..prover.num_inputs]
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );
    let h_aux_assignment = Arc::new(
        h[prover.num_inputs..]
            .into_par_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );
    drop(h);

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");
    let (a_inputs_source, a_aux_source) = params.get_a_query(prover.num_inputs)?;
    let a_inputs_acc = VariableBaseMSM::multi_scalar_mul(a_inputs_source, &input_assignment);
    let a_aux_acc = VariableBaseMSM::multi_scalar_mul(a_aux_source, &aux_assignment);

    let r_g1 = params.delta_g1.mul(r);

    let mut g_a = r_g1;
    g_a.add_assign(&params.get_a_query_full()?[0].into_projective());
    g_a.add_assign(&a_inputs_acc);
    g_a.add_assign(&a_aux_acc);
    g_a.add_assign(&params.alpha_g1.into());
    end_timer!(a_acc_time);

    // Compute B in G1
    let b_g1_acc_time = start_timer!(|| "Compute B in G1");

    let (b_inputs_source, b_aux_source) = params.get_b_g1_query(prover.num_inputs)?;
    let b_inputs_acc = VariableBaseMSM::multi_scalar_mul(b_inputs_source, &input_assignment);
    let b_aux_acc = VariableBaseMSM::multi_scalar_mul(b_aux_source, &aux_assignment);

    let s_g1 = params.delta_g1.mul(s.clone());

    let mut g1_b = s_g1;
    g1_b.add_assign(&params.get_b_g1_query_full()?[0].into_projective());
    g1_b.add_assign(&b_inputs_acc);
    g1_b.add_assign(&b_aux_acc);
    g1_b.add_assign(&params.beta_g1.into());
    end_timer!(b_g1_acc_time);

    // Compute B in G2
    let b_g2_acc_time = start_timer!(|| "Compute B in G2");

    let (b_inputs_source, b_aux_source) = params.get_b_g2_query(prover.num_inputs)?;
    let b_inputs_acc = VariableBaseMSM::multi_scalar_mul(b_inputs_source, &input_assignment);
    let b_aux_acc = VariableBaseMSM::multi_scalar_mul(b_aux_source, &aux_assignment);

    let s_g2 = params.delta_g2.mul(s.clone());

    let mut g2_b = s_g2;
    g2_b.add_assign(&params.get_b_g2_query_full()?[0].into_projective());
    g2_b.add_assign(&b_inputs_acc);
    g2_b.add_assign(&b_aux_acc);
    g2_b.add_assign(&params.beta_g2.into());
    end_timer!(b_g2_acc_time);

    // Compute C
    let c_acc_time = start_timer!(|| "Compute C");

    let (h_inputs_source, h_aux_source) = params.get_h_query(prover.num_inputs)?;
    let h_inputs_acc = VariableBaseMSM::multi_scalar_mul(h_inputs_source, &h_input_assignment);
    let h_aux_acc = VariableBaseMSM::multi_scalar_mul(h_aux_source, &h_aux_assignment);

    let l_aux_source = params.get_l_query_full()?;
    let l_aux_acc = VariableBaseMSM::multi_scalar_mul(l_aux_source, &aux_assignment);

    let s_g_a = g_a.clone().mul(&s);
    let r_g1_b = g1_b.clone().mul(&r);
    let r_s_delta_g1 = params.delta_g1.into_projective().mul(&r).mul(&s);

    let mut g_c = s_g_a;
    g_c.add_assign(&r_g1_b);
    g_c.sub_assign(&r_s_delta_g1);
    g_c.add_assign(&l_aux_acc);
    g_c.add_assign(&h_inputs_acc);
    g_c.add_assign(&h_aux_acc);
    end_timer!(c_acc_time);

    end_timer!(prover_time);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g2_b.into_affine(),
        c: g_c.into_affine(),
    })
}
