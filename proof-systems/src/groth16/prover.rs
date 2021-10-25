use rand::Rng;
use rayon::prelude::*;

use algebra::msm::VariableBaseMSM;
use algebra::{
    groups::Group, AffineCurve, Field, PairingEngine, PrimeField, ProjectiveCurve, UniformRand,
};

use crate::groth16::{push_constraints, r1cs_to_qap::R1CStoQAP, Parameters, Proof};

use r1cs_core::{
    ConstraintSynthesizer, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable,
};

use std::{
    ops::{AddAssign, SubAssign},
    sync::Arc,
};

pub struct ProvingAssignment<E: PairingEngine> {
    // Constraints
    pub(crate) at: Vec<Vec<(E::Fr, Index)>>,
    pub(crate) bt: Vec<Vec<(E::Fr, Index)>>,
    pub(crate) ct: Vec<Vec<(E::Fr, Index)>>,

    // Assignments of variables
    pub(crate) input_assignment: Vec<E::Fr>,
    pub(crate) aux_assignment: Vec<E::Fr>,
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
        let index = self.aux_assignment.len();

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
        let index = self.input_assignment.len();

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
        let num_constraints = self.num_constraints();

        self.at.push(Vec::new());
        self.bt.push(Vec::new());
        self.ct.push(Vec::new());

        push_constraints(a(LinearCombination::zero()), &mut self.at, num_constraints);
        push_constraints(b(LinearCombination::zero()), &mut self.bt, num_constraints);
        push_constraints(c(LinearCombination::zero()), &mut self.ct, num_constraints);
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
        self.at.len()
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
    let r = E::Fr::rand(rng);
    let s = E::Fr::rand(rng);

    create_proof::<E, C>(circuit, params, r, s)
}

pub fn create_proof_no_zk<E, C>(
    circuit: C,
    params: &Parameters<E>,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    create_proof::<E, C>(circuit, params, E::Fr::zero(), E::Fr::zero())
}

pub fn create_proof<E, C>(
    circuit: C,
    params: &Parameters<E>,
    r: E::Fr,
    s: E::Fr,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    let prover_time = start_timer!(|| "Prover");
    let mut prover = ProvingAssignment {
        at: vec![],
        bt: vec![],
        ct: vec![],
        input_assignment: vec![],
        aux_assignment: vec![],
    };

    // Allocate the "one" input variable
    prover.alloc_input(|| "", || Ok(E::Fr::one()))?;

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(&mut prover)?;
    end_timer!(synthesis_time);

    let witness_map_time = start_timer!(|| "R1CS to QAP witness map");
    let h = R1CStoQAP::witness_map::<E>(&prover)?;
    end_timer!(witness_map_time);

    let input_assignment = Arc::new(
        prover.input_assignment[1..]
            .iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let aux_assignment = Arc::new(
        prover
            .aux_assignment
            .into_par_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let assignment = [&input_assignment[..], &aux_assignment[..]].concat();

    let h_assignment = h.into_par_iter().map(|s| s.into_repr()).collect::<Vec<_>>();

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");

    let a_query = params.get_a_query_full()?;
    let r_g1 = params.delta_g1.mul(r);
    let g_a = calculate_coeff(r_g1, a_query, &params.alpha_g1, &assignment)?;

    end_timer!(a_acc_time);

    // Compute B in G1 if needed
    let g1_b = if r != E::Fr::zero() {
        let b_g1_acc_time = start_timer!(|| "Compute B in G1");

        let s_g1 = params.delta_g1.mul(s);
        let b_query = params.get_b_g1_query_full()?;
        let g1_b = calculate_coeff(s_g1, b_query, &params.beta_g1, &assignment)?;

        end_timer!(b_g1_acc_time);
        g1_b
    } else {
        <E::G1Projective as ProjectiveCurve>::zero()
    };

    // Compute B in G2
    let b_g2_acc_time = start_timer!(|| "Compute B in G2");

    let b_query = params.get_b_g2_query_full()?;
    let s_g2 = params.delta_g2.mul(s);
    let g2_b = calculate_coeff(s_g2, b_query, &params.beta_g2, &assignment)?;

    end_timer!(b_g2_acc_time);

    // Compute C
    let c_acc_time = start_timer!(|| "Compute C");

    let h_query = params.get_h_query_full()?;
    let h_acc = VariableBaseMSM::multi_scalar_mul(&h_query, &h_assignment)?;

    let l_aux_source = params.get_l_query_full()?;
    let l_aux_acc = VariableBaseMSM::multi_scalar_mul(l_aux_source, &aux_assignment)?;

    let s_g_a = g_a.mul(&s);
    let r_g1_b = g1_b.mul(&r);
    let r_s_delta_g1 = params.delta_g1.into_projective().mul(&r).mul(&s);

    let mut g_c = s_g_a;
    g_c.add_assign(&r_g1_b);
    g_c.sub_assign(&r_s_delta_g1);
    g_c.add_assign(&l_aux_acc);
    g_c.add_assign(&h_acc);

    end_timer!(c_acc_time);

    end_timer!(prover_time);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g2_b.into_affine(),
        c: g_c.into_affine(),
    })
}

fn calculate_coeff<G: AffineCurve>(
    initial: G::Projective,
    query: &[G],
    vk_param: &G,
    assignment: &[<G::ScalarField as PrimeField>::BigInt],
) -> Result<G::Projective, SynthesisError> {
    let el = query[0];
    let acc = VariableBaseMSM::multi_scalar_mul(&query[1..], assignment)?;

    let mut res = initial;
    res.add_assign_mixed(&el);
    res += &acc;
    res.add_assign_mixed(vk_param);

    Ok(res)
}
