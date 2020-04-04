use rand::Rng;

use algebra_core::{
    msm::VariableBaseMSM, AffineCurve, One, PairingEngine, PrimeField, ProjectiveCurve,
    UniformRand, Zero,
};

use crate::{push_constraints, r1cs_to_qap::R1CStoQAP, Parameters, Proof, String, Vec};

use r1cs_core::{
    ConstraintSynthesizer, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable,
};

use ff_fft::cfg_into_iter;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct ProvingAssignment<E: PairingEngine> {
    // Constraints
    pub(crate) at: Vec<Vec<(E::Fr, Index)>>,
    pub(crate) bt: Vec<Vec<(E::Fr, Index)>>,
    pub(crate) ct: Vec<Vec<(E::Fr, Index)>>,

    // Assignments of variables
    pub(crate) input_assignment: Vec<E::Fr>,
    pub(crate) aux_assignment:   Vec<E::Fr>,
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
        at:               vec![],
        bt:               vec![],
        ct:               vec![],
        input_assignment: vec![],
        aux_assignment:   vec![],
    };

    // Allocate the "one" input variable
    prover.alloc_input(|| "", || Ok(E::Fr::one()))?;

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(&mut prover)?;
    end_timer!(synthesis_time);

    let witness_map_time = start_timer!(|| "R1CS to QAP witness map");
    let (full_input_assignment, h, _) = R1CStoQAP::witness_map::<E>(&prover)?;
    end_timer!(witness_map_time);

    let num_inputs = prover.input_assignment.len();

    let input_assignment = cfg_into_iter!(full_input_assignment[1..num_inputs])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    let aux_assignment = cfg_into_iter!(full_input_assignment[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    drop(full_input_assignment);

    let h_input_assignment = cfg_into_iter!(h[0..num_inputs])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    let h_aux_assignment = cfg_into_iter!(h[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    drop(h);

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");
    let (a_inputs_source, a_aux_source) = params.get_a_query(num_inputs)?;
    let a_query = params.get_a_query_full()?[0];
    let r_g1 = params.delta_g1.mul(r);

    let g_a = calculate_coeff(
        r_g1,
        a_inputs_source,
        a_aux_source,
        &input_assignment,
        &aux_assignment,
        a_query,
        params.vk.alpha_g1,
    );

    end_timer!(a_acc_time);

    // Compute B in G1 if needed
    let g1_b = if r != E::Fr::zero() {
        let b_g1_acc_time = start_timer!(|| "Compute B in G1");
        let s_g1 = params.delta_g1.mul(s);
        let b_query = params.get_b_g1_query_full()?[0];
        let (b_inputs_source, b_aux_source) = params.get_b_g1_query(num_inputs)?;

        let g1_b = calculate_coeff(
            s_g1,
            b_inputs_source,
            b_aux_source,
            &input_assignment,
            &aux_assignment,
            b_query,
            params.beta_g1,
        );

        end_timer!(b_g1_acc_time);

        g1_b
    } else {
        E::G1Projective::zero()
    };

    // Compute B in G2
    let b_g2_acc_time = start_timer!(|| "Compute B in G2");
    let (b_inputs_source, b_aux_source) = params.get_b_g2_query(num_inputs)?;
    let b_query = params.get_b_g2_query_full()?[0];
    let s_g2 = params.vk.delta_g2.mul(s);
    let g2_b = calculate_coeff(
        s_g2,
        b_inputs_source,
        b_aux_source,
        &input_assignment,
        &aux_assignment,
        b_query,
        params.vk.beta_g2,
    );

    end_timer!(b_g2_acc_time);

    // Compute C
    let c_acc_time = start_timer!(|| "Compute C");

    let (h_inputs_source, h_aux_source) = params.get_h_query(num_inputs)?;
    let h_inputs_acc = VariableBaseMSM::multi_scalar_mul(h_inputs_source, &h_input_assignment);
    let h_aux_acc = VariableBaseMSM::multi_scalar_mul(h_aux_source, &h_aux_assignment);

    let l_aux_source = params.get_l_query_full()?;
    let l_aux_acc = VariableBaseMSM::multi_scalar_mul(l_aux_source, &aux_assignment);

    let s_g_a = g_a.mul(s);
    let r_g1_b = g1_b.mul(r);
    let r_s_delta_g1 = params.delta_g1.into_projective().mul(r).mul(s);

    let mut g_c = s_g_a;
    g_c += &r_g1_b;
    g_c -= &r_s_delta_g1;
    g_c += &l_aux_acc;
    g_c += &h_inputs_acc;
    g_c += &h_aux_acc;
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
    inputs_source: &[G],
    aux_source: &[G],
    input_assignment: &[<G::ScalarField as PrimeField>::BigInt],
    aux_assignment: &[<G::ScalarField as PrimeField>::BigInt],
    query: G,
    vk_param: G,
) -> G::Projective {
    let mut res = initial;

    let inputs_acc = VariableBaseMSM::multi_scalar_mul(inputs_source, input_assignment);
    let aux_acc = VariableBaseMSM::multi_scalar_mul(aux_source, aux_assignment);

    res.add_assign_mixed(&query);
    res += &inputs_acc;
    res += &aux_acc;
    res.add_assign_mixed(&vk_param);

    res
}
