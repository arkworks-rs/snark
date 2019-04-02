use algebra::{
    fft::domain::{EvaluationDomain, Scalar},
    msm::FixedBaseMSM,
    AffineCurve, Field, PairingEngine, PrimeField, ProjectiveCurve,
};

use rand::Rng;
use rayon::prelude::*;

use super::{Parameters, VerifyingKey};

use crate::{Circuit, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable};

use crate::gm17::r1cs_to_sap::R1CStoSAP;

/// Generates a random common reference string for
/// a circuit.
pub fn generate_random_parameters<E, C, R>(
    circuit: C,
    rng: &mut R,
) -> Result<Parameters<E>, SynthesisError>
where
    E: PairingEngine,
    C: Circuit<E>,
    R: Rng,
{
    let alpha = rng.gen();
    let beta = rng.gen();
    let gamma = E::Fr::one();
    let g = rng.gen();
    let h = rng.gen();

    generate_parameters::<E, C, R>(circuit, alpha, beta, gamma, g, h, rng)
}

/// This is our assembly structure that we'll use to synthesize the
/// circuit into a SAP.
pub struct KeypairAssembly<E: PairingEngine> {
    pub(crate) num_inputs:      usize,
    pub(crate) num_aux:         usize,
    pub(crate) num_constraints: usize,
    pub(crate) at:              Vec<Vec<(E::Fr, Index)>>,
    pub(crate) bt:              Vec<Vec<(E::Fr, Index)>>,
    pub(crate) ct:              Vec<Vec<(E::Fr, Index)>>,
}

impl<E: PairingEngine> ConstraintSystem<E> for KeypairAssembly<E> {
    type Root = Self;

    #[inline]
    fn alloc<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't invoke the
        // function for obtaining one.

        let index = self.num_aux;
        self.num_aux += 1;

        Ok(Variable::new_unchecked(Index::Aux(index)))
    }

    #[inline]
    fn alloc_input<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't invoke the
        // function for obtaining one.

        let index = self.num_inputs;
        self.num_inputs += 1;

        Ok(Variable::new_unchecked(Index::Input(index)))
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        fn eval<E: PairingEngine>(
            l: LinearCombination<E>,
            constraints: &mut [Vec<(E::Fr, Index)>],
            this_constraint: usize,
        ) {
            for (var, coeff) in l.as_ref() {
                match var.get_unchecked() {
                    Index::Input(i) => constraints[this_constraint].push((*coeff, Index::Input(i))),
                    Index::Aux(i) => constraints[this_constraint].push((*coeff, Index::Aux(i))),
                }
            }
        }

        self.at.push(vec![]);
        self.bt.push(vec![]);
        self.ct.push(vec![]);

        eval(
            a(LinearCombination::zero()),
            &mut self.at,
            self.num_constraints,
        );
        eval(
            b(LinearCombination::zero()),
            &mut self.bt,
            self.num_constraints,
        );
        eval(
            c(LinearCombination::zero()),
            &mut self.ct,
            self.num_constraints,
        );

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
        self.num_constraints
    }
}

/// Create parameters for a circuit, given some toxic waste.
pub fn generate_parameters<E, C, R>(
    circuit: C,
    alpha: E::Fr,
    beta: E::Fr,
    gamma: E::Fr,
    g: E::G1Projective,
    h: E::G2Projective,
    rng: &mut R,
) -> Result<Parameters<E>, SynthesisError>
where
    E: PairingEngine,
    C: Circuit<E>,
    R: Rng,
{
    let mut assembly = KeypairAssembly {
        num_inputs:      0,
        num_aux:         0,
        num_constraints: 0,
        at:              vec![],
        bt:              vec![],
        ct:              vec![],
    };

    // Allocate the "one" input variable
    assembly.alloc_input(|| "", || Ok(E::Fr::one()))?;

    // Synthesize the circuit.
    let synthesis_time = timer_start!(|| "Circuit synthesis");
    circuit.synthesize(&mut assembly)?;
    timer_end!(synthesis_time);

    ///////////////////////////////////////////////////////////////////////////
    let domain_time = timer_start!(|| "Constructing evaluation domain");

    let domain = vec![
        Scalar::<E>(E::Fr::zero());
        2 * assembly.num_constraints + 2 * assembly.num_inputs - 1
    ];

    let domain = EvaluationDomain::<E, _>::from_coeffs(domain)
        .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
    let mut t = rng.gen();

    while domain.z(&t).is_zero() {
        t = rng.gen();
    }

    timer_end!(domain_time);
    ///////////////////////////////////////////////////////////////////////////

    let reduction_time = timer_start!(|| "R1CS to SAP Instance Map with Evaluation");
    let (a, c, zt, sap_num_variables, m_raw) =
        R1CStoSAP::instance_map_with_evaluation(&assembly, &t)?;
    timer_end!(reduction_time);

    // Compute query densities
    let non_zero_a = (0..sap_num_variables)
        .into_par_iter()
        .map(|i| (!a[i].is_zero()) as usize)
        .sum();
    let scalar_bits = E::Fr::size_in_bits();

    // Compute G window table
    let g_window_time = timer_start!(|| "Compute G window table");
    let g_window = FixedBaseMSM::get_mul_window_size(
        // Verifier query
        assembly.num_inputs
        // A query
        + non_zero_a
        // C query 1
        + (sap_num_variables - (assembly.num_inputs - 1))
        // C query 2
        + sap_num_variables + 1
        // G gamma2 Z t
        + m_raw + 1,
    );
    let g_table = FixedBaseMSM::get_window_table::<E::G1Projective>(scalar_bits, g_window, g);
    timer_end!(g_window_time);

    // Compute H_gamma window table
    let h_gamma_time = timer_start!(|| "Compute H table");
    let h_gamma = h.into_affine().mul(gamma.into_repr());
    let h_gamma_window = FixedBaseMSM::get_mul_window_size(non_zero_a);
    let h_gamma_table =
        FixedBaseMSM::get_window_table::<E::G2Projective>(scalar_bits, h_gamma_window, h_gamma);
    timer_end!(h_gamma_time);

    // Generate the R1CS proving key
    let proving_key_time = timer_start!(|| "Generate the R1CS proving key");

    // Compute the A-query
    let a_time = timer_start!(|| "Calculate A");
    let mut a_query = FixedBaseMSM::batch_mul::<E, E::G1Projective>(
        scalar_bits,
        g_window,
        &g_table,
        &a.par_iter().map(|a| *a * &gamma).collect(),
    );
    timer_end!(a_time);

    // Compute the B-query
    let b_time = timer_start!(|| "Calculate B");
    let mut b_query = FixedBaseMSM::batch_mul::<E, E::G2Projective>(
        scalar_bits,
        h_gamma_window,
        &h_gamma_table,
        &a,
    );
    timer_end!(b_time);

    // Compute the G_gamma-query
    let g_gamma_time = timer_start!(|| "Calculate G gamma");
    let gamma_z = zt * &gamma;
    let alpha_beta = alpha + &beta;
    let ab_gamma_z = alpha_beta * &gamma * &zt;
    let g_gamma = g.into_affine().mul(gamma.into_repr());
    let g_gamma_z = g.into_affine().mul(gamma_z.into_repr());
    let h_gamma_z = h_gamma.into_affine().mul(zt.into_repr());
    let g_ab_gamma_z = g.into_affine().mul(ab_gamma_z.into_repr());
    let g_gamma2_z2 = g.into_affine().mul(gamma_z.square().into_repr());

    // Compute the vector G_gamma2_z_t := Z(t) * t^i * gamma^2 * G
    let gamma2_z_t = gamma_z * &gamma;
    let mut g_gamma2_z_t = FixedBaseMSM::batch_mul::<E, E::G1Projective>(
        scalar_bits,
        g_window,
        &g_table,
        &(0..m_raw + 1)
            .into_par_iter()
            .map(|i| gamma2_z_t * &(t.pow([i as u64])))
            .collect(),
    );
    timer_end!(g_gamma_time);

    // Compute the C_1-query
    let c1_time = timer_start!(|| "Calculate C1");
    let result = FixedBaseMSM::batch_mul::<E, E::G1Projective>(
        scalar_bits,
        g_window,
        &g_table,
        &(0..sap_num_variables + 1)
            .into_par_iter()
            .map(|i| c[i] * &gamma + &(a[i] * &alpha_beta))
            .collect(),
    );
    let (verifier_query, c_query_1) = result.split_at(assembly.num_inputs);
    timer_end!(c1_time);

    // Compute the C_2-query
    let c2_time = timer_start!(|| "Calculate C2");
    let double_gamma2_z = (zt * &gamma.square()).double();
    let mut c_query_2 = FixedBaseMSM::batch_mul::<E, E::G1Projective>(
        scalar_bits,
        g_window,
        &g_table,
        &(0..sap_num_variables + 1)
            .into_par_iter()
            .map(|i| a[i] * &double_gamma2_z)
            .collect(),
    );
    timer_end!(c2_time);

    timer_end!(proving_key_time);

    // Generate R1CS verification key
    let verifying_key_time = timer_start!(|| "Generate the R1CS verification key");
    let g_alpha = g.into_affine().mul(alpha.into_repr());
    let h_beta = h.into_affine().mul(beta.into_repr());
    timer_end!(verifying_key_time);

    let vk = VerifyingKey::<E> {
        h_g2:       h.into_affine(),
        g_alpha_g1: g_alpha.into_affine(),
        h_beta_g2:  h_beta.into_affine(),
        g_gamma_g1: g_gamma.into_affine(),
        h_gamma_g2: h_gamma.into_affine(),
        query:      verifier_query
            .into_par_iter()
            .map(|e| e.into_affine())
            .collect(),
    };

    let mut c_query_1 = c_query_1.to_vec();

    let batch_normalization_time = timer_start!(|| "Convert proving key elements to affine");
    E::G1Projective::batch_normalization(a_query.as_mut_slice());
    E::G2Projective::batch_normalization(b_query.as_mut_slice());
    E::G1Projective::batch_normalization(c_query_1.as_mut_slice());
    E::G1Projective::batch_normalization(c_query_2.as_mut_slice());
    E::G1Projective::batch_normalization(g_gamma2_z_t.as_mut_slice());
    timer_end!(batch_normalization_time);

    Ok(Parameters {
        vk,
        a_query: a_query.into_iter().map(Into::into).collect(),
        b_query: b_query.into_iter().map(Into::into).collect(),
        c_query_1: c_query_1.into_iter().map(Into::into).collect(),
        c_query_2: c_query_2.into_iter().map(Into::into).collect(),
        g_gamma_z: g_gamma_z.into_affine(),
        h_gamma_z: h_gamma_z.into_affine(),
        g_ab_gamma_z: g_ab_gamma_z.into_affine(),
        g_gamma2_z2: g_gamma2_z2.into_affine(),
        g_gamma2_z_t: g_gamma2_z_t.into_iter().map(Into::into).collect(),
    })
}
