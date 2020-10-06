use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use algebra_core::{
    msm::VariableBaseMSM, AffineCurve, PairingEngine, PrimeField, ProjectiveCurve, UniformRand,
};

use crate::{r1cs_to_sap::R1CStoSAP, Parameters, Proof, Vec};
use ff_fft::{cfg_into_iter, EvaluationDomain};

use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

use core::ops::{AddAssign, MulAssign};

pub fn create_random_proof<E, C, D, R>(
    circuit: C,
    params: &Parameters<E>,
    rng: &mut R,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    D: EvaluationDomain<E::Fr>,
    R: Rng,
{
    let d1 = E::Fr::rand(rng);
    let d2 = E::Fr::rand(rng);
    let r = E::Fr::rand(rng);

    create_proof::<E, C, D>(circuit, params, d1, d2, r)
}

pub fn create_proof<E, C, D>(
    circuit: C,
    params: &Parameters<E>,
    d1: E::Fr,
    d2: E::Fr,
    r: E::Fr,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    D: EvaluationDomain<E::Fr>,
{
    let prover_time = start_timer!(|| "GrothMaller17::Prover");
    let cs = ConstraintSystem::new_ref();

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(cs.clone())?;
    end_timer!(synthesis_time);

    let lc_time = start_timer!(|| "Inlining LCs");
    cs.inline_all_lcs();
    end_timer!(lc_time);

    let witness_map_time = start_timer!(|| "R1CS to SAP witness map");
    let (full_input_assignment, h, _) = R1CStoSAP::witness_map::<E, D>(cs.clone(), &d1, &d2)?;
    end_timer!(witness_map_time);
    let prover = cs.borrow().unwrap();
    let num_inputs = prover.num_instance_variables;

    let input_assignment = full_input_assignment[1..num_inputs]
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    let aux_assignment = cfg_into_iter!(full_input_assignment[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    drop(full_input_assignment);

    let h_input = h[0..num_inputs]
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    let h_aux = cfg_into_iter!(h[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    drop(h);

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");
    let (a_inputs_source, a_aux_source) = params.get_a_query(num_inputs)?;
    let a_inputs_acc = VariableBaseMSM::multi_scalar_mul(a_inputs_source, &input_assignment);
    let a_aux_acc = VariableBaseMSM::multi_scalar_mul(a_aux_source, &aux_assignment);

    let r_g = params.get_g_gamma_z()?.mul(r);
    let d1_g = params.get_g_gamma_z()?.mul(d1);

    let mut g_a = r_g;
    g_a.add_assign(&params.get_a_query_full()?[0].into_projective());
    g_a.add_assign(&d1_g);
    g_a.add_assign(&a_inputs_acc);
    g_a.add_assign(&a_aux_acc);
    end_timer!(a_acc_time);

    // Compute B
    let b_acc_time = start_timer!(|| "Compute B");

    let (b_inputs_source, b_aux_source) = params.get_b_query(num_inputs)?;
    let b_inputs_acc = VariableBaseMSM::multi_scalar_mul(b_inputs_source, &input_assignment);
    let b_aux_acc = VariableBaseMSM::multi_scalar_mul(b_aux_source, &aux_assignment);

    let r_h = params.get_h_gamma_z()?.mul(r);
    let d1_h = params.get_h_gamma_z()?.mul(d1);

    let mut g_b = r_h;
    g_b.add_assign(&params.get_b_query_full()?[0].into_projective());
    g_b.add_assign(&d1_h);
    g_b.add_assign(&b_inputs_acc);
    g_b.add_assign(&b_aux_acc);
    end_timer!(b_acc_time);

    // Compute C
    let c_acc_time = start_timer!(|| "Compute C");
    let r_2 = r + &r;
    let r2 = r * &r;
    let d1_r_2 = d1 * &r_2;

    let c1_acc_time = start_timer!(|| "Compute C1");
    let (_, c1_aux_source) = params.get_c_query_1(0)?;
    let c1_acc = VariableBaseMSM::multi_scalar_mul(c1_aux_source, &aux_assignment);
    end_timer!(c1_acc_time);

    let c2_acc_time = start_timer!(|| "Compute C2");

    let (c2_inputs_source, c2_aux_source) = params.get_c_query_2(num_inputs)?;
    let c2_inputs_acc = VariableBaseMSM::multi_scalar_mul(c2_inputs_source, &input_assignment);
    let c2_aux_acc = VariableBaseMSM::multi_scalar_mul(c2_aux_source, &aux_assignment);

    let c2_acc = c2_inputs_acc + &c2_aux_acc;
    end_timer!(c2_acc_time);

    // Compute G
    let g_acc_time = start_timer!(|| "Compute G");

    let (g_inputs_source, g_aux_source) = params.get_g_gamma2_z_t(num_inputs)?;
    let g_inputs_acc = VariableBaseMSM::multi_scalar_mul(g_inputs_source, &h_input);
    let g_aux_acc = VariableBaseMSM::multi_scalar_mul(g_aux_source, &h_aux);

    let g_acc = g_inputs_acc + &g_aux_acc;
    end_timer!(g_acc_time);

    let r2_g_gamma2_z2 = params.get_g_gamma2_z2()?.mul(r2);
    let r_g_ab_gamma_z = params.get_g_ab_gamma_z()?.mul(r);
    let d1_g_ab_gamma_z = params.get_g_ab_gamma_z()?.mul(d1);
    let r_c0 = params.get_c_query_2_full()?[0].mul(r);
    let r2_d1_g_gamma2_z2 = params.get_g_gamma2_z2()?.mul(d1_r_2);
    let d2_g_gamma2_z_t0 = params.get_g_gamma2_z_t_full()?[0].mul(d2);
    let mut r_c2_exp = c2_acc;
    r_c2_exp.mul_assign(r);

    let mut g_c = c1_acc;
    g_c.add_assign(&r2_g_gamma2_z2);
    g_c.add_assign(&r_g_ab_gamma_z);
    g_c.add_assign(&d1_g_ab_gamma_z);
    g_c.add_assign(&r_c0);
    g_c.add_assign(&r2_d1_g_gamma2_z2);
    g_c.add_assign(&r_c2_exp);
    g_c.add_assign(&d2_g_gamma2_z_t0);
    g_c.add_assign(&g_acc);
    end_timer!(c_acc_time);

    end_timer!(prover_time);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_c.into_affine(),
    })
}
