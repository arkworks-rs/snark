use algebra_core::{AffineCurve, One, PairingEngine, ProjectiveCurve};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use crate::SynthesisError;

use core::ops::Neg;

pub fn prepare_verifying_key<E: PairingEngine>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        g_alpha: vk.g_alpha_g1,
        h_beta: vk.h_beta_g2,
        g_alpha_h_beta_ml: E::miller_loop([(vk.g_alpha_g1.into(), vk.h_beta_g2.into())].iter()),
        g_gamma_pc: vk.g_gamma_g1.into(),
        h_gamma_pc: vk.h_gamma_g2.into(),
        h_pc: vk.h_g2.into(),
        query: vk.query.clone(),
    }
}

pub fn verify_proof<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::Fr],
) -> Result<bool, SynthesisError> {
    if (public_inputs.len() + 1) != pvk.query.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
    // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]

    let mut g_psi = pvk.query[0].into_projective();
    for (i, b) in public_inputs.iter().zip(pvk.query.iter().skip(1)) {
        g_psi += &b.mul(*i);
    }

    let mut test1_a_g_alpha = proof.a.into_projective();
    test1_a_g_alpha.add_assign_mixed(&pvk.g_alpha);
    let test1_a_g_alpha = test1_a_g_alpha.into_affine();

    let mut test1_b_h_beta = proof.b.into_projective();
    test1_b_h_beta.add_assign_mixed(&pvk.h_beta);
    let test1_b_h_beta = test1_b_h_beta.into_affine();

    let test1_r1 = pvk.g_alpha_h_beta_ml;

    let test1_r2_input: &[(E::G1Prepared, E::G2Prepared)] = &[
        (test1_a_g_alpha.neg().into(), test1_b_h_beta.into()),
        (g_psi.into().into(), pvk.h_gamma_pc.clone()),
        (proof.c.into(), pvk.h_pc.clone()),
    ];

    let test1_r2 = E::miller_loop(test1_r2_input.iter());
    let test1 = E::final_exponentiation(&(test1_r2 * &test1_r1)).unwrap();

    // e(A, H^{gamma}) = e(G^{gamma}, B)
    let test2_exp_input: &[(E::G1Prepared, E::G2Prepared)] = &[
        (proof.a.into(), pvk.h_gamma_pc.clone()),
        (pvk.g_gamma_pc.clone(), proof.b.neg().into()),
    ];

    let test2_exp = E::miller_loop(test2_exp_input.iter());
    let test2 = E::final_exponentiation(&test2_exp).unwrap();

    Ok(test1.is_one() && test2.is_one())
}
