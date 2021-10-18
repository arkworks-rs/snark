use crate::{r1cs_to_qap::R1CStoQAP, Proof, ProvingKey};
use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField, UniformRand, Zero, One, to_bytes};
use ark_poly::GeneralEvaluationDomain;
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, OptimizationGoal, Result as R1CSResult,
};
use ark_std::rand::Rng;
use ark_std::{cfg_into_iter, cfg_iter, vec::Vec};

use blake2::{Blake2b, Digest};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Create a proof that is zero-knowledge.
/// This method samples randomness for zero knowledges via `rng`.
#[inline]
pub fn create_random_proof<E, C, R>(
    circuit: C,
    pk: &ProvingKey<E>,
    rng: &mut R,
) -> R1CSResult<Proof<E>>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    let r = E::Fr::rand(rng);
    let s = E::Fr::rand(rng);
    let mut zeta = E::Fr::zero();
    while zeta.is_zero(){
        zeta = E::Fr::rand(rng);
    }

    create_proof::<E, C>(circuit, pk, r, s, zeta)
}

/// Create a proof that is *not* zero-knowledge.
#[inline]
pub fn create_proof_no_zk<E, C>(circuit: C, pk: &ProvingKey<E>) -> R1CSResult<Proof<E>>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    create_proof::<E, C>(circuit, pk, E::Fr::zero(), E::Fr::zero(), E::Fr::one())
}

/// Create a proof using randomness `r` and `s`.
#[inline]
pub fn create_proof<E, C>(
    circuit: C,
    pk: &ProvingKey<E>,
    r: E::Fr,
    s: E::Fr,
    zeta: E::Fr,
) -> R1CSResult<Proof<E>>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    type D<F> = GeneralEvaluationDomain<F>;

    let prover_time = start_timer!(|| "BPR20::Prover");
    let cs = ConstraintSystem::new_ref();

    // Set the optimization goal
    cs.set_optimization_goal(OptimizationGoal::Constraints);

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(cs.clone())?;
    debug_assert!(cs.is_satisfied().unwrap());
    end_timer!(synthesis_time);

    let lc_time = start_timer!(|| "Inlining LCs");
    cs.finalize();
    end_timer!(lc_time);

    let witness_map_time = start_timer!(|| "R1CS to QAP witness map");
    let h = R1CStoQAP::witness_map::<E::Fr, D<E::Fr>>(cs.clone())?;
    end_timer!(witness_map_time);

    let c_acc_time = start_timer!(|| "Compute C");

    // Compute C
    let prover = cs.borrow().unwrap();

    let aux_assignment = cfg_iter!(prover.witness_assignment)
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    
 
    
    let delta_prime_g1 = pk.delta_g1.clone().mul(zeta).into_affine();
    let delta_prime_g2 = pk.vk.delta_g2.clone().mul(zeta).into_affine();


    let r_s_delta_g1 = delta_prime_g1.mul((&(r*s)).into_repr());
    
    end_timer!(c_acc_time);

    let input_assignment = prover.instance_assignment[1..]
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    
    

    let assignment = [&input_assignment[..], &aux_assignment[..]].concat();
    drop(aux_assignment);

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");
    let r_g1 = delta_prime_g1.mul(r);

    let g_a = calculate_coeff(r_g1, &pk.a_query, pk.vk.alpha_g1, &assignment);

    let s_g_a = g_a.mul(&s.into_repr());
    end_timer!(a_acc_time);

    // Compute B in G1 if needed
    let g1_b = if !r.is_zero() {
        let b_g1_acc_time = start_timer!(|| "Compute B in G1");
        let s_g1 = delta_prime_g1.mul(s);
        let g1_b = calculate_coeff(s_g1, &pk.b_g1_query, pk.beta_g1, &assignment);

        end_timer!(b_g1_acc_time);

        g1_b
    } else {
        E::G1Projective::zero()
    };

    // Compute B in G2
    let b_g2_acc_time = start_timer!(|| "Compute B in G2");
    let s_g2 = delta_prime_g2.mul(s);
    let g2_b = calculate_coeff(s_g2, &pk.b_g2_query, pk.vk.beta_g2, &assignment);
    let r_g1_b = g1_b.mul(&r.into_repr());
    drop(assignment);

    end_timer!(b_g2_acc_time);

    let c_time = start_timer!(|| "Finish C");

    //Compute the hash message
    let hash = Blake2b::new()
    .chain(to_bytes!(&g_a.into_affine()).unwrap())
    .chain(to_bytes!(&g2_b.into_affine()).unwrap())
    .chain(to_bytes!(&delta_prime_g2).unwrap());
    let mut output = [0u8; 64];
    output.copy_from_slice(&hash.finalize());
    
    let m_fr = E::Fr::from_le_bytes_mod_order(&output);
    //println!("m_fr prover {0}", m_fr);
    let factor = zeta * (zeta + m_fr).inverse().unwrap();   
    let zeta_m_inv = (zeta + m_fr).inverse().unwrap();

    let h_assignment = cfg_into_iter!(h).map(|s| (s*zeta_m_inv).into()).collect::<Vec<_>>();
    let h_acc = VariableBaseMSM::multi_scalar_mul(&pk.h_query, &h_assignment);
    let aux_assignment_unscaled = cfg_iter!(prover.witness_assignment)
        .map(|s| (*s*zeta_m_inv).into()) 
        .collect::<Vec<_>>();

    let l_aux_acc = VariableBaseMSM::multi_scalar_mul(&pk.l_query, &aux_assignment_unscaled);
    
    let mut g_c = s_g_a.mul(&factor.into_repr());
    g_c += &r_g1_b.mul(&factor.into_repr());
    g_c -= &r_s_delta_g1.mul(&factor.into_repr());
    g_c += &l_aux_acc;
    g_c += &h_acc;
    end_timer!(c_time);

    end_timer!(prover_time);
    
    
    Ok(Proof {
        a: g_a.into_affine(),
        b: g2_b.into_affine(),
        c: g_c.into_affine(),
        delta_prime: delta_prime_g2,
    })
}


fn calculate_coeff<G: AffineCurve>(
    initial: G::Projective,
    query: &[G],
    vk_param: G,
    assignment: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    let el = query[0];
    let acc = VariableBaseMSM::multi_scalar_mul(&query[1..], assignment);

    let mut res = initial;
    res.add_assign_mixed(&el);
    res += &acc;
    res.add_assign_mixed(&vk_param);

    res
}
