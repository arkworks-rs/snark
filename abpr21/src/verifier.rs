use ark_ec::msm::FixedBaseMSM;
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField, to_bytes,};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};

use core::ops::{AddAssign, Neg};

use blake2::{Blake2b, Digest};

/// Prepare the verifying key `vk` for use in proof verification.
pub fn prepare_verifying_key<E: PairingEngine>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        gamma_g2_neg_pc: vk.gamma_g2.neg().into(),
    }
}

/// Prepare proof inputs for use with [`verify_proof_with_prepared_inputs`], wrt the prepared
/// verification key `pvk` and instance public inputs.
pub fn prepare_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    public_inputs: &[E::Fr],
) -> R1CSResult<E::G1Projective> {
    if (public_inputs.len() + 1) != pvk.vk.gamma_abc_g1.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    let mut g_ic = pvk.vk.gamma_abc_g1[0].into_projective();
    for (i, b) in public_inputs.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
        g_ic.add_assign(&b.mul(i.into_repr()));
    }

    Ok(g_ic)
}


/// Verify a proof `proof` against the prepared verification key `pvk` and prepared public
/// inputs. This should be preferred over [`verify_proof`] if the instance's public inputs are
/// known in advance.
pub fn verify_proof_with_prepared_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    prepared_inputs: &E::G1Projective,
) -> R1CSResult<bool> {

    let hash = Blake2b::new()
    .chain(to_bytes!(&proof.a).unwrap())
    .chain(to_bytes!(&proof.b).unwrap())
    .chain(to_bytes!(&proof.delta_prime).unwrap());
    let mut output = [0u8; 64];
    output.copy_from_slice(&hash.finalize());

    let m_fr = E::Fr::from_le_bytes_mod_order(&output);
    //println!("m_fr verifier {0}", m_fr);
    let mut delta_prime_delta_m = pvk.vk.delta_g2.mul(m_fr);
    delta_prime_delta_m.add_assign_mixed(&proof.delta_prime);


    let qap = E::miller_loop(
        [
            (proof.a.into(), proof.b.into()),
            (
                prepared_inputs.into_affine().into(),
                pvk.gamma_g2_neg_pc.clone(),
            ),
            (proof.c.into(), delta_prime_delta_m.into_affine().neg().into()),
        ]
        .iter(),
    );

    let test = E::final_exponentiation(&qap).ok_or(SynthesisError::UnexpectedIdentity)?;

    Ok(test == pvk.vk.alpha_g1_beta_g2)
}

/// Verify a proof `proof` against the prepared verification key `pvk`,
/// with respect to the instance `public_inputs`.
pub fn verify_proof<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::Fr],
) -> R1CSResult<bool> {
    let prepared_inputs = prepare_inputs(pvk, public_inputs)?;
    verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
}





/// Verify a vector of proofs `proofs` against the prepared verification key `pvk` and prepared public
/// inputs vector.
pub fn vec_verify_proof_with_prepared_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proofs: &Vec<Proof<E>>,
    prepared_inputs: &Vec<E::G1Projective>,
) -> R1CSResult<bool> {


    let num_proofs = proofs.len();
    let mut m_fr: Vec<E::Fr> = Vec::with_capacity(num_proofs as usize);
    
    let start = ark_std::time::Instant::now();
    for proof in proofs.iter() {
        let hash = Blake2b::new()
        .chain(to_bytes!(&proof.a).unwrap())
        .chain(to_bytes!(&proof.b).unwrap())
        .chain(to_bytes!(&proof.delta_prime).unwrap());
        let mut output = [0u8; 64];
        output.copy_from_slice(&hash.finalize());
        m_fr.push(E::Fr::from_le_bytes_mod_order(&output));
        
    }

    let scalar_bits = E::Fr::size_in_bits();

    let delta_g2_window = FixedBaseMSM::get_mul_window_size(num_proofs);
    let delta_g2_table =
        FixedBaseMSM::get_window_table::<E::G2Projective>(scalar_bits, delta_g2_window, pvk.vk.delta_g2.into_projective());
    let elem_g2 =
        FixedBaseMSM::multi_scalar_mul::<E::G2Projective>(scalar_bits, delta_g2_window, &delta_g2_table, &m_fr);

    println!(
        "Hashing + Exponentiation (G2) time is {}ns per proof doing {} exponentiations",
        start.elapsed().as_nanos() / num_proofs as u128,
        num_proofs
    );
   


    
    let mut bool_results: Vec<_> = Vec::new();
    for ((x,y),z) in elem_g2.iter().zip(proofs.iter()).zip(prepared_inputs.iter()){
       // x -> [m_fr * delta]_2    ;;  y -> proof  ;;  z -> prepared_inputs 

        let tmp1 = E::final_exponentiation(&E::miller_loop(
            [
                (y.a.into(), y.b.into()),
                (
                    z.into_affine().into(),
                    pvk.gamma_g2_neg_pc.clone(),
                ),
                (y.c.into(), (*x + y.delta_prime.into_projective()).neg().into_affine().into()),
            ]
            .iter(),
        )).unwrap();
        let tmp2 = pvk.vk.alpha_g1_beta_g2;
        let tmp =  tmp1 == tmp2 ;

        bool_results.push(tmp);

    }
    
    let result = bool_results.iter().fold(true, |total, next| {total && *next});
    //println!("result is {:?}", result);
    

    
    Ok(result)
    
    

}

/// Verify a vector of proofs `proofs` against the prepared verification key `pvk`,
/// with respect to the instances `public_inputs`'s.
pub fn vec_verify_proof<E: PairingEngine>(
    vk: &VerifyingKey<E>,
    proofs: &Vec<Proof<E>>,
    public_inputs: &Vec<Vec<E::Fr>>,
) -> R1CSResult<bool> {
    //let pvk = prepare_verifying_key(vk);
    let mut prepared_inputs: Vec<_> = Vec::new();
    for (_,pub_input) in public_inputs.iter().enumerate(){
        let pvk = prepare_verifying_key(vk);
        prepared_inputs.push(prepare_inputs(&pvk, pub_input)?);
    }
    let pvk = prepare_verifying_key(vk);
    vec_verify_proof_with_prepared_inputs(&pvk, proofs, &prepared_inputs)
}
