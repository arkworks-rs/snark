use ark_ec::msm::FixedBaseMSM;
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{to_bytes, PrimeField};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};

use ark_std::{ops::{AddAssign, Neg}, vec::Vec};

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
            (
                proof.c.into(),
                delta_prime_delta_m.into_affine().neg().into(),
            ),
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

/// Verify a list of proofs `proofs` against the prepared verification key `pvk` and prepared public
/// inputs vector.
pub fn verify_proofs_with_prepared_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proofs: &[Proof<E>],
    prepared_inputs: &[E::G1Projective],
) -> R1CSResult<bool> {
    let verifier_time = start_timer!(|| "ABPR21::BatchVerifier w/ prepared inputs");

    let num_proofs = proofs.len();
    let m_fr: Vec<_> = proofs.iter().map(|proof| {
        let hash = Blake2b::new()
            .chain(to_bytes!(&proof.a).unwrap())
            .chain(to_bytes!(&proof.b).unwrap())
            .chain(to_bytes!(&proof.delta_prime).unwrap());
        let mut output = [0u8; 64];
        output.copy_from_slice(&hash.finalize());
        E::Fr::from_le_bytes_mod_order(&output)
    }).collect();

    let scalar_bits = E::Fr::size_in_bits();

    let delta_g2_window = FixedBaseMSM::get_mul_window_size(num_proofs);
    let delta_g2_table = FixedBaseMSM::get_window_table::<E::G2Projective>(
        scalar_bits,
        delta_g2_window,
        pvk.vk.delta_g2.into_projective(),
    );
    let elem_g2 = FixedBaseMSM::multi_scalar_mul::<E::G2Projective>(
        scalar_bits,
        delta_g2_window,
        &delta_g2_table,
        &m_fr,
    );

    let result = cfg_iter!(elem_g2)
        .zip(proofs.iter())
        .zip(prepared_inputs.iter())
        .map(|((x, y), z)| {
            // x -> [m_fr * delta]_2    ;;  y -> proof  ;;  z -> prepared_inputs
            let tmp1 = E::final_exponentiation(&E::miller_loop(
                [
                    (y.a.into(), y.b.into()),
                    (z.into_affine().into(), pvk.gamma_g2_neg_pc.clone()),
                    (
                        y.c.into(),
                        (*x + y.delta_prime.into_projective())
                            .neg()
                            .into_affine()
                            .into(),
                    ),
                ]
                .iter(),
            ))
            .unwrap();
            tmp1 == pvk.vk.alpha_g1_beta_g2
        }).all(ark_std::convert::identity);
    end_timer!(verifier_time);
    Ok(result)
}

/// Verify a vector of proofs `proofs` against the prepared verification key `pvk`,
/// with respect to the instances `public_inputs`'s.
pub fn verify_proofs<E: PairingEngine>(
    vk: &VerifyingKey<E>,
    proofs: &[Proof<E>],
    public_inputs: &[Vec<E::Fr>],
) -> R1CSResult<bool> {
    let pvk = prepare_verifying_key(vk);
    let prepared_inputs = cfg_iter!(public_inputs).map(|inp| prepare_inputs(&pvk, inp)).collect::<Result<Vec<_>, _>>()?;
    verify_proofs_with_prepared_inputs(&pvk, proofs, &prepared_inputs)
}
