use algebra::{Field, FromBytesChecked};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
struct MySillyCircuit<F: Field> {
    a: Option<F>,
    b: Option<F>,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for MySillyCircuit<ConstraintF> {
    fn generate_constraints<CS: ConstraintSystem<ConstraintF>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let a = cs.alloc(|| "a", || self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.alloc(|| "b", || self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.alloc_input(
            || "c",
            || {
                let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                a.mul_assign(&b);
                Ok(a)
            },
        )?;

        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);
        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);
        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);
        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);
        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);
        cs.enforce(|| "a*b=c", |lc| lc + a, |lc| lc + b, |lc| lc + c);

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::groth16::{Parameters, Proof, VerifyingKey, PreparedVerifyingKey, create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof, create_proof_no_zk};

    use algebra::{UniformRand, ToBytes, FromBytes, to_bytes, PairingEngine};
    use rand::thread_rng;
    use std::ops::MulAssign;

    fn prove_and_verify<E: PairingEngine>(zk: bool) {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<E, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<E>(&params.vk);

        for _ in 0..100 {
            let a = E::Fr::rand(rng);
            let b = E::Fr::rand(rng);
            let mut c = a;
            c.mul_assign(&b);

            let proof = if zk {
                create_random_proof(
                    MySillyCircuit {
                        a: Some(a),
                        b: Some(b),
                    },
                    &params,
                    rng,
                ).unwrap()
            } else {
                create_proof_no_zk(
                    MySillyCircuit {
                        a: Some(a),
                        b: Some(b),
                    },
                    &params,
                ).unwrap()
            };

            assert!(verify_proof(&pvk, &proof, &[c]).unwrap());
            assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
        }
    }

    fn serialize_deserialize<E: PairingEngine>() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<E, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<E>::read_checked(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<E>::read_checked(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);

        let a = E::Fr::rand(rng);
        let b = E::Fr::rand(rng);
        let c = a * &b;

        let proof = create_random_proof(
            MySillyCircuit {
                a: Some(a),
                b: Some(b),
            },
            &params_deserialized,
            rng,
        )
            .unwrap();

        let proof_deserialized = Proof::<E>::read_checked(to_bytes!(proof).unwrap().as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);
        drop(proof);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_deserialized = PreparedVerifyingKey::<E>::read(to_bytes!(pvk).unwrap().as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap())
    }

    #[test]
    fn bls12_377_groth16_test() {
        prove_and_verify::<algebra::curves::bls12_377::Bls12_377>(true);
        prove_and_verify::<algebra::curves::bls12_377::Bls12_377>(false);
        serialize_deserialize::<algebra::curves::bls12_377::Bls12_377>();
    }

    #[test]
    fn sw6_groth16_test() {
        prove_and_verify::<algebra::curves::sw6::SW6>(true);
        prove_and_verify::<algebra::curves::sw6::SW6>(false);
        serialize_deserialize::<algebra::curves::sw6::SW6>();
    }

    #[test]
    fn mnt4753_groth16_test() {
        prove_and_verify::<algebra::curves::mnt4753::MNT4>(true);
        prove_and_verify::<algebra::curves::mnt4753::MNT4>(false);
        serialize_deserialize::<algebra::curves::mnt4753::MNT4>();
    }

    #[test]
    fn mnt6753_groth16_test() {
        prove_and_verify::<algebra::curves::mnt6753::MNT6>(true);
        prove_and_verify::<algebra::curves::mnt6753::MNT6>(false);
        serialize_deserialize::<algebra::curves::mnt6753::MNT6>();
    }

    #[test]
    fn bn_382_groth16_test() {
        prove_and_verify::<algebra::curves::bn_382::Bn382>(true);
        prove_and_verify::<algebra::curves::bn_382::Bn382>(false);
        serialize_deserialize::<algebra::curves::bn_382::Bn382>();
    }
}