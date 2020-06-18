use algebra::Field;
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

mod test{
    use super::*;
    use crate::groth16::{
        Parameters, Proof, VerifyingKey, PreparedVerifyingKey,
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use algebra::{
        curves::{bls12_377::Bls12_377, sw6::SW6, mnt4753::MNT4, mnt6753::MNT6},
        UniformRand, ToBytes, FromBytes, to_bytes, PairingEngine
    };
    use rand::thread_rng;
    use std::ops::MulAssign;

    fn test_prove_and_verify<E: PairingEngine>() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<E, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<E>(&params.vk);

        for _ in 0..100 {
            let a = <E as PairingEngine>::Fr::rand(rng);
            let b = <E as PairingEngine>::Fr::rand(rng);
            let mut c = a;
            c.mul_assign(&b);

            let proof = create_random_proof(
                MySillyCircuit {
                    a: Some(a),
                    b: Some(b),
                },
                &params,
                rng,
            )
                .unwrap();

            assert!(verify_proof(&pvk, &proof, &[c]).unwrap());
            assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
        }
    }

    fn test_serialize_deserialize<E: PairingEngine>() {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<E, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_deserialized = Parameters::<E>::read(to_bytes!(params).unwrap().as_slice()).unwrap();
        assert_eq!(params, params_deserialized);
        drop(params);

        let vk_deserialized = VerifyingKey::<E>::read(to_bytes!(vk).unwrap().as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);
        drop(vk);

        let a = <E as PairingEngine>::Fr::rand(rng);
        let b = <E as PairingEngine>::Fr::rand(rng);
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

        let proof_deserialized = Proof::<E>::read(to_bytes!(proof).unwrap().as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);
        drop(proof);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_deserialized = PreparedVerifyingKey::<E>::read(to_bytes!(pvk).unwrap().as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap())
    }

    #[test]
    fn prove_verify() {
        test_prove_and_verify::<Bls12_377>();
        test_prove_and_verify::<SW6>();
        test_prove_and_verify::<MNT4>();
        test_prove_and_verify::<MNT6>();
    }

    #[test]
    fn serialize_deserialize() {
        test_serialize_deserialize::<Bls12_377>();
        test_serialize_deserialize::<SW6>();
        test_serialize_deserialize::<MNT4>();
        test_serialize_deserialize::<MNT6>();
    }
}