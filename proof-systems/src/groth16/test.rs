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

mod bls12_377 {
    use super::*;
    use crate::groth16::{
        Parameters, Proof, VerifyingKey, PreparedVerifyingKey,
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use algebra::{curves::bls12_377::Bls12_377, fields::bls12_377::Fr, UniformRand,
            ToBytes, FromBytes, to_bytes,
    };
    use rand::thread_rng;
    use std::ops::MulAssign;

    #[test]
    fn prove_and_verify() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<Bls12_377, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<Bls12_377>(&params.vk);

        for _ in 0..100 {
            let a = Fr::rand(rng);
            let b = Fr::rand(rng);
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

    #[test]
    fn serialize_deserialize() {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<Bls12_377, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<Bls12_377>::read(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<Bls12_377>::read(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);


        let a = Fr::rand(rng);
        let b = Fr::rand(rng);
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

        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = Proof::<Bls12_377>::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_serialized = to_bytes!(pvk).unwrap();
        let pvk_deserialized = PreparedVerifyingKey::<Bls12_377>::read(pvk_serialized.as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap());
    }
}

mod sw6 {
    use super::*;
    use crate::groth16::{
        Parameters, Proof, VerifyingKey, PreparedVerifyingKey,
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use rand::thread_rng;

    use algebra::{curves::sw6::SW6, fields::sw6::Fr as SW6Fr, Field, UniformRand,
                  ToBytes, FromBytes, to_bytes,
    };

    #[test]
    fn prove_and_verify() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<SW6, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<SW6>(&params.vk);

        let a = SW6Fr::rand(rng);
        let b = SW6Fr::rand(rng);
        let c = a * &b;

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
        assert!(!verify_proof(&pvk, &proof, &[SW6Fr::zero()]).unwrap());
    }

    #[test]
    fn serialize_deserialize() {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<SW6, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<SW6>::read(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<SW6>::read(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);


        let a = SW6Fr::rand(rng);
        let b = SW6Fr::rand(rng);
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

        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = Proof::<SW6>::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_serialized = to_bytes!(pvk).unwrap();
        let pvk_deserialized = PreparedVerifyingKey::<SW6>::read(pvk_serialized.as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap());
    }
}

mod mnt4753 {
    use super::*;
    use crate::groth16::{
        Parameters, Proof, VerifyingKey, PreparedVerifyingKey,
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use rand::thread_rng;

    use algebra::{curves::mnt4753::MNT4, fields::mnt4753::Fr as MNT4Fr, Field, UniformRand,
                  ToBytes, FromBytes, to_bytes};

    #[test]
    fn prove_and_verify() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<MNT4, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<MNT4>(&params.vk);

        let a = MNT4Fr::rand(rng);
        let b = MNT4Fr::rand(rng);
        let c = a * &b;

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
        assert!(!verify_proof(&pvk, &proof, &[MNT4Fr::zero()]).unwrap());
    }

    #[test]
    fn serialize_deserialize() {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<MNT4, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<MNT4>::read(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<MNT4>::read(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);


        let a = MNT4Fr::rand(rng);
        let b = MNT4Fr::rand(rng);
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

        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = Proof::<MNT4>::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_serialized = to_bytes!(pvk).unwrap();
        let pvk_deserialized = PreparedVerifyingKey::<MNT4>::read(pvk_serialized.as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap());
    }
}

mod mnt6753 {
    use super::*;
    use crate::groth16::{
        Parameters, Proof, VerifyingKey, PreparedVerifyingKey,
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use rand::thread_rng;

    use algebra::{curves::mnt6753::MNT6, fields::mnt6753::Fr as MNT6Fr, Field, UniformRand,
                  ToBytes, FromBytes, to_bytes,};

    #[test]
    fn prove_and_verify() {
        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<MNT6, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<MNT6>(&params.vk);

        let a = MNT6Fr::rand(rng);
        let b = MNT6Fr::rand(rng);
        let c = a * &b;

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
        assert!(!verify_proof(&pvk, &proof, &[MNT6Fr::zero()]).unwrap());
    }

    #[test]
    fn serialize_deserialize() {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<MNT6, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<MNT6>::read(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<MNT6>::read(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);


        let a = MNT6Fr::rand(rng);
        let b = MNT6Fr::rand(rng);
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

        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = Proof::<MNT6>::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_serialized = to_bytes!(pvk).unwrap();
        let pvk_deserialized = PreparedVerifyingKey::<MNT6>::read(pvk_serialized.as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap());
    }
}
