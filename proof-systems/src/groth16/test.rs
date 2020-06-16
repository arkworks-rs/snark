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

    fn test_serialize_deserialize<E: PairingEngine>() -> (Vec<u8>, Vec<u8>, Vec<u8>) {

        let rng = &mut thread_rng();

        let params =
            generate_random_parameters::<E, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let vk = params.vk.clone();

        let params_serialized = to_bytes!(params).unwrap();
        let params_deserialized = Parameters::<E>::read(params_serialized.as_slice()).unwrap();
        assert_eq!(params, params_deserialized);

        let vk_serialized = to_bytes!(vk).unwrap();
        let vk_deserialized = VerifyingKey::<E>::read(vk_serialized.as_slice()).unwrap();
        assert_eq!(vk, vk_deserialized);


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

        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = Proof::<E>::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);

        let pvk = prepare_verifying_key(&vk_deserialized);
        let pvk_serialized = to_bytes!(pvk).unwrap();
        let pvk_deserialized = PreparedVerifyingKey::<E>::read(pvk_serialized.as_slice()).unwrap();
        assert_eq!(pvk, pvk_deserialized);

        assert!(verify_proof(&pvk_deserialized, &proof_deserialized, &[c]).unwrap());

        (params_serialized, proof_serialized, vk_serialized)
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

        // Basic checks
        let (params_bls, proof_bls, vk_bls) = test_serialize_deserialize::<Bls12_377>();
        let (params_sw6, proof_sw6, vk_sw6) = test_serialize_deserialize::<SW6>();
        let (params_mnt4, proof_mnt4, vk_mnt4) = test_serialize_deserialize::<MNT4>();
        let (params_mnt6, proof_mnt6, vk_mnt6) = test_serialize_deserialize::<MNT6>();

        // Cross checks: let's assert that deserialization fails when composing fields/points are not valid

        // BLS12-377
        assert!(Parameters::<Bls12_377>::read(params_sw6.as_slice()).is_err());
        assert!(Parameters::<Bls12_377>::read(params_mnt4.as_slice()).is_err());
        assert!(Parameters::<Bls12_377>::read(params_mnt6.as_slice()).is_err());

        assert!(Proof::<Bls12_377>::read(proof_sw6.as_slice()).is_err());
        assert!(Proof::<Bls12_377>::read(proof_mnt4.as_slice()).is_err());
        assert!(Proof::<Bls12_377>::read(proof_mnt6.as_slice()).is_err());

        assert!(VerifyingKey::<Bls12_377>::read(vk_sw6.as_slice()).is_err());
        assert!(VerifyingKey::<Bls12_377>::read(vk_mnt4.as_slice()).is_err());
        assert!(VerifyingKey::<Bls12_377>::read(vk_mnt6.as_slice()).is_err());

        // SW6
        assert!(Parameters::<SW6>::read(params_bls.as_slice()).is_err());
        assert!(Parameters::<SW6>::read(params_mnt4.as_slice()).is_err());
        assert!(Parameters::<SW6>::read(params_mnt6.as_slice()).is_err());

        assert!(Proof::<SW6>::read(proof_bls.as_slice()).is_err());
        assert!(Proof::<SW6>::read(proof_mnt4.as_slice()).is_err());
        assert!(Proof::<SW6>::read(proof_mnt6.as_slice()).is_err());

        assert!(VerifyingKey::<SW6>::read(vk_bls.as_slice()).is_err());
        assert!(VerifyingKey::<SW6>::read(vk_mnt4.as_slice()).is_err());
        assert!(VerifyingKey::<SW6>::read(vk_mnt6.as_slice()).is_err());

        //MNT4-753
        assert!(Parameters::<MNT4>::read(params_bls.as_slice()).is_err());
        assert!(Parameters::<MNT4>::read(params_sw6.as_slice()).is_err());
        assert!(Parameters::<MNT4>::read(params_mnt6.as_slice()).is_err());

        assert!(Proof::<MNT4>::read(proof_bls.as_slice()).is_err());
        assert!(Proof::<MNT4>::read(proof_sw6.as_slice()).is_err());
        assert!(Proof::<MNT4>::read(proof_mnt6.as_slice()).is_err());

        assert!(VerifyingKey::<MNT4>::read(vk_bls.as_slice()).is_err());
        assert!(VerifyingKey::<MNT4>::read(vk_sw6.as_slice()).is_err());
        assert!(VerifyingKey::<MNT4>::read(vk_mnt6.as_slice()).is_err());

        //MNT6-753
        assert!(Parameters::<MNT6>::read(params_bls.as_slice()).is_err());
        assert!(Parameters::<MNT6>::read(params_sw6.as_slice()).is_err());
        assert!(Parameters::<MNT6>::read(params_mnt4.as_slice()).is_err());

        assert!(Proof::<MNT6>::read(proof_bls.as_slice()).is_err());
        assert!(Proof::<MNT6>::read(proof_sw6.as_slice()).is_err());
        assert!(Proof::<MNT6>::read(proof_mnt4.as_slice()).is_err());

        assert!(VerifyingKey::<MNT6>::read(vk_bls.as_slice()).is_err());
        assert!(VerifyingKey::<MNT6>::read(vk_sw6.as_slice()).is_err());
        assert!(VerifyingKey::<MNT6>::read(vk_mnt4.as_slice()).is_err());
    }

    #[test]
    fn read_proof_checked_unchecked(){

        // BLS12-377
        let mut proof = vec![0u8; 387];
        assert!(Proof::<Bls12_377>::read_unchecked(proof.as_slice()).is_ok());
        assert!(Proof::<Bls12_377>::read(proof.as_slice()).is_err());

        // SW6
        proof = vec![0u8; 1083];
        assert!(Proof::<SW6>::read_unchecked(proof.as_slice()).is_ok());
        assert!(Proof::<SW6>::read(proof.as_slice()).is_err());

        // MNT4-753
        proof = vec![0u8; 771];
        assert!(Proof::<MNT4>::read_unchecked(proof.as_slice()).is_ok());
        assert!(Proof::<MNT4>::read(proof.as_slice()).is_err());

        // MNT6-753
        proof = vec![0u8; 963];
        Proof::<MNT6>::read_unchecked(proof.as_slice()).unwrap();
        assert!(Proof::<MNT6>::read_unchecked(proof.as_slice()).is_ok());
        assert!(Proof::<MNT6>::read(proof.as_slice()).is_err());
    }
}