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

        Ok(())
    }
}

mod bls12_377 {
    use super::*;
    use crate::{
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use algebra::{UniformRand, curves::bls12_377::Bls12_377, fields::bls12_377::Fr};
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
}

mod sw6 {
    use super::*;
    use crate::{
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use rand::thread_rng;

    use algebra::{UniformRand, curves::sw6::SW6, fields::sw6::Fr as SW6Fr, Field};

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
}
