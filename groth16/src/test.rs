use algebra_core::Field;
use r1cs_core::{lc, ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
struct MySillyCircuit<F: Field> {
    a: Option<F>,
    b: Option<F>,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for MySillyCircuit<ConstraintF> {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_input_variable(|| {
            let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

            a.mul_assign(&b);
            Ok(a)
        })?;

        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;

        Ok(())
    }
}

mod bls12_377 {
    use super::*;
    use crate::{
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use algebra::bls12_377::{Bls12_377, Fr};
    use algebra_core::{test_rng, UniformRand};
    use core::ops::MulAssign;

    #[test]
    fn prove_and_verify() {
        let rng = &mut test_rng();

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

mod cp6_782 {
    use super::*;
    use crate::{
        create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    };

    use algebra::cp6_782::{Fr, CP6_782};
    use algebra_core::{test_rng, UniformRand, Zero};

    #[test]
    fn prove_and_verify() {
        let rng = &mut test_rng();

        let params =
            generate_random_parameters::<CP6_782, _, _>(MySillyCircuit { a: None, b: None }, rng)
                .unwrap();

        let pvk = prepare_verifying_key::<CP6_782>(&params.vk);

        let a = Fr::rand(rng);
        let b = Fr::rand(rng);
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
        assert!(!verify_proof(&pvk, &proof, &[Fr::zero()]).unwrap());
    }
}
