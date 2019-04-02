use crate::crypto_primitives::{
    commitment::pedersen::{PedersenCommitment, PedersenParameters, PedersenRandomness},
    crh::pedersen::PedersenWindow,
};
use algebra::{to_bytes, Group, ToBytes};
use snark::{ConstraintSystem, SynthesisError};

use crate::gadgets::CommitmentGadget;
use algebra::PairingEngine;
use snark_gadgets::{groups::GroupGadget, uint8::UInt8, utils::AllocGadget};
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(Clone(bound = "G: Group, W: PedersenWindow, E: PairingEngine"))]
pub struct PedersenCommitmentGadgetParameters<G: Group, W: PedersenWindow, E: PairingEngine> {
    params:  PedersenParameters<G>,
    _group:  PhantomData<G>,
    _engine: PhantomData<E>,
    _window: PhantomData<W>,
}

#[derive(Clone, Debug)]
pub struct PedersenRandomnessGadget(Vec<UInt8>);

pub struct PedersenCommitmentGadget<G: Group, E: PairingEngine, GG: GroupGadget<G, E>>(
    PhantomData<*const G>,
    PhantomData<*const GG>,
    PhantomData<E>,
);

impl<E, G, GG, W> CommitmentGadget<PedersenCommitment<G, W>, E>
    for PedersenCommitmentGadget<G, E, GG>
where
    E: PairingEngine,
    G: Group,
    GG: GroupGadget<G, E>,
    W: PedersenWindow,
{
    type OutputGadget = GG;
    type ParametersGadget = PedersenCommitmentGadgetParameters<G, W, E>;
    type RandomnessGadget = PedersenRandomnessGadget;

    fn check_commitment_gadget<CS: ConstraintSystem<E>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
        r: &Self::RandomnessGadget,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        assert!((input.len() * 8) <= (W::WINDOW_SIZE * W::NUM_WINDOWS));

        let mut padded_input = input.to_vec();
        // Pad if input length is less than `W::WINDOW_SIZE * W::NUM_WINDOWS`.
        if (input.len() * 8) < W::WINDOW_SIZE * W::NUM_WINDOWS {
            let current_length = input.len();
            for _ in current_length..((W::WINDOW_SIZE * W::NUM_WINDOWS) / 8) {
                padded_input.push(UInt8::constant(0u8));
            }
        }

        assert_eq!(padded_input.len() * 8, W::WINDOW_SIZE * W::NUM_WINDOWS);
        assert_eq!(parameters.params.generators.len(), W::NUM_WINDOWS);

        // Allocate new variable for commitment output.
        let input_in_bits: Vec<_> = padded_input
            .iter()
            .flat_map(|byte| byte.into_bits_le())
            .collect();
        let input_in_bits = input_in_bits.chunks(W::WINDOW_SIZE);
        let mut result = GG::precomputed_base_multiscalar_mul(
            cs.ns(|| "multiexp"),
            &parameters.params.generators,
            input_in_bits,
        )?;

        // Compute h^r
        let rand_bits: Vec<_> = r.0.iter().flat_map(|byte| byte.into_bits_le()).collect();
        result.precomputed_base_scalar_mul(
            cs.ns(|| "Randomizer"),
            rand_bits
                .iter()
                .zip(&parameters.params.randomness_generator),
        )?;

        Ok(result)
    }
}

impl<G, W, E> AllocGadget<PedersenParameters<G>, E> for PedersenCommitmentGadgetParameters<G, W, E>
where
    G: Group,
    W: PedersenWindow,
    E: PairingEngine,
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(_cs: CS, value_gen: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PedersenParameters<G>>,
    {
        let temp = value_gen()?;
        let parameters = temp.borrow().clone();

        Ok(PedersenCommitmentGadgetParameters {
            params:  parameters,
            _group:  PhantomData,
            _engine: PhantomData,
            _window: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(
        _cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PedersenParameters<G>>,
    {
        let temp = value_gen()?;
        let parameters = temp.borrow().clone();

        Ok(PedersenCommitmentGadgetParameters {
            params:  parameters,
            _group:  PhantomData,
            _engine: PhantomData,
            _window: PhantomData,
        })
    }
}

impl<G, E> AllocGadget<PedersenRandomness<G>, E> for PedersenRandomnessGadget
where
    G: Group,
    E: PairingEngine,
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(cs: CS, value_gen: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PedersenRandomness<G>>,
    {
        let temp = value_gen()?;
        let randomness = to_bytes![temp.borrow().0].unwrap();
        Ok(PedersenRandomnessGadget(UInt8::alloc_vec(cs, &randomness)?))
    }

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(
        cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PedersenRandomness<G>>,
    {
        let temp = value_gen()?;
        let randomness = to_bytes![temp.borrow().0].unwrap();
        Ok(PedersenRandomnessGadget(UInt8::alloc_input_vec(
            cs,
            &randomness,
        )?))
    }
}

#[cfg(test)]
mod test {
    use algebra::{curves::bls12_381::Bls12_381, fields::jubjub::fr::Fr};
    use rand::{thread_rng, Rand};

    use crate::{
        crypto_primitives::{
            commitment::{
                pedersen::{PedersenCommitment, PedersenRandomness},
                CommitmentScheme,
            },
            crh::pedersen::PedersenWindow,
        },
        gadgets::commitment::{pedersen::PedersenCommitmentGadget, CommitmentGadget},
    };
    use algebra::curves::{jubjub::JubJubProjective as JubJub, ProjectiveCurve};
    use snark::ConstraintSystem;
    use snark_gadgets::{
        groups::jubjub::JubJubGadget, test_constraint_system::TestConstraintSystem, uint8::UInt8,
        utils::AllocGadget,
    };

    #[test]
    fn commitment_gadget_test() {
        let mut cs = TestConstraintSystem::<Bls12_381>::new();

        #[derive(Clone, PartialEq, Eq, Hash)]
        pub(super) struct Window;

        impl PedersenWindow for Window {
            const WINDOW_SIZE: usize = 4;
            const NUM_WINDOWS: usize = 8;
        }

        let input = [1u8; 4];

        let rng = &mut thread_rng();

        type TestCOMM = PedersenCommitment<JubJub, Window>;
        type TestCOMMGadget = PedersenCommitmentGadget<JubJub, Bls12_381, JubJubGadget>;

        let randomness = PedersenRandomness(Fr::rand(rng));

        let parameters = PedersenCommitment::<JubJub, Window>::setup(rng).unwrap();
        let primitive_result =
            PedersenCommitment::<JubJub, Window>::commit(&parameters, &input, &randomness).unwrap();

        let mut input_bytes = vec![];
        for (byte_i, input_byte) in input.into_iter().enumerate() {
            let cs = cs.ns(|| format!("input_byte_gadget_{}", byte_i));
            input_bytes.push(UInt8::alloc(cs, || Ok(*input_byte)).unwrap());
        }

        let randomness =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Bls12_381>>::RandomnessGadget::alloc(
                &mut cs.ns(|| "gadget_randomness"),
                || Ok(&randomness),
            )
            .unwrap();
        let gadget_parameters =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Bls12_381>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let gadget_result =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Bls12_381>>::check_commitment_gadget(
                &mut cs.ns(|| "gadget_evaluation"),
                &gadget_parameters,
                &input_bytes,
                &randomness,
            )
            .unwrap();

        let primitive_result = primitive_result.into_affine();
        assert_eq!(primitive_result.x, gadget_result.x.value.unwrap());
        assert_eq!(primitive_result.y, gadget_result.y.value.unwrap());
        assert!(cs.is_satisfied());
    }
}
