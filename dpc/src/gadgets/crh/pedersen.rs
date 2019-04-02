use algebra::PairingEngine;
use std::hash::Hash;

use crate::{
    crypto_primitives::crh::pedersen::{PedersenCRH, PedersenParameters, PedersenWindow},
    gadgets::crh::FixedLengthCRHGadget,
};
use algebra::groups::Group;
use snark::{ConstraintSystem, SynthesisError};
use snark_gadgets::{groups::GroupGadget, uint8::UInt8, utils::AllocGadget};

use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(Clone(
    bound = "G: Group, W: PedersenWindow, E: PairingEngine, GG: GroupGadget<G, E>"
))]
pub struct PedersenCRHGadgetParameters<
    G: Group,
    W: PedersenWindow,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
> {
    params:   PedersenParameters<G>,
    _group_g: PhantomData<GG>,
    _engine:  PhantomData<E>,
    _window:  PhantomData<W>,
}

pub struct PedersenCRHGadget<G: Group, E: PairingEngine, GG: GroupGadget<G, E>> {
    _group:        PhantomData<*const G>,
    _group_gadget: PhantomData<*const GG>,
    _engine:       PhantomData<E>,
}

impl<E, G, GG, W> FixedLengthCRHGadget<PedersenCRH<G, W>, E> for PedersenCRHGadget<G, E, GG>
where
    E: PairingEngine,
    G: Group + Hash,
    GG: GroupGadget<G, E>,
    W: PedersenWindow,
{
    type OutputGadget = GG;
    type ParametersGadget = PedersenCRHGadgetParameters<G, W, E, GG>;

    fn check_evaluation_gadget<CS: ConstraintSystem<E>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError> {
        let mut padded_input = input.to_vec();
        // Pad the input if it is not the current length.
        if input.len() * 8 < W::WINDOW_SIZE * W::NUM_WINDOWS {
            let current_length = input.len();
            for _ in current_length..(W::WINDOW_SIZE * W::NUM_WINDOWS / 8) {
                padded_input.push(UInt8::constant(0u8));
            }
        }
        assert_eq!(padded_input.len() * 8, W::WINDOW_SIZE * W::NUM_WINDOWS);
        assert_eq!(parameters.params.generators.len(), W::NUM_WINDOWS);

        // Allocate new variable for the result.
        let input_in_bits: Vec<_> = padded_input
            .iter()
            .flat_map(|byte| byte.into_bits_le())
            .collect();
        let input_in_bits = input_in_bits.chunks(W::WINDOW_SIZE);
        let result =
            GG::precomputed_base_multiscalar_mul(cs, &parameters.params.generators, input_in_bits)?;

        Ok(result)
    }

    fn cost() -> usize {
        use snark_gadgets::utils::CondSelectGadget;
        W::NUM_WINDOWS * W::WINDOW_SIZE * (GG::cost_of_add() + <GG as CondSelectGadget<E>>::cost())
    }
}

impl<G: Group, W: PedersenWindow, E: PairingEngine, GG: GroupGadget<G, E>>
    AllocGadget<PedersenParameters<G>, E> for PedersenCRHGadgetParameters<G, W, E, GG>
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(_cs: CS, value_gen: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PedersenParameters<G>>,
    {
        let params = value_gen()?.borrow().clone();
        Ok(PedersenCRHGadgetParameters {
            params,
            _group_g: PhantomData,
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
        let params = value_gen()?.borrow().clone();
        Ok(PedersenCRHGadgetParameters {
            params,
            _group_g: PhantomData,
            _engine: PhantomData,
            _window: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use algebra::curves::bls12_381::Bls12_381;
    use rand::{thread_rng, Rng};

    use crate::{
        crypto_primitives::crh::{
            pedersen::{PedersenCRH, PedersenWindow},
            FixedLengthCRH,
        },
        gadgets::crh::{pedersen::PedersenCRHGadget, FixedLengthCRHGadget},
    };
    use algebra::curves::{jubjub::JubJubProjective as JubJub, ProjectiveCurve};
    use snark::ConstraintSystem;
    use snark_gadgets::{
        groups::curves::twisted_edwards::jubjub::JubJubGadget,
        test_constraint_system::TestConstraintSystem, uint8::UInt8, utils::AllocGadget,
    };

    type TestCRH = PedersenCRH<JubJub, Window>;
    type TestCRHGadget = PedersenCRHGadget<JubJub, Bls12_381, JubJubGadget>;

    #[derive(Clone, PartialEq, Eq, Hash)]
    pub(super) struct Window;

    impl PedersenWindow for Window {
        const WINDOW_SIZE: usize = 128;
        const NUM_WINDOWS: usize = 8;
    }

    #[test]
    fn num_constraints() {
        let rng = &mut thread_rng();
        let mut cs = TestConstraintSystem::<Bls12_381>::new();

        let (_input, input_bytes) = generate_input(&mut cs, rng);
        let input_constraints = cs.num_constraints();
        println!("number of constraints for input: {}", cs.num_constraints());

        let parameters = TestCRH::setup(rng).unwrap();

        let gadget_parameters =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let param_constraints = cs.num_constraints() - input_constraints;
        println!(
            "number of constraints for input + params: {}",
            cs.num_constraints()
        );

        let _ =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::check_evaluation_gadget(
                &mut cs.ns(|| "gadget_evaluation"),
                &gadget_parameters,
                &input_bytes,
            )
            .unwrap();

        println!("number of constraints total: {}", cs.num_constraints());
        let eval_constraints = cs.num_constraints() - param_constraints - input_constraints;
        assert_eq!(
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::cost(),
            eval_constraints
        );
        assert_eq!(
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::cost(),
            256 * (6 + 4 * (6 + 2))
        );
    }

    fn generate_input<CS: ConstraintSystem<Bls12_381>>(
        mut cs: CS,
        rng: &mut dyn Rng,
    ) -> ([u8; 128], Vec<UInt8>) {
        let mut input = [1u8; 128];
        rng.fill_bytes(&mut input);

        let mut input_bytes = vec![];
        for (byte_i, input_byte) in input.into_iter().enumerate() {
            let cs = cs.ns(|| format!("input_byte_gadget_{}", byte_i));
            input_bytes.push(UInt8::alloc(cs, || Ok(*input_byte)).unwrap());
        }
        (input, input_bytes)
    }

    #[test]
    fn crh_primitive_gadget_test() {
        let rng = &mut thread_rng();
        let mut cs = TestConstraintSystem::<Bls12_381>::new();

        let (input, input_bytes) = generate_input(&mut cs, rng);
        println!("number of constraints for input: {}", cs.num_constraints());

        let parameters = TestCRH::setup(rng).unwrap();
        let primitive_result = TestCRH::evaluate(&parameters, &input).unwrap();

        let gadget_parameters =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        println!(
            "number of constraints for input + params: {}",
            cs.num_constraints()
        );

        let gadget_result =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Bls12_381>>::check_evaluation_gadget(
                &mut cs.ns(|| "gadget_evaluation"),
                &gadget_parameters,
                &input_bytes,
            )
            .unwrap();

        println!("number of constraints total: {}", cs.num_constraints());

        let primitive_result = primitive_result.into_affine();
        assert_eq!(primitive_result.x, gadget_result.x.value.unwrap());
        assert_eq!(primitive_result.y, gadget_result.y.value.unwrap());
        assert!(cs.is_satisfied());
    }
}
