use crate::crh::FixedLengthCRHGadget;
use algebra::Field;
use std::hash::Hash;

use algebra::groups::Group;
use primitives::{
    bowe_hopwood::{BoweHopwoodPedersenCRH, BoweHopwoodPedersenParameters, CHUNK_SIZE},
    crh::pedersen::PedersenWindow,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{alloc::AllocGadget, groups::GroupGadget, uint8::UInt8};

use r1cs_std::bits::boolean::Boolean;
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(Clone(
    bound = "G: Group, W: PedersenWindow, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"
))]
pub struct BoweHopwoodPedersenCRHGadgetParameters<
    G: Group,
    W: PedersenWindow,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
> {
    params: BoweHopwoodPedersenParameters<G>,
    _group_g: PhantomData<GG>,
    _engine: PhantomData<ConstraintF>,
    _window: PhantomData<W>,
}

pub struct BoweHopwoodPedersenCRHGadget<
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
> {
    _group: PhantomData<*const G>,
    _group_gadget: PhantomData<*const GG>,
    _engine: PhantomData<ConstraintF>,
}

impl<ConstraintF, G, GG, W> FixedLengthCRHGadget<BoweHopwoodPedersenCRH<G, W>, ConstraintF>
    for BoweHopwoodPedersenCRHGadget<G, ConstraintF, GG>
where
    ConstraintF: Field,
    G: Group + Hash,
    GG: GroupGadget<G, ConstraintF>,
    W: PedersenWindow,
{
    type OutputGadget = GG;
    type ParametersGadget = BoweHopwoodPedersenCRHGadgetParameters<G, W, ConstraintF, GG>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError> {
        // Pad the input if it is not the current length.
        let mut input_in_bits: Vec<_> = input.iter().flat_map(|byte| byte.into_bits_le()).collect();
        if (input_in_bits.len()) % CHUNK_SIZE != 0 {
            let current_length = input_in_bits.len();
            for _ in 0..(CHUNK_SIZE - current_length % CHUNK_SIZE) {
                input_in_bits.push(Boolean::constant(false));
            }
        }
        if input_in_bits.len() % CHUNK_SIZE != 0 {
            return Err(SynthesisError::Other(format!(
                "Input is not multiple of the chunk size. Input len: {}, chunk size: {}",
                input_in_bits.len(),
                CHUNK_SIZE,
            )));
        }
        if parameters.params.generators.len() != W::NUM_WINDOWS {
            return Err(SynthesisError::Other(format!(
                "Incorrect pp of size {:?} for window params {:?}x{:?}x{}",
                parameters.params.generators.len(),
                W::WINDOW_SIZE,
                W::NUM_WINDOWS,
                CHUNK_SIZE
            )));
        }
        for generators in parameters.params.generators.iter() {
            if generators.len() != W::WINDOW_SIZE {
                return Err(SynthesisError::Other(format!(
                    "Number of generators: {} not enough for the selected window size: {}",
                    parameters.params.generators.len(),
                    W::WINDOW_SIZE
                )));
            }
        }

        // Allocate new variable for the result.
        let input_in_bits = input_in_bits
            .chunks(W::WINDOW_SIZE * CHUNK_SIZE)
            .map(|x| x.chunks(CHUNK_SIZE).into_iter().collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let result = GG::precomputed_base_3_bit_signed_digit_scalar_mul(
            cs,
            &parameters.params.generators,
            &input_in_bits,
        )?;

        Ok(result)
    }
}

impl<G: Group, W: PedersenWindow, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>>
    AllocGadget<BoweHopwoodPedersenParameters<G>, ConstraintF>
    for BoweHopwoodPedersenCRHGadgetParameters<G, W, ConstraintF, GG>
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        _cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<BoweHopwoodPedersenParameters<G>>,
    {
        let params = value_gen()?.borrow().clone();
        Ok(BoweHopwoodPedersenCRHGadgetParameters {
            params,
            _group_g: PhantomData,
            _engine: PhantomData,
            _window: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        _cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<BoweHopwoodPedersenParameters<G>>,
    {
        let params = value_gen()?.borrow().clone();
        Ok(BoweHopwoodPedersenCRHGadgetParameters {
            params,
            _group_g: PhantomData,
            _engine: PhantomData,
            _window: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::crh::{bowe_hopwood::BoweHopwoodPedersenCRHGadget, FixedLengthCRHGadget};
    use algebra::fields::sw6::fr::Fr;
    use algebra::{curves::edwards_sw6::EdwardsProjective as Edwards, ProjectiveCurve};
    use primitives::crh::{
        bowe_hopwood::BoweHopwoodPedersenCRH, pedersen::PedersenWindow, FixedLengthCRH,
    };
    use r1cs_core::ConstraintSystem;
    use r1cs_std::{
        alloc::AllocGadget, instantiated::edwards_sw6::EdwardsSWGadget,
        test_constraint_system::TestConstraintSystem, uint8::UInt8,
    };
    use rand::{thread_rng, Rng};

    type TestCRH = BoweHopwoodPedersenCRH<Edwards, Window>;
    type TestCRHGadget = BoweHopwoodPedersenCRHGadget<Edwards, Fr, EdwardsSWGadget>;

    #[derive(Clone, PartialEq, Eq, Hash)]
    pub(super) struct Window;

    impl PedersenWindow for Window {
        const WINDOW_SIZE: usize = 90;
        const NUM_WINDOWS: usize = 8;
    }

    fn generate_input<CS: ConstraintSystem<Fr>, R: Rng>(
        mut cs: CS,
        rng: &mut R,
    ) -> ([u8; 270], Vec<UInt8>) {
        let mut input = [1u8; 270];
        rng.fill_bytes(&mut input);

        let mut input_bytes = vec![];
        for (byte_i, input_byte) in input.iter().enumerate() {
            let cs = cs.ns(|| format!("input_byte_gadget_{}", byte_i));
            input_bytes.push(UInt8::alloc(cs, || Ok(*input_byte)).unwrap());
        }
        (input, input_bytes)
    }

    #[test]
    fn crh_primitive_gadget_test() {
        let rng = &mut thread_rng();
        let mut cs = TestConstraintSystem::<Fr>::new();

        let (input, input_bytes) = generate_input(&mut cs, rng);
        println!("number of constraints for input: {}", cs.num_constraints());

        let parameters = TestCRH::setup(rng).unwrap();
        let primitive_result = TestCRH::evaluate(&parameters, &input).unwrap();

        let gadget_parameters =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Fr>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        println!(
            "number of constraints for input + params: {}",
            cs.num_constraints()
        );

        let gadget_result =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Fr>>::check_evaluation_gadget(
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
