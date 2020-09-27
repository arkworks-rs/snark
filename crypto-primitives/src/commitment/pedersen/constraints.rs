use crate::{
    commitment::pedersen::{Commitment, Parameters, Randomness},
    crh::pedersen::Window,
    Vec,
};
use algebra_core::{
    fields::{Field, PrimeField},
    to_bytes, ProjectiveCurve, Zero,
};
use r1cs_core::{Namespace, SynthesisError};

use core::{borrow::Borrow, marker::PhantomData};
use r1cs_std::prelude::*;

type ConstraintF<C> = <<C as ProjectiveCurve>::BaseField as Field>::BasePrimeField;

#[derive(Derivative)]
#[derivative(Clone(bound = "C: ProjectiveCurve, GG: CurveVar<C, ConstraintF<C>>"))]
pub struct ParametersVar<C: ProjectiveCurve, GG: CurveVar<C, ConstraintF<C>>>
where
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
{
    params: Parameters<C>,
    #[doc(hidden)]
    _group_var: PhantomData<GG>,
}

#[derive(Clone, Debug)]
pub struct RandomnessVar<F: Field>(Vec<UInt8<F>>);

pub struct CommGadget<C: ProjectiveCurve, GG: CurveVar<C, ConstraintF<C>>, W: Window>
where
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
{
    #[doc(hidden)]
    _curve: PhantomData<*const C>,
    #[doc(hidden)]
    _group_var: PhantomData<*const GG>,
    #[doc(hidden)]
    _window: PhantomData<*const W>,
}

impl<C, GG, W> crate::commitment::CommitmentGadget<Commitment<C, W>, ConstraintF<C>>
    for CommGadget<C, GG, W>
where
    C: ProjectiveCurve,
    GG: CurveVar<C, ConstraintF<C>>,
    W: Window,
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
    ConstraintF<C>: PrimeField,
{
    type OutputVar = GG;
    type ParametersVar = ParametersVar<C, GG>;
    type RandomnessVar = RandomnessVar<ConstraintF<C>>;

    #[tracing::instrument(target = "r1cs", skip(parameters, r))]
    fn commit(
        parameters: &Self::ParametersVar,
        input: &[UInt8<ConstraintF<C>>],
        r: &Self::RandomnessVar,
    ) -> Result<Self::OutputVar, SynthesisError> {
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
        let input_in_bits: Vec<Boolean<_>> = padded_input
            .iter()
            .flat_map(|byte| byte.to_bits_le().unwrap())
            .collect();
        let input_in_bits = input_in_bits.chunks(W::WINDOW_SIZE);
        let mut result =
            GG::precomputed_base_multiscalar_mul_le(&parameters.params.generators, input_in_bits)?;

        // Compute h^r
        let rand_bits: Vec<_> =
            r.0.iter()
                .flat_map(|byte| byte.to_bits_le().unwrap())
                .collect();
        result.precomputed_base_scalar_mul_le(
            rand_bits
                .iter()
                .zip(&parameters.params.randomness_generator),
        )?;

        Ok(result)
    }
}

impl<C, GG> AllocVar<Parameters<C>, ConstraintF<C>> for ParametersVar<C, GG>
where
    C: ProjectiveCurve,
    GG: CurveVar<C, ConstraintF<C>>,
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
{
    fn new_variable<T: Borrow<Parameters<C>>>(
        _cs: impl Into<Namespace<ConstraintF<C>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        _mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let params = f()?.borrow().clone();
        Ok(ParametersVar {
            params,
            _group_var: PhantomData,
        })
    }
}

impl<C, F> AllocVar<Randomness<C>, F> for RandomnessVar<F>
where
    C: ProjectiveCurve,
    F: PrimeField,
{
    fn new_variable<T: Borrow<Randomness<C>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let r = to_bytes![&f().map(|b| b.borrow().0).unwrap_or(C::ScalarField::zero())].unwrap();
        match mode {
            AllocationMode::Constant => Ok(Self(UInt8::constant_vec(&r))),
            AllocationMode::Input => UInt8::new_input_vec(cs, &r).map(Self),
            AllocationMode::Witness => UInt8::new_witness_vec(cs, &r).map(Self),
        }
    }
}

#[cfg(test)]
mod test {
    use algebra::{
        ed_on_bls12_381::{EdwardsProjective as JubJub, Fq, Fr},
        test_rng, UniformRand,
    };

    use crate::{
        commitment::{
            pedersen::{constraints::CommGadget, Commitment, Randomness},
            CommitmentGadget, CommitmentScheme,
        },
        crh::pedersen,
    };
    use r1cs_core::ConstraintSystem;
    use r1cs_std::{ed_on_bls12_381::EdwardsVar, prelude::*};

    #[test]
    fn commitment_gadget_test() {
        let cs = ConstraintSystem::<Fq>::new_ref();

        #[derive(Clone, PartialEq, Eq, Hash)]
        pub(super) struct Window;

        impl pedersen::Window for Window {
            const WINDOW_SIZE: usize = 4;
            const NUM_WINDOWS: usize = 8;
        }

        let input = [1u8; 4];

        let rng = &mut test_rng();

        type TestCOMM = Commitment<JubJub, Window>;
        type TestCOMMGadget = CommGadget<JubJub, EdwardsVar, Window>;

        let randomness = Randomness(Fr::rand(rng));

        let parameters = Commitment::<JubJub, Window>::setup(rng).unwrap();
        let primitive_result =
            Commitment::<JubJub, Window>::commit(&parameters, &input, &randomness).unwrap();

        let mut input_var = vec![];
        for input_byte in input.iter() {
            input_var.push(UInt8::new_witness(cs.clone(), || Ok(*input_byte)).unwrap());
        }

        let randomness_var =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Fq>>::RandomnessVar::new_witness(
                r1cs_core::ns!(cs, "gadget_randomness"),
                || Ok(&randomness),
            )
            .unwrap();
        let parameters_var =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Fq>>::ParametersVar::new_witness(
                r1cs_core::ns!(cs, "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let result_var =
            TestCOMMGadget::commit(&parameters_var, &input_var, &randomness_var).unwrap();

        let primitive_result = primitive_result;
        assert_eq!(primitive_result, result_var.value().unwrap());
        assert!(cs.is_satisfied().unwrap());
    }
}
