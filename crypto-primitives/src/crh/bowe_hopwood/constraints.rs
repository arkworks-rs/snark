use core::{borrow::Borrow, marker::PhantomData};

use crate::{
    crh::{
        bowe_hopwood::{Parameters, CHUNK_SIZE, CRH},
        pedersen::Window,
        FixedLengthCRHGadget,
    },
    Vec,
};
use algebra_core::{
    curves::{ModelParameters, TEModelParameters},
    Field,
};
use r1cs_core::{Namespace, SynthesisError};
use r1cs_std::{
    alloc::AllocVar, groups::curves::twisted_edwards::AffineVar, prelude::*, uint8::UInt8,
};

use r1cs_std::bits::boolean::Boolean;

type ConstraintF<P> = <<P as ModelParameters>::BaseField as Field>::BasePrimeField;

#[derive(Derivative)]
#[derivative(Clone(bound = "P: TEModelParameters, W: Window"))]
pub struct ParametersVar<P: TEModelParameters, W: Window> {
    params: Parameters<P>,
    #[doc(hidden)]
    _window: PhantomData<W>,
}

pub struct CRHGadget<P: TEModelParameters, F: FieldVar<P::BaseField, ConstraintF<P>>>
where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[doc(hidden)]
    _params: PhantomData<P>,
    #[doc(hidden)]
    _base_field: PhantomData<F>,
}

impl<P, F, W> FixedLengthCRHGadget<CRH<P, W>, ConstraintF<P>> for CRHGadget<P, F>
where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
    F: FieldVar<P::BaseField, ConstraintF<P>>,
    F: TwoBitLookupGadget<ConstraintF<P>, TableConstant = P::BaseField>
        + ThreeBitCondNegLookupGadget<ConstraintF<P>, TableConstant = P::BaseField>,
    P: TEModelParameters,
    W: Window,
{
    type OutputVar = AffineVar<P, F>;
    type ParametersVar = ParametersVar<P, W>;

    #[tracing::instrument(target = "r1cs", skip(parameters, input))]
    fn evaluate(
        parameters: &Self::ParametersVar,
        input: &[UInt8<ConstraintF<P>>],
    ) -> Result<Self::OutputVar, SynthesisError> {
        // Pad the input if it is not the current length.
        let mut input_in_bits: Vec<Boolean<_>> = input
            .iter()
            .flat_map(|byte| byte.to_bits_le().unwrap())
            .collect();
        if (input_in_bits.len()) % CHUNK_SIZE != 0 {
            let current_length = input_in_bits.len();
            for _ in 0..(CHUNK_SIZE - current_length % CHUNK_SIZE) {
                input_in_bits.push(Boolean::constant(false));
            }
        }
        assert!(input_in_bits.len() % CHUNK_SIZE == 0);
        assert_eq!(parameters.params.generators.len(), W::NUM_WINDOWS);
        for generators in parameters.params.generators.iter() {
            assert_eq!(generators.len(), W::WINDOW_SIZE);
        }

        // Allocate new variable for the result.
        let input_in_bits = input_in_bits
            .chunks(W::WINDOW_SIZE * CHUNK_SIZE)
            .map(|x| x.chunks(CHUNK_SIZE).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let result = AffineVar::precomputed_base_3_bit_signed_digit_scalar_mul(
            &parameters.params.generators,
            &input_in_bits,
        )?;

        Ok(result)
    }
}

impl<P, W> AllocVar<Parameters<P>, ConstraintF<P>> for ParametersVar<P, W>
where
    P: TEModelParameters,
    W: Window,
{
    #[tracing::instrument(target = "r1cs", skip(_cs, f))]
    fn new_variable<T: Borrow<Parameters<P>>>(
        _cs: impl Into<Namespace<ConstraintF<P>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        _mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let params = f()?.borrow().clone();
        Ok(ParametersVar {
            params,
            _window: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use rand::Rng;

    use crate::crh::{
        bowe_hopwood::{constraints::CRHGadget, CRH},
        pedersen::Window as PedersenWindow,
        FixedLengthCRH, FixedLengthCRHGadget,
    };
    use algebra::{
        ed_on_bls12_381::{EdwardsParameters, Fq as Fr},
        test_rng, ProjectiveCurve,
    };
    use r1cs_core::{ConstraintSystem, ConstraintSystemRef};
    use r1cs_std::{alloc::AllocVar, ed_on_bls12_381::FqVar, uint8::UInt8, R1CSVar};

    type TestCRH = CRH<EdwardsParameters, Window>;
    type TestCRHGadget = CRHGadget<EdwardsParameters, FqVar>;

    #[derive(Clone, PartialEq, Eq, Hash)]
    pub(super) struct Window;

    impl PedersenWindow for Window {
        const WINDOW_SIZE: usize = 63;
        const NUM_WINDOWS: usize = 8;
    }

    fn generate_input<R: Rng>(
        cs: ConstraintSystemRef<Fr>,
        rng: &mut R,
    ) -> ([u8; 189], Vec<UInt8<Fr>>) {
        let mut input = [1u8; 189];
        rng.fill_bytes(&mut input);

        let mut input_bytes = vec![];
        for byte in input.iter() {
            input_bytes.push(UInt8::new_witness(cs.clone(), || Ok(byte)).unwrap());
        }
        (input, input_bytes)
    }

    #[test]
    fn test_native_equality() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<Fr>::new_ref();

        let (input, input_var) = generate_input(cs.clone(), rng);
        println!("number of constraints for input: {}", cs.num_constraints());

        let parameters = TestCRH::setup(rng).unwrap();
        let primitive_result = TestCRH::evaluate(&parameters, &input).unwrap();

        let parameters_var =
            <TestCRHGadget as FixedLengthCRHGadget<TestCRH, Fr>>::ParametersVar::new_witness(
                r1cs_core::ns!(cs, "parameters_var"),
                || Ok(&parameters),
            )
            .unwrap();
        println!(
            "number of constraints for input + params: {}",
            cs.num_constraints()
        );

        let result_var = TestCRHGadget::evaluate(&parameters_var, &input_var).unwrap();

        println!("number of constraints total: {}", cs.num_constraints());

        let primitive_result = primitive_result.into_affine();
        assert_eq!(primitive_result, result_var.value().unwrap().into_affine());
        assert!(cs.is_satisfied().unwrap());
    }
}
