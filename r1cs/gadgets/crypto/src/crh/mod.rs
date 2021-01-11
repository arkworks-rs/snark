use algebra::{
    Field, PrimeField
};
use std::fmt::Debug;

use primitives::crh::{
    FieldBasedHash, FixedLengthCRH
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use r1cs_std::prelude::*;

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;

pub mod sbox;
pub use self::sbox::*;

pub mod poseidon;
pub use self::poseidon::*;
use primitives::AlgebraicSponge;

pub trait FixedLengthCRHGadget<H: FixedLengthCRH, ConstraintF: Field>: Sized {
    type OutputGadget: EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + CondSelectGadget<ConstraintF>
        + AllocGadget<H::Output, ConstraintF>
        + Debug
        + Clone
        + Sized;
    type ParametersGadget: AllocGadget<H::Parameters, ConstraintF> + Clone;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError>;
}

pub trait FieldBasedHashGadget<H: FieldBasedHash<Data = ConstraintF>, ConstraintF: Field>: Sized {
    type DataGadget: FieldGadget<ConstraintF, ConstraintF>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>;
}

pub trait FieldHasherGadget<
    H: FieldBasedHash<Data = ConstraintF>,
    ConstraintF: Field,
    HG: FieldBasedHashGadget<H, ConstraintF>
>
{
    fn enforce_hash<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        personalization: Option<&[HG::DataGadget]>
    ) -> Result<HG::DataGadget, SynthesisError>;
}

pub trait AlgebraicSpongeGadget<H: AlgebraicSponge<ConstraintF>, ConstraintF: PrimeField>: Sized {
    type DataGadget: FieldGadget<ConstraintF, ConstraintF>;

    fn new<CS: ConstraintSystem<ConstraintF>>(cs: CS) -> Result<Self, SynthesisError>;

    fn enforce_absorb<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<(), SynthesisError>;

    fn enforce_squeeze<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        num: usize,
    ) -> Result<Vec<Self::DataGadget>, SynthesisError>;
}

#[cfg(test)]
mod test {
    use algebra::PrimeField;
    use primitives::{
        FieldBasedHash, AlgebraicSponge
    };
    use crate::{
        FieldBasedHashGadget, AlgebraicSpongeGadget,
    };
    use r1cs_std::{
        fields::fp::FpGadget,
        test_constraint_system::TestConstraintSystem,
        alloc::AllocGadget,
    };
    use r1cs_core::ConstraintSystem;

    pub(crate) fn field_based_hash_gadget_native_test<
        F: PrimeField,
        H: FieldBasedHash<Data = F>,
        HG: FieldBasedHashGadget<H, F, DataGadget = FpGadget<F>>
    >(inputs: Vec<F>)
    {
        let mut cs = TestConstraintSystem::<F>::new();

        let primitive_result = {
            let mut digest = H::init(None);
            inputs.iter().for_each(|elem| { digest.update(*elem);});
            digest.finalize()
        };

        let mut input_gadgets = Vec::with_capacity(inputs.len());
        inputs.into_iter().enumerate().for_each(|(i, elem)|{
            let elem_gadget = HG::DataGadget::alloc(
                cs.ns(|| format!("alloc input {}", i)),
                || Ok(elem)
            ).unwrap();
            input_gadgets.push(elem_gadget);
        });

        let gadget_result = HG::check_evaluation_gadget(
            cs.ns(||"check_poseidon_gadget"),
            input_gadgets.as_slice()
        ).unwrap();

        assert_eq!(primitive_result, gadget_result.value.unwrap());
        assert!(cs.is_satisfied());
    }

    pub(crate) fn algebraic_sponge_gadget_native_test<
        F: PrimeField,
        H: AlgebraicSponge<F>,
        HG: AlgebraicSpongeGadget<H, F, DataGadget = FpGadget<F>>
    >(inputs: Vec<F>)
    {
        let mut cs = TestConstraintSystem::<F>::new();

        let mut primitive_sponge = H::new();
        primitive_sponge.absorb(inputs.clone());

        let mut input_gadgets = Vec::with_capacity(inputs.len());
        inputs.iter().enumerate().for_each(|(i, elem)|{
            let elem_gadget = HG::DataGadget::alloc(
                cs.ns(|| format!("alloc input {}", i)),
                || Ok(elem.clone())
            ).unwrap();
            input_gadgets.push(elem_gadget);
        });
        
        let mut sponge_gadget = HG::new(cs.ns(|| "new poseidon sponge")).unwrap();
        sponge_gadget.enforce_absorb(cs.ns(|| "absorb inputs"), input_gadgets.as_slice()).unwrap();

        for i in 0..inputs.len() {
            let output_gadgets = sponge_gadget.enforce_squeeze(
                cs.ns(|| format!("squeeze {} field elements",  i + 1)),
                i + 1
            ).unwrap().iter().map(|fe_gadget| fe_gadget.value.unwrap()).collect::<Vec<_>>();
            assert_eq!(output_gadgets, primitive_sponge.squeeze(i + 1));
        }

        assert!(cs.is_satisfied());
    }
}