use algebra::Field;
use std::fmt::Debug;

use crate::crh::FixedLengthCRH;
use r1cs_core::{ConstraintSystem, SynthesisError};

use r1cs_std::prelude::*;

pub trait FixedLengthCRHGadget<H: FixedLengthCRH, ConstraintF: Field>: Sized {
    type OutputGadget: ConditionalEqGadget<ConstraintF>
        + EqGadget<ConstraintF>
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

//Temporary mock of Poseidon interfaces
use algebra::PrimeField;
use crate::crh::{PoseidonHash, FieldBasedHash, FieldBasedHashParameters};
use r1cs_std::fields::fp::FpGadget;
use std::borrow::Borrow;
use std::marker::PhantomData;

pub trait FieldBasedHashGadget<H: FieldBasedHash<Data = ConstraintF>, ConstraintF: Field>: Sized {
    type DataGadget: FieldGadget<ConstraintF, ConstraintF>;
    type ParametersGadget: AllocGadget<H::Parameters, ConstraintF> + Clone;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>;
}

pub struct PoseidonHashGadget
<
    ConstraintF: Field,
    P:           FieldBasedHashParameters<Fr = ConstraintF>,
>
{
    _field:      PhantomData<ConstraintF>,
    _parameters: PhantomData<P>,
}

#[derive(Derivative)]
#[derivative(Clone)]
pub struct PoseidonHashParametersGadget<
    ConstraintF: PrimeField,
    P: FieldBasedHashParameters<Fr = ConstraintF>
>
{
    params: P,
    _field: PhantomData<ConstraintF>,
}

impl<ConstraintF, P> AllocGadget<P, ConstraintF> for PoseidonHashParametersGadget<ConstraintF, P>
    where
        ConstraintF: PrimeField,
        P: FieldBasedHashParameters<Fr = ConstraintF>,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(_cs: CS, f: F) -> Result<Self, SynthesisError> where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<P>
    {
        let params = f()?.borrow().clone();
        Ok(PoseidonHashParametersGadget {
            params,
            _field: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(_cs: CS, f: F) -> Result<Self, SynthesisError> where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<P>
    {
        let params = f()?.borrow().clone();
        Ok(PoseidonHashParametersGadget {
            params,
            _field: PhantomData,
        })
    }
}

impl<ConstraintF, P> FieldBasedHashGadget<PoseidonHash<ConstraintF, P>, ConstraintF>
for PoseidonHashGadget<ConstraintF, P>
    where
        ConstraintF: PrimeField,
        P:           FieldBasedHashParameters<Fr = ConstraintF>,
{
    type DataGadget = FpGadget<ConstraintF>;
    type ParametersGadget = PoseidonHashParametersGadget<ConstraintF, P>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>
    (
        _cs: CS,
        input: &[Self::DataGadget]
    ) -> Result<Self::DataGadget, SynthesisError>
    {
        Ok(input[0].clone())
    }
}

use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use crate::crh::{MNT4HashParameters, MNT6HashParameters};

pub type MNT4PoseidonHashGadget = PoseidonHashGadget<MNT4753Fr, MNT4HashParameters>;
pub type MNT6PoseidonHashGadget = PoseidonHashGadget<MNT6753Fr, MNT6HashParameters>;