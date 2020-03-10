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
use crate::crh::{PoseidonParameters, PoseidonHash, FieldBasedHash};
use r1cs_std::fields::fp::FpGadget;
use std::marker::PhantomData;

pub trait FieldBasedHashGadget<H: FieldBasedHash<Data = ConstraintF>, ConstraintF: Field>: Sized {
    type DataGadget: FieldGadget<ConstraintF, ConstraintF>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>;
}

pub struct PoseidonHashGadget
<
    ConstraintF: Field,
    P:           PoseidonParameters<Fr = ConstraintF>,
>
{
    _field:      PhantomData<ConstraintF>,
    _parameters: PhantomData<P>,
}

impl<ConstraintF, P> FieldBasedHashGadget<PoseidonHash<ConstraintF, P>, ConstraintF>
for PoseidonHashGadget<ConstraintF, P>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>
    (
        mut cs: CS,
        input: &[Self::DataGadget]
    ) -> Result<Self::DataGadget, SynthesisError>

    {
        //Dummy impl, just for test
        let mut res = Self::DataGadget::zero(cs.ns(|| "alloc result"))?;
        for (i, fg) in input.iter().enumerate() {
            res = res.add(cs.ns(|| format!("add_{}", i)), fg)?;
        }
        Ok(res)
    }
}

use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::bls12_381::Fr as BLS12Fr;
use crate::crh::{MNT4HashParameters, MNT6HashParameters, Bls12_381HashParameters};

pub type MNT4PoseidonHashGadget = PoseidonHashGadget<MNT4753Fr, MNT4HashParameters>;
pub type MNT6PoseidonHashGadget = PoseidonHashGadget<MNT6753Fr, MNT6HashParameters>;
pub type BLS12PoseidonHashGadget = PoseidonHashGadget<BLS12Fr, Bls12_381HashParameters>;