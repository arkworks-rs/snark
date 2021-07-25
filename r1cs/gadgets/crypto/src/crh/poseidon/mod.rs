use algebra::PrimeField;
use primitives::crh::poseidon::{
        PoseidonHash, PoseidonParameters
};
use crate::crh::{
    SBoxGadget, FieldBasedHashGadget
};
use r1cs_std::{
    fields::{
        FieldGadget, fp::FpGadget
    },
    alloc::ConstantGadget,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::marker::PhantomData;

#[cfg(feature = "mnt4_753")]
pub mod mnt4753;
#[cfg(feature = "mnt4_753")]
pub use self::mnt4753::*;

#[cfg(feature = "mnt6_753")]
pub mod mnt6753;
#[cfg(feature = "mnt6_753")]
pub use self::mnt6753::*;

#[cfg(feature = "tweedle")]
pub mod tweedle;
#[cfg(feature = "tweedle")]
pub use self::tweedle::*;

#[cfg(feature = "bn_382")]
pub mod bn382;
#[cfg(feature = "bn_382")]
pub use self::bn382::*;

use primitives::SBox;

pub struct PoseidonHashGadget
<
    ConstraintF: PrimeField,
    P:           PoseidonParameters<Fr = ConstraintF>,
    SB:          SBox<Field = ConstraintF, Parameters = P>,
    SBG:         SBoxGadget<ConstraintF, SB>,
>
{
    _field:             PhantomData<ConstraintF>,
    _parameters:        PhantomData<P>,
    _sbox:              PhantomData<SB>,
    _sbox_gadget:       PhantomData<SBG>,
}

impl<
    ConstraintF: PrimeField,
    P:   PoseidonParameters<Fr = ConstraintF>,
    SB:  SBox<Field = ConstraintF, Parameters = P>,
    SBG: SBoxGadget<ConstraintF, SB>
> PoseidonHashGadget<ConstraintF, P, SB, SBG>
{

    fn poseidon_perm<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        state: &mut [FpGadget<ConstraintF>],
    ) -> Result<(), SynthesisError>
    {
        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(||format!("S-Box_1_{}_{}",i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_first_full_round_{}", i)), state)?;

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply S-Box only to the first element of the state vector
            SBG::apply(
                cs.ns(||format!("S-Box_2_{}_{}",_i, 0)),
                &mut state[0]
            )?;

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_partial_round_{}", _i)), state)?;
        }

        // Second full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(|| format!("S-Box_3_{}_{}", _i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix(cs.ns(|| format!("poseidon_mix_matrix_second_full_round_{}", _i)), state)?;
        }
        Ok(())
    }

    // Function that does the dot product for the mix matrix
    fn dot_prod<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        res: &mut FpGadget<ConstraintF>,
        state: &mut [FpGadget<ConstraintF>],
        mut start_idx_cst: usize,
    ) -> Result<(), SynthesisError>
    {
        for x in state.iter() {
            let elem = x.mul_by_constant(cs.ns(|| format!("partial_product_{}", start_idx_cst)), &P::MDS_CST[start_idx_cst])?;
            start_idx_cst += 1;
            (*res).add_in_place(cs.ns(|| format!("add_partial_product_{}", start_idx_cst)), &elem)?;
        }

        Ok(())
    }

    // Function that does the mix matrix
    fn matrix_mix<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        state: &mut [FpGadget<ConstraintF>],
    ) -> Result<(), SynthesisError>
    {

        // Check that the length of the state vector is t
        assert_eq!(state.len(), P::T);

        // Destination state vector
        let mut new_state = Vec::new();

        // Initialize new destination state vector with zero elements
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(cs.ns(|| format!("hardcode_new_state_elem_{}", i)), &P::ZERO);
            new_state.push(elem);
        }

        // Performs the dot products
        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::dot_prod(cs.ns(|| format!("poseidon_dot_product_{}", i)), &mut new_state[i], state, idx_cst)?;
            idx_cst += P::T;
        }

        // Copy result to the state vector
        for i in 0..P::T {
            state[i] = new_state[i].clone();
        }

        Ok(())
    }
}

impl<ConstraintF, P, SB, SBG> FieldBasedHashGadget<PoseidonHash<ConstraintF, P, SB>, ConstraintF>
    for PoseidonHashGadget<ConstraintF, P, SB, SBG>
        where
            ConstraintF: PrimeField,
            P:           PoseidonParameters<Fr = ConstraintF>,
            SB:          SBox<Field = ConstraintF, Parameters = P>,
            SBG:         SBoxGadget<ConstraintF, SB>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn enforce_hash_constant_length<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>
    // Assumption:
    //     capacity c = 1
    {
        if input.len() == 0 {
            Err(SynthesisError::Other("Input data array does not contain any data".to_owned()))?
        }

        let mut state = Vec::new();
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(
                cs.ns(|| format!("hardcode_state_{}",i)),
                &P::AFTER_ZERO_PERM[i]
            );
            state.push(elem);
        }

        // calculate the number of cycles to process the input dividing in portions of rate elements
        let num_cycles = input.len() / P::R;
        // check if the input is a multiple of the rate by calculating the remainder of the division
        // the remainder of dividing the input length by the rate can be 1 or 0 because we are assuming
        // that the rate is 2
        let rem = input.len() % P::R;

        // index to process the input
        let mut input_idx = 0;
        // iterate of the portions of rate elements
        for i in 0..num_cycles {
            // add the elements to the state vector. Add rate elements
            for j in 0..P::R {
                state[j].add_in_place(cs.ns(|| format!("add_input_{}_{}", i, j)), &input[input_idx])?;
                input_idx += 1;
            }
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| format!("poseidon_perm_{}", i)), &mut state)?;
        }

        // in case the input is not a multiple of the rate, process the remainder part padding zeros
        if rem != 0 {
            for j in 0..rem {
                state[j].add_in_place(cs.ns(|| format!("poseidon_padding_add_{}",j)), &input[input_idx])?;
                input_idx += 1;
            }
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| "poseidon_padding_perm"), &mut state)?;
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0].clone())
    }
}

#[cfg(test)]
mod test {

    use algebra::PrimeField;
    use crate::crh::test::constant_length_field_based_hash_gadget_native_test;

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    #[cfg(feature = "mnt4_753")]
    #[test]
    fn poseidon_mnt4_753_gadget_native_test() {
        use crate::MNT4PoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, MNT4PoseidonHashGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "mnt6_753")]
    #[test]
    fn poseidon_mnt6_753_gadget_native_test() {
        use crate::MNT6PoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, MNT6PoseidonHashGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn crh_bn382_fr_primitive_gadget_test() {
        use crate::BN382FrPoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, BN382FrPoseidonHashGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn crh_bn382_fq_primitive_gadget_test() {
        use crate::BN382FqPoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, BN382FqPoseidonHashGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn crh_tweedle_fr_primitive_gadget_test() {
        use crate::TweedleFrPoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, TweedleFrPoseidonHashGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn crh_tweedle_fq_primitive_gadget_test() {
        use crate::TweedleFqPoseidonHashGadget;

        for ins in 1..=3 {
            constant_length_field_based_hash_gadget_native_test::<_, _, TweedleFqPoseidonHashGadget>(generate_inputs(ins));
        }
    }
}