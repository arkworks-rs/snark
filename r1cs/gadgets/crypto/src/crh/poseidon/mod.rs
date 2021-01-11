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

#[cfg(feature = "bn_382")]
pub mod bn382;
#[cfg(feature = "bn_382")]
pub use self::bn382::*;
use primitives::PoseidonSBox;
use crate::AlgebraicSpongeGadget;

pub struct PoseidonHashGadget
<
    ConstraintF: PrimeField,
    P:           PoseidonParameters<Fr = ConstraintF>,
    SB:          PoseidonSBox<P>,
    SBG:         SBoxGadget<ConstraintF, SB>,
>
{
    state:   Vec<FpGadget<ConstraintF>>,
    pending: Vec<FpGadget<ConstraintF>>,

    _parameters:        PhantomData<P>,
    _sbox:              PhantomData<SB>,
    _sbox_gadget:       PhantomData<SBG>,
}

impl<
    ConstraintF: PrimeField,
    P:   PoseidonParameters<Fr = ConstraintF>,
    SB:  PoseidonSBox<P>,
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

        {
            // Add initial round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }
        }

        // First full rounds
        for i in 0..P::R_F {

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(||format!("mod_inv_S-Box_1_{}_{}",i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_first_full_round_{}", i)), state)?;

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_1_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Apply S-Box only to the first element of the state vector
            SBG::apply(
                cs.ns(||format!("mod_inv_S-Box_2_{}_{}",_i, 0)),
                &mut state[0]
            )?;

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_partial_round_{}", _i)), state)?;

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_2_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }
        }

        // Second full rounds
        // Process only to R_F -1 iterations. The last iteration does not contain a matrix mix
        for _i in 0..(P::R_F-1) {

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(||format!("mod_inv_S-Box_3_{}_{}",_i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_second_full_round_{}", _i)), state)?;

            // Add the round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_3_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }
        }

        // Last full round does not perform the matrix_mix
        {
            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(|| format!("mod_inv_S-Box_4_{}_{}", P::R_F-1, j)), d)?;
            }
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
            SB:          PoseidonSBox<P>,
            SBG:         SBoxGadget<ConstraintF, SB>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>
    // Assumption:
    //     capacity c = 1
    {
        assert_ne!(input.len(), 0, "Input data array does not contain any data.");
        assert_eq!(P::T - P::R, 1, "The assumption that the capacity is one field element is not satisfied.");

        let mut state = Vec::new();
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(cs.ns(|| format!("hardcode_state_{}",i)), &P::AFTER_ZERO_PERM[i]);
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
            // for application to a 2-1 Merkle tree, add the constant 3 to the third state vector
            state[P::R].add_constant_in_place(cs.ns(|| format!("add_constant_C2_{}", i)), &P::C2)?;
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| format!("poseidon_perm_{}", i)), &mut state)?;
        }

        // in case the input is not a multiple of the rate, process the remainder part padding zeros
        if rem != 0 {
            for j in 0..rem {
                state[j].add_in_place(cs.ns(|| format!("poseidon_padding_add_{}",j)), &input[input_idx])?;
                input_idx += 1;
            }
            // add the constant associated to the m-ary Merkle tree
            // assumption capacity = 1
            state[P::R].add_constant_in_place(cs.ns(|| "add_constant_C2_last_chunk"), &P::C2)?;
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| "poseidon_padding_perm"), &mut state)?;
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0].clone())
    }
}

impl<ConstraintF, P, SB, SBG> AlgebraicSpongeGadget<PoseidonHash<ConstraintF, P, SB>, ConstraintF>
for PoseidonHashGadget<ConstraintF, P, SB, SBG>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
        SB:          PoseidonSBox<P>,
        SBG:         SBoxGadget<ConstraintF, SB>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn new<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        assert_eq!(P::T - P::R, 1, "The assumption that the capacity is one field element is not satisfied.");

        let mut state = Vec::new();
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(
                cs.ns(|| format!("hardcode_state_{}",i)),
                &P::AFTER_ZERO_PERM[i]
            );
            state.push(elem);
        }

        Ok(Self {
            state,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
            _sbox: PhantomData,
            _sbox_gadget: PhantomData,
        })
    }

    fn enforce_absorb<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        elems: &[Self::DataGadget]
    ) -> Result<(), SynthesisError> {
        elems.iter().enumerate().map(|(i, f)| {
            self.pending.push(f.clone());
            if self.pending.len() == P::R {
                // add the elements to the state vector. Add rate elements
                for (j, (input, state)) in self.pending.iter().zip(self.state.iter_mut()).enumerate() {
                    state.add_in_place(cs.ns(|| format!("add_input_{}_{}_to_state", i, j)), input)?;
                }
                // apply permutation after adding the input vector
                Self::poseidon_perm(cs.ns(|| format!("poseidon_perm_{}", i)), &mut self.state)?;

                self.pending.clear();
            }
            Ok(())
        }).collect::<Result<(), SynthesisError>>()?;
        Ok(())
    }

    fn enforce_squeeze<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        num: usize
    ) -> Result<Vec<Self::DataGadget>, SynthesisError> {
        let mut state = self.state.clone();
        for (i, (input, state)) in self.pending.iter().zip(state.iter_mut()).enumerate() {
            state.add_in_place(cs.ns(|| format!("add_input_{}_to_state", i)), input)?;
        }
        let mut output = Vec::with_capacity(num);
        for i in 0..num {
            Self::poseidon_perm(cs.ns(|| format!("squeeze field {}", i)), &mut state)?;
            output.push(state[0].clone())
        }
        Ok(output)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        MNT4PoseidonHashGadget, MNT6PoseidonHashGadget,
        BN382FqPoseidonHashGadget, BN382FrPoseidonHashGadget,
    };

    use algebra::fields::PrimeField;
    use crate::crh::test::{
        field_based_hash_gadget_native_test, algebraic_sponge_gadget_native_test
    };

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    #[test]
    fn poseidon_mnt4_753_gadget_native_test() {
        field_based_hash_gadget_native_test::<_, _, MNT4PoseidonHashGadget>(generate_inputs(2));
        algebraic_sponge_gadget_native_test::<_, _, MNT4PoseidonHashGadget>(generate_inputs(5));
    }

    #[test]
    fn poseidon_mnt6_753_gadget_native_test() {
        field_based_hash_gadget_native_test::<_, _, MNT6PoseidonHashGadget>(generate_inputs(2));
        algebraic_sponge_gadget_native_test::<_, _, MNT6PoseidonHashGadget>(generate_inputs(5));
    }

    #[test]
    fn poseidon_bn382_fr_gadget_native_test() {
        field_based_hash_gadget_native_test::<_, _, BN382FrPoseidonHashGadget>(generate_inputs(2));
        algebraic_sponge_gadget_native_test::<_, _, BN382FrPoseidonHashGadget>(generate_inputs(5));
    }

    #[test]
    fn poseidon_bn382_fq_gadget_native_test() {
        field_based_hash_gadget_native_test::<_, _, BN382FqPoseidonHashGadget>(generate_inputs(2));
        algebraic_sponge_gadget_native_test::<_, _, BN382FqPoseidonHashGadget>(generate_inputs(5));
    }
}