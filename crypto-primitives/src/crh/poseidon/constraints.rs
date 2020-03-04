use algebra::{Field, PrimeField};
use crate::crh::poseidon::PoseidonParameters;
use std::marker::PhantomData;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::fields::fp::FpGadget;
use r1cs_std::fields::FieldGadget;
use r1cs_std::bits::boolean::Boolean;
use r1cs_std::alloc::AllocGadget;
use r1cs_std::Assignment;
use r1cs_std::eq::ConditionalEqGadget;

pub struct PoseidonHashGadget
<
    ConstraintF: PrimeField,
    P:           PoseidonParameters<Fr = ConstraintF>,
>
{
    _field:      PhantomData<ConstraintF>,
    _parameters: PhantomData<P>,
}

impl<ConstraintF, P> PoseidonHashGadget<ConstraintF, P>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>
{
    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        input: &[FpGadget<ConstraintF>],
    ) -> Result<FpGadget<ConstraintF>, SynthesisError>
    // Assumption: The input vector has always a multiple of rate elements
    {

        let zero = FpGadget::<ConstraintF>::zero()?;
        // state is a vector of 3 elements. They are initialized to zero elements
        let mut state = vec![zero, zero.clone(), zero.clone()];

        // calculate the number of cycles to process the input dividing in portions of rate elements
        let num_cycles = input.len() / P::R;
        // check if the input is a multiple of the rate by calculating the remainder of the division
        let rem = input.len() % P::R;

        // apply permutation to all zeros state vector
        Self::poseidon_perm(cs.ns(|| "initial_poseidon_perm"),  &mut state)?;

        // index to process the input
        let mut input_idx = 0;
        // iterate of the portions of rate elements
        for i in 0..num_cycles {
            // add the elements to the state vector. Add rate elements
            for j in 0..P::R {
                state[j].add_in_place(cs.ns(|| format!("add_{}_{}", i, j)), &input[input_idx])?;
                input_idx += 1;
            }
            // for application to a 2-1 Merkle tree, add the constant 3 to the third state vector
            state[P::R].add_constant_in_place(&P::C2)?;
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| format!("poseidon_perm_{}", i)),&mut state)?;
        }

        // in case the input is not a multiple of the rate process the remainder part padding a zero
        // rem_input   0       C2
        // state[0] state[1] state[2]
        if rem != 0 {
            state[0].add_in_place(cs.ns(|| "poseidon_padding_add"), &input[input_idx])?;
            state[P::R].add_constant_in_place(&P::C2)?;
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| "poseidon_padding_perm"), &mut state)?;
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0])
    }

    fn mod_inv_sbox<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        x: &mut FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError>
    {
        let b = Boolean::alloc(cs.ns(|| "alloc b"), || {
            let x_val = x.get_value().get()?;
            if x_val == F::zero() {
                Ok(false)
            } else {
                Ok(true)
            }
        })?;
        let y = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc y"), || {
            let x_val = x.get_value().get()?;
            if x_val == ConstraintF::zero() {
                Ok(x_val)
            } else {
                let inv = x_val.inverse().get()?;
                Ok(inv)
            }
        })?;
        cs.enforce(
            || "b=x*y",
            |lc| &x.variable + lc,
            |lc| &y.variable + lc,
            |lc| lc + &b.lc(CS::one(), ConstraintF::one()),
        );
        x.conditional_enforce_equal(cs.ns(|| "0=(1-b)*(x-y)"), &y, &b.not());

        *x=y;
        Ok(())
    }

    fn poseidon_perm<CS: ConstraintSystem<ConstraintF>>(
         mut cs: CS,
         state: &mut [FpGadget<ConstraintF>],
        ) -> Result<(), SynthesisError>
    {

        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                Self::mod_inv_sbox(cs.ns(||format!("mod_inv_S-Box_{}_{}",_i, j)), d)?;h
            }

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_first_full_round_{}", _i)), &mut state)?;
        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }

            //==============================================================
            //==============================================================
            // Apply S-BOX only to the first element of the state vector
            let t1 = (state[0]).inverse();
            match t1 {
                None => println!("Field inversion error"),
                Some(inv) => {
                    state[0] = inv;
                },
            }
            //==============================================================
            //==============================================================

            // Apply the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_partial_round_{}", _i)), &mut state)?;
        }

        // Second full rounds
        // Process only to R_F -1 iterations. The last iteration does not contain a matrix mix
        for _i in 0..(P::R_F-1) {

            // Add the round constants
            for d in state.iter_mut() {
                //let rc = MNT4753Fr::from_str(ROUND_CST[round_cst_idx]).map_err(|_| ()).unwrap();
                let rc = P::ROUND_CST[round_cst_idx];
                (*d).add_constant_in_place(cs.ns(|| format!("add_constant_{}", round_cst_idx)), &rc)?;
                round_cst_idx += 1;
            }

            //==============================================================
            //==============================================================
            // Apply the S-BOX to each element of the state vector
            for d in state.iter_mut() {
                let elem_state = (*d).inverse();
                match elem_state {
                    None => println!("Field inversion error"),
                    Some(inv) => {
                        *d = inv;
                    },
                }
            }
            //==============================================================
            //==============================================================


            // Apply matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_second_full_round_{}", _i)), &mut state)?;
        }

        // Last full round does not perform the matrix_mix
        // Add the round constants
        for d in state.iter_mut() {
            let rc = P::ROUND_CST[round_cst_idx];
            (*d).add_constant_in_place(cs.ns(|| format!("add_constant_{}", round_cst_idx)), &rc)?;
            round_cst_idx += 1;
        }

        //==============================================================
        //==============================================================
        // Apply the S-BOX to each element of the state vector
        for d in state.iter_mut() {
            let elem_state = (*d).inverse();
            match elem_state {
                None => println!("Field inversion error"),
                Some(inv) => {
                    *d = inv;
                },
            }
        }
        //==============================================================
        //==============================================================

        Ok(())
    }

    // Function that does the mix matrix
    fn matrix_mix<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        state: &mut [FpGadget<ConstraintF>],
    ) -> Result<(), SynthesisError>
    {
    // fn matrix_mix<P:PoseidonParameters> (state: &mut Vec<P::Fr>) {

        // the new state where the result will be stored initialized to zero elements
//        let mut new_state = vec![P::Fr::zero(); P::T];

        //  (s1, s2, s3) * [[ m11 m12 m13]
                            [ m21 m22 m23]
                            [ m31 m32 m33]]

        s1' = s1 * m11 + s2 * m21 + s3 * m3
        s2' = s1 *


        let zero = FpGadget::<ConstraintF>::zero()?;
        // state is a vector of 3 elements. They are initialized to zero elements
        let mut new_state = vec![zero, zero.clone(), zero.clone()];


        let m_11 = P::MDS_CST[0];
        let m_12 = P::MDS_CST[1];
        let m_13 = P::MDS_CST[2];

        // scalar multiplication for position 0 of the state vector
        let elem_0 = state[0].mul(&m_11);
        let elem_1 = state[1].mul(&m_12);
        let elem_2 = state[2].mul(&m_13);

        new_state[0] = elem_0;
        new_state[0] += &elem_1;
        new_state[0] += &elem_2;

        // scalar multiplication for position 1 of the state vector
        let m_21 = P::MDS_CST[3];
        let m_22 = P::MDS_CST[4];
        let m_23 = P::MDS_CST[5];

        let elem_3 = state[0].mul(&m_21);
        let elem_4 = state[1].mul(&m_22);
        let elem_5 = state[2].mul(&m_23);

        new_state[1] = elem_3;
        new_state[1] += &elem_4;
        new_state[1] += &elem_5;

        // scalar multiplication for the position 2 of the state vector
        let m_31 = P::MDS_CST[6];
        let m_32 = P::MDS_CST[7];
        let m_33 = P::MDS_CST[8];

        let elem_6 = state[0].mul(&m_31);
        let elem_7 = state[1].mul(&m_32);
        let elem_8 = state[2].mul(&m_33);

        new_state[2] = elem_6;
        new_state[2] += &elem_7;
        new_state[2] += &elem_8;

        // copy the result to the state vector
        state[0] = new_state[0];
        state[1] = new_state[1];
        state[2] = new_state[2];

        Ok(())
    }

}