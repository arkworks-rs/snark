extern crate rand;
extern crate rayon;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{PrimeField, MulShort};

use std::marker::PhantomData;

use crate::crh::{
    FieldBasedHashParameters, poseidon::{
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters}
    }
};

use crate::crh::FieldBasedHash;
use crate::Error;

pub mod parameters;
pub mod updatable;
pub mod batched_crh;

pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

// pub struct PoseidonBatchHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
//     _field:      PhantomData<F>,
//     _parameters: PhantomData<P>,
// }

pub trait PoseidonParameters: 'static + FieldBasedHashParameters{

    const T: usize;  // Number of S-Boxes
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const R:usize;   // The rate of the hash function
    const ZERO:Self::Fr;   // The zero element in the field
    const C2:Self::Fr;     // The constant to add in the position corresponding to the capacity
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix
    const MDS_CST_SHORT: &'static[Self::Fr];  // The MDS matrix for fast matrix multiplication

}

// Function that does the scalar multiplication
// It uses Montgomery multiplication
// Constants are defined such that the result is x * t * 2^n mod M,
// that is the Montgomery representation of the operand x * t mod M, and t is the 64-bit constant
#[allow(dead_code)]
#[inline]
pub fn scalar_mul<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {

    state.iter().for_each(|x| {
        let elem = x.mul(&P::MDS_CST[start_idx_cst]);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix
#[allow(dead_code)]
#[inline]
pub fn matrix_mix<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>>  (state: &mut Vec<F>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![F::zero(); P::T];

    let mut idx_cst = 0;
    for i in 0..P::T {
        scalar_mul::<F,P>(&mut new_state[i], state, idx_cst);
        idx_cst += P::T;
    }
    *state = new_state;
}

// Function that does the scalar multiplication
// It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
// t is a 64-bit matrix constant. In the algorithm, the constants are represented in
// partial Montgomery representation, i.e. t * 2^64 mod M
#[inline]
pub fn scalar_mul_fast<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {
    state.iter().for_each(|x| {
        let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(&x);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix with fast algorithm
#[inline]
pub fn matrix_mix_short<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> (state: &mut Vec<F>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![F::zero(); P::T];

    let mut idx_cst = 0;
    for i in 0..P::T {
        scalar_mul_fast::<F,P>(&mut new_state[i], state, idx_cst);
        idx_cst += P::T;
    }
    *state = new_state;
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonHash<F, P> {

    pub(crate) fn poseidon_perm (state: &mut Vec<F>) {


        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Apply Montomgery's simulateneous inversion
            let w2 = state[0] * &state[1];
            let w = state[2] * &w2;
            if w == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                // Calculate inverses individually
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = w.inverse().unwrap();

                let z_2 = w_bar * &w2;
                w_bar = w_bar * &state[2];
                state[2] = z_2;
                let z_1 = w_bar * &state[0];
                state[0] = w_bar * &state[1];
                state[1] = z_1;
            }

            // Perform the matrix mix
            matrix_mix_short::<F,P>(state);

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply S-BOX only to the first element of the state vector
            if state[0]!=P::Fr::zero() {
                state[0] = state[0].inverse().unwrap();
            }

            // Apply the matrix mix
            matrix_mix_short::<F,P>(state);
        }

        // Second full rounds
        // Process only to R_F - 1 iterations. The last iteration does not contain a matrix mix
        for _i in 0..(P::R_F-1) {

            // Add the round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Apply Montgomery's simulatenous inversion
            let w2 = state[0] * &state[1];
            let w = state[2] * &w2;
            if w == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = w.inverse().unwrap();

                let z_2 = w_bar * &w2;
                w_bar = w_bar * &state[2];
                state[2] = z_2;
                let z_1 = w_bar * &state[0];
                state[0] = w_bar * &state[1];
                state[1] = z_1;
            }

            // Apply matrix mix
            matrix_mix_short::<F,P>(state);
        }

        // Last full round does not perform the matrix_mix
        // Add the round constants
        for d in state.iter_mut() {
            let rc = P::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        // Apply Montgomery's simultaneous inversion
        let w2 = state[0] * &state[1];
        let w = state[2] * &w2;
        if w == P::Fr::zero() {
            for d in state.iter_mut() {
                // The S-BOX is an inversion function
                if *d != P::Fr::zero() {
                    *d = (*d).inverse().unwrap();
                }
            }
        } else {
            let mut w_bar = w.inverse().unwrap();

            let z_2 = w_bar * &w2;
            w_bar = w_bar * &state[2];
            state[2] = z_2;
            let z_1 = w_bar * &state[0];
            state[0] = w_bar * &state[1];
            state[1] = z_1;
        }
    }
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>> FieldBasedHash for PoseidonHash<F, P> {
    type Data = F;
    type Parameters = P;

    // Assumption:
    //     capacity c = 1
    fn evaluate(input: &[F]) -> Result<F, Error> {

        assert_ne!(input.len(), 0, "Input data array does not contain any data.");
        assert_eq!(P::T - P::R, 1, "The assumption that the capacity is one field element is not satisfied.");

        // state is a vector of P::T elements.
        // They are initialized to constants that are obtained after applying a permutation to a zero elements vector.
        let mut state = Vec::new();
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }

        // calculate the number of cycles to process the input dividing in portions of rate elements
        let num_cycles = input.len() / P::R;
        // check if the input is a multiple of the rate by calculating the remainder of the division
        let rem = input.len() % P::R;

        // index to process the input
        let mut input_idx = 0;
        // iterate of the portions of rate elements
        for _i in 0..num_cycles {
            // add the elements to the state vector. Add rate elements
            for j in 0..P::R {
                state[j] += &input[input_idx];
                input_idx += 1;
            }
            // add the constant associated to the m-ary Merkle tree
            state[P::R] += &P::C2;

            // apply permutation after adding the input vector
            Self::poseidon_perm(&mut state);
        }

        // in case the input is not a multiple of the rate process the remainder part padding a zero
        if rem != 0 {
            for j in 0..rem {
                state[j] += &input[input_idx];
                input_idx += 1;
            }
            // add the constant associated to the m-ary Merkle tree
            // assumption capacity = 1
            state[P::R] += &P::C2;
            // apply permutation after adding the input vector
            Self::poseidon_perm(&mut state);
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0])
    }
}

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

#[cfg(test)]
mod test {
    use super::*;
    use std::str::FromStr;
    use crate::FieldBasedHash;
    use algebra::biginteger::BigInteger768;

    #[test]
    fn test_poseidon_hash_mnt4() {
        let expected_output = MNT4753Fr::new(BigInteger768([120759599714708995, 15132412086599307425, 1270378153255747692, 3280164418217209635, 5680179791594071572, 2475152338055275001, 9455820118751334058, 6363436228419696186, 3538976751580678769, 14987158621073838958, 10703097083485496843, 48481977539350]));
        let mut input = Vec::new();
        input.push(MNT4753Fr::from_str("1").unwrap());
        input.push(MNT4753Fr::from_str("2").unwrap());
        let output = MNT4PoseidonHash::evaluate(&input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    #[should_panic]
    fn test_poseidon_hash_mnt4_null_element() {
        let input = Vec::new();
        let output = MNT4PoseidonHash::evaluate(&input);
        println!("{:?}", output);
    }

    #[test]
    fn test_poseidon_hash_mnt4_single_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([10133114337753187244, 13011129467758174047, 14520750556687040981, 911508844858788085, 1859877757310385382, 9602832310351473622, 8300303689130833769, 981323167857397563, 5760566649679562093, 8644351468476031499, 10679665778836668809, 404482168782668]));
        let mut input = Vec::new();
        input.push(MNT4753Fr::from_str("1").unwrap());
        let output = MNT4PoseidonHash::evaluate(&input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt4_three_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([5991160601160569512, 9804741598782512164, 8257389273544061943, 15170134696519047397, 9908596892162673198, 7815454566429677811, 9000639780203615183, 8443915450757188195, 1987926952117715938, 17724141978374492147, 13890449093436164383, 191068391234529]));
        let mut input = Vec::new();
        input.push(MNT4753Fr::from_str("1").unwrap());
        input.push(MNT4753Fr::from_str("2").unwrap());
        input.push(MNT4753Fr::from_str("3").unwrap());
        let output = MNT4PoseidonHash::evaluate(&input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([8195238283171732026, 13694263410588344527, 1885103367289967816, 17142467091011072910, 13844754763865913168, 14332001103319040991, 8911700442280604823, 6452872831806760781, 17467681867740706391, 5384727593134901588, 2343350281633109128, 244405261698305]));
        let mut input = Vec::new();
        input.push(MNT6753Fr::from_str("1").unwrap());
        input.push(MNT6753Fr::from_str("2").unwrap());
        let output = MNT6PoseidonHash::evaluate(&mut input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT6753.");
    }

    #[test]
    #[should_panic]
    fn test_poseidon_hash_mnt6_null_element() {
        let input = Vec::new();
        let output = MNT6PoseidonHash::evaluate(&input);
        println!("{:?}", output);
    }

    #[test]
    fn test_poseidon_hash_mnt6_single_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([9820480440897423048, 13953114361017832007, 6124683910518350026, 12198883805142820977, 16542063359667049427, 16554395404701520536, 6092728884107650560, 1511127385771028618, 14755502041894115317, 9806346309586473535, 5880260960930089738, 191119811429922]));
        let mut input = Vec::new();
        input.push(MNT6753Fr::from_str("1").unwrap());
        let output = MNT6PoseidonHash::evaluate(&input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT6753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6_three_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([13800884891843937189, 3814452749758584714, 14612220153016028606, 15886322817426727111, 12444362646204085653, 5214641378156871899, 4248022398370599899, 5982332416470364372, 3842784910369906888, 11445718704595887413, 5723531295320926061, 101830932453997]));
        let mut input = Vec::new();
        input.push(MNT6753Fr::from_str("1").unwrap());
        input.push(MNT6753Fr::from_str("2").unwrap());
        input.push(MNT6753Fr::from_str("3").unwrap());
        let output = MNT6PoseidonHash::evaluate(&input);
        assert_eq!(output.unwrap(), expected_output, "Outputs do not match for MNT6753.");
    }

    #[test]
    fn gen_smt_empty_hash_mnt4() {
        let mut input = Vec::new();
        input.push(MNT4753Fr::from_str("0").unwrap());
        input.push(MNT4753Fr::from_str("0").unwrap());
        let hash_level0 = input.clone();
        println!("level 0: {:?}", hash_level0[0]);

        let level1 = MNT4PoseidonHash::evaluate(&input);
        let mut hash_level = level1.unwrap().clone();
        println!("level 1: {:?}", hash_level);

        for i in 0..40 {
            input.clear();
            let input_level = hash_level.clone();
            input.push(input_level);
            input.push(input_level);
            let next_level = MNT4PoseidonHash::evaluate(&input);
            hash_level = next_level.unwrap().clone();

            let next_hash = hash_level.clone();
            println!("level {}: {:?}", i, next_hash);
        }
    }

    #[test]
    fn gen_smt_empty_hash_mnt6() {
        let mut input = Vec::new();
        input.push(MNT6753Fr::from_str("0").unwrap());
        input.push(MNT6753Fr::from_str("0").unwrap());
        let hash_level0 = input.clone();
        println!("level 0: {:?}", hash_level0[0]);

        let level1 = MNT6PoseidonHash::evaluate(&input);
        let mut hash_level = level1.unwrap().clone();
        println!("level 1: {:?}", hash_level);

        for i in 0..40 {
            input.clear();
            let input_level = hash_level.clone();
            input.push(input_level);
            input.push(input_level);
            let next_level = MNT6PoseidonHash::evaluate(&input);
            hash_level = next_level.unwrap().clone();

            let next_hash = hash_level.clone();
            println!("level {}: {:?}", i, next_hash);
        }
    }

}