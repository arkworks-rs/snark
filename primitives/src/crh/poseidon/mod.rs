extern crate rand;
extern crate rayon;

use algebra::{PrimeField, MulShort};

use std::marker::PhantomData;

use crate::{
    crh::{
        FieldBasedHash,
        FieldBasedHashParameters,
    }, CryptoError, Error
};

pub mod batched_crh;

pub mod parameters;
pub use self::parameters::*;

#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
)]
pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    state: Vec<F>,
    pending: Vec<F>,
    input_size: Option<usize>,
    updates_ctr: usize,
    mod_rate: bool,
    _parameters: PhantomData<P>,
}

pub trait PoseidonParameters: 'static + FieldBasedHashParameters + Clone {

    const T: usize;  // Number of S-Boxes
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const ZERO:Self::Fr;   // The zero element in the field
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
pub fn scalar_mul<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {

    state.iter().for_each(|x| {
        let elem = x.mul(&P::MDS_CST[start_idx_cst]);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix
#[allow(dead_code)]
#[inline]
pub fn matrix_mix<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>>  (state: &mut Vec<F>) {

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
pub fn scalar_mul_fast<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {
    state.iter().for_each(|&x| {
        let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(x);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix with fast algorithm
#[inline]
pub fn matrix_mix_short<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (state: &mut Vec<F>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![F::zero(); P::T];

    let mut idx_cst = 0;
    for i in 0..P::T {
        scalar_mul_fast::<F,P>(&mut new_state[i], state, idx_cst);
        idx_cst += P::T;
    }
    *state = new_state;
}

impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> PoseidonHash<F, P> {

    fn _init(constant_size: Option<usize>, mod_rate: bool, personalization: Option<&[F]>) -> Self
    {
        assert_eq!(P::T - P::R, 1, "The assumption that the capacity is one field element is not satisfied.");

        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }


        let mut instance = Self {
            state,
            pending: Vec::with_capacity(P::R),
            input_size: constant_size,
            updates_ctr: 0,
            mod_rate,
            _parameters: PhantomData,
        };

        // If personalization Vec is not multiple of the rate, we pad it with zero field elements.
        // This will allow eventually to precompute the constants of the initial state. This
        // is exactly as doing H(personalization, padding, ...). NOTE: this way of personalizing
        // the hash is not mentioned in https://eprint.iacr.org/2019/458.pdf
        if personalization.is_some(){
            let personalization = personalization.unwrap();

            for &p in personalization.into_iter(){
                instance.update(p);
            }

            let padding = if personalization.len() % P::R != 0 {
                P::R - ( personalization.len() % P::R )
            } else {
                0
            };

            for _ in 0..padding {
                instance.update(F::zero());
            }
            assert_eq!(instance.pending.len(), 0);
            instance.updates_ctr = 0;
        }
        instance
    }

    #[inline]
    fn apply_permutation_if_needed(&mut self) {
        if self.pending.len() == P::R {
            for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
                *state += input;
            }
            Self::poseidon_perm(&mut self.state);
            self.pending.clear();
        }
    }

    #[inline]
    fn get_hash(mut state: Vec<F>, inputs: Vec<F>) -> F {
        for (input, s) in inputs.iter().zip(state.iter_mut()) {
            *s += input;
        }
        Self::poseidon_perm(&mut state);
        state[0]
    }

    #[inline]
    // Padding strategy is described in https://eprint.iacr.org/2019/458.pdf(Section 4.2)
    fn pad_and_finalize(&self) -> F
    {
        // Constant input length instance
        if self.input_size.is_some() {

            // The constant size is modulus rate, so we already have the hash output in state[0]
            // as a permutation is applied each time pending reaches P::R length
            if self.pending.is_empty() {
                self.state[0].clone()
            }

            // Pending is not empty: pad with 0s up to rate then compute the hash
            else {
                let mut pending = self.pending.clone();
                pending.append(&mut vec![F::zero(); P::R - (pending.len() % P::R)]);
                Self::get_hash(self.state.clone(), pending)
            }
        }

        // Variable input length instance
        else {

            // The input is of variable length, but always modulus rate: result is already available
            // in state[0] as a permutation is applied each time pending reaches P::R length
            if self.mod_rate {
                self.state[0].clone()
            }

            // The input is of variable length, but not modulus rate: we always need to apply
            // padding. Pad with a single 1 and then 0s up to rate. Compute hash.
            else {
                let mut pending = self.pending.clone();
                pending.push(F::one());
                pending.append(&mut vec![F::zero(); P::R - (pending.len() % P::R)]);
                Self::get_hash(self.state.clone(), pending)
            }
        }
    }

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
            // Use Montgomery simultaneous inversion
            let mut w: Vec<P::Fr> = Vec::new();
            let mut accum_prod = P::Fr::one();

            w.push(accum_prod);

            // Calculate the intermediate partial products
            for d in state.iter() {
                accum_prod = accum_prod * &d;
                w.push(accum_prod);
            }

            if accum_prod == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                // Calculate inverses individually
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = accum_prod.inverse().unwrap();

                // Extract the individual inversions
                let mut idx: i64 = w.len() as i64 - P::R as i64;
                for d in state.iter_mut().rev() {
                    let tmp = d.clone();
                    *d = w_bar * &w[idx as usize];
                    w_bar = w_bar * &tmp;
                    idx -= 1;
                }
            }

            // Perform the matrix mix
            matrix_mix_short::<F, P>(state);

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
            if state[0]!= P::Fr::zero() {
                state[0] = state[0].inverse().unwrap();
            }

            // Apply the matrix mix
            matrix_mix_short::<F, P>(state);
        }

        // Second full rounds
        // Process only to R_F - 1 iterations.
        for _i in 0..(P::R_F - 1) {

            // Add the round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Use Montgomery simultaneous inversion
            let mut w: Vec<P::Fr> = Vec::new();
            let mut accum_prod = P::Fr::one();

            w.push(accum_prod);

            // Calculate the intermediate partial products
            for d in state.iter() {
                accum_prod = accum_prod * &d;
                w.push(accum_prod);
            }

            if accum_prod == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                // Calculate inverses individually
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = accum_prod.inverse().unwrap();

                // Extract the individual inversions
                let mut idx: i64 = w.len() as i64 - P::R as i64;
                for d in state.iter_mut().rev() {
                    let tmp = d.clone();
                    *d = w_bar * &w[idx as usize];
                    w_bar = w_bar * &tmp;
                    idx -= 1;
                }
            }

            // Perform the matrix mix
            matrix_mix_short::<F, P>(state);
        }

        // Last full round does not perform the matrix_mix
        // Add the round constants
        for d in state.iter_mut() {
            let rc = P::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        // Use Montgomery simultaneous inversion
        let mut w: Vec<P::Fr> = Vec::new();
        let mut accum_prod = P::Fr::one();

        w.push(accum_prod);

        // Calculate the intermediate partial products
        for d in state.iter() {
            accum_prod = accum_prod * &d;
            w.push(accum_prod);
        }

        if accum_prod == P::Fr::zero() {
            // At least one of the S-Boxes is zero
            // Calculate inverses individually
            for d in state.iter_mut() {
                // The S-BOX is an inversion function
                if *d != P::Fr::zero() {
                    *d = (*d).inverse().unwrap();
                }
            }
        } else {
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - P::R as i64;
            for d in state.iter_mut().rev() {
                let tmp = d.clone();
                *d = w_bar * &w[idx as usize];
                w_bar = w_bar * &tmp;
                idx -= 1;
            }
        }
    }
}

impl<F, P> FieldBasedHash for PoseidonHash<F, P>
    where
        F: PrimeField + MulShort<F, Output = F>,
        P: PoseidonParameters<Fr = F>,
{
    type Data = F;
    type Parameters = P;

    fn init_constant_length(input_size: usize, personalization: Option<&[Self::Data]>) -> Self {
        Self::_init(
            Some(input_size),
            input_size % P::R == 0, // Not taken into account, can be any
            personalization
        )
    }

    fn init_variable_length(mod_rate: bool, personalization: Option<&[Self::Data]>) -> Self {
        Self::_init(None, mod_rate, personalization)
    }

    fn update(&mut self, input: Self::Data) -> &mut Self {
        self.pending.push(input);
        self.updates_ctr += 1;
        self.apply_permutation_if_needed();
        self
    }

    fn finalize(&self) -> Result<Self::Data, Error> {

        let error_condition =

            // Constant input length instance, but the size of the input is different from the declared one
            (self.input_size.is_some() && self.updates_ctr != self.input_size.unwrap())

            ||

            // Variable modulus rate input length instance, but the size of the input is not modulus rate
            (self.input_size.is_none() && self.mod_rate && self.updates_ctr % P::R != 0);

        // If one of the conditions above is true, we must throw an error
        if error_condition
        {
            Err(Box::new(CryptoError::HashingError("attempt to finalize with an input of invalid size".to_owned())))
        }

        // Otherwise pad if needed (according to the Self instance type) and return the hash output
        else
        {
            Ok(self.pad_and_finalize())
        }
    }

    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self {
        let new_instance = Self::_init(self.input_size, self.mod_rate, personalization);
        *self = new_instance;
        self
    }
}

#[cfg(test)]
mod test {
    use algebra::{Field, PrimeField};
    use crate::crh::{
        FieldBasedHash,
        test::{
            field_based_hash_regression_test, constant_length_field_based_hash_test, variable_length_field_based_hash_test
        }
    };
    use crate::FieldBasedHashParameters;

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    fn test_routine<F: PrimeField, H: FieldBasedHash<Data = F>>(
        expected_outputs: Vec<(F, F)>
    )
    {
        let rate = <H::Parameters as FieldBasedHashParameters>::R;
        expected_outputs.into_iter().enumerate().for_each(|(i, (expected_output_constant, expected_output_variable))| {

            let ins = generate_inputs::<F>(i + 1);

            // Constant length
            {
                let mut digest = H::init_constant_length(i + 1, None);

                /*field_based_hash_regression_test::<H>(
                    &mut digest,
                    ins.clone(),
                    expected_output_constant,
                );*/

                constant_length_field_based_hash_test::<H>(
                    &mut digest,
                    ins.clone()
                );
            }

            // Variable length
            {
                let mod_rate = (i + 1) % rate == 0;
                let mut digest = H::init_variable_length(mod_rate, None);

                /*field_based_hash_regression_test::<H>(
                    &mut digest,
                    ins.clone(),
                    expected_output_variable
                );*/

                variable_length_field_based_hash_test::<H>(
                    &mut digest,
                    ins,
                    mod_rate
                );
            }
        });
    }

    #[cfg(feature = "mnt4_753")]
    #[test]
    fn test_poseidon_hash_mnt4() {
        use algebra::{
            biginteger::BigInteger768,
            fields::mnt4753::Fr as MNT4753Fr
        };
        use crate::crh::poseidon::parameters::mnt4753::MNT4PoseidonHash;

        let expected_outputs = vec![
            (
                MNT4753Fr::zero(),
                MNT4753Fr::zero(),
            ),
            (
                MNT4753Fr::zero(),
                MNT4753Fr::zero(),
            ),
            (
                MNT4753Fr::zero(),
                MNT4753Fr::zero(),
            ),
        ];

        test_routine::<MNT4753Fr, MNT4PoseidonHash>(expected_outputs)
    }

    #[cfg(feature = "mnt6_753")]
    #[test]
    fn test_poseidon_hash_mnt6() {
        use algebra::{
            biginteger::BigInteger768,
            fields::mnt6753::Fr as MNT6753Fr
        };
        use crate::crh::poseidon::parameters::mnt6753::MNT6PoseidonHash;

        let expected_outputs = vec![
            (
                MNT6753Fr::zero(),
                MNT6753Fr::zero(),
            ),
            (
                MNT6753Fr::zero(),
                MNT6753Fr::zero(),
            ),
            (
                MNT6753Fr::zero(),
                MNT6753Fr::zero(),
            ),
        ];

        test_routine::<MNT6753Fr, MNT6PoseidonHash>(expected_outputs)
    }
}