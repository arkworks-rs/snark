use algebra::{Field, PrimeField, MulShort};
use crate::{PoseidonParameters, SBox, BatchSBox};
use std::marker::PhantomData;

pub trait PoseidonInverseParameters: PoseidonParameters {
    const MDS_CST_SHORT: &'static [Self::Fr];  // The MDS matrix for fast matrix multiplication
}

pub trait PoseidonSBox<P: PoseidonParameters>: SBox<Field = P::Fr, Parameters = P> {

    // Function that does the scalar multiplication
    fn scalar_mul(res: &mut <Self as SBox>::Field, state: &mut [<Self as SBox>::Field], start_idx_cst: usize);

    // Function that does the mix matrix
    fn matrix_mix(state: &mut Vec<<Self as SBox>::Field>)
    {
        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![<Self as SBox>::Field::zero(); P::T];

        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::scalar_mul(&mut new_state[i], state, idx_cst);
            idx_cst += P::T;
        }
        *state = new_state;
    }
}

pub trait PoseidonBatchSBox<P: PoseidonParameters>: PoseidonSBox<P> + BatchSBox {}

/// S-Box: S(x) = x^-1
#[derive(Debug)]
pub struct PoseidonInverseSBox<F: PrimeField + MulShort<F, Output = F>, P: PoseidonInverseParameters<Fr = F>> {
    _field: PhantomData<F>,
    _parameters: PhantomData<P>,
}

impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonInverseParameters<Fr = F>> SBox for PoseidonInverseSBox<F, P> {
    type Field = F;
    type Parameters = P;

    #[inline]
    fn apply_full(state: &mut Vec<F>, last: bool) {
        // Apply Montomgery's simulateneous inversion
        let w2 = state[0] * &state[1];
        let w = state[2] * &w2;
        if w == F::zero() {
            // At least one of the S-Boxes is zero
            // Calculate inverses individually
            for d in state.iter_mut() {
                // The S-BOX is an inversion function
                if *d != F::zero() {
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

        if !last {
            // Perform the matrix mix
            Self::matrix_mix(state);
        }
    }

    #[inline]
    fn apply_partial(state: &mut Vec<F>) {
        if state[0]!= F::zero() {
            state[0] = state[0].inverse().unwrap();
        }

        // Apply the matrix mix
        Self::matrix_mix(state);
    }
}

impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonInverseParameters<Fr = F>> PoseidonSBox<P> for PoseidonInverseSBox<F, P> {

    // It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
    // t is a 64-bit matrix constant. In the algorithm, the constants are represented in
    // partial Montgomery representation, i.e. t * 2^64 mod M
    #[inline]
    fn scalar_mul(res: &mut F, state: &mut [F], mut start_idx_cst: usize) {
        state.iter().for_each(|&x| {
            let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(x);
            start_idx_cst += 1;
            *res += &elem;
        });
    }
}


impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonInverseParameters<Fr = F>> BatchSBox for PoseidonInverseSBox<F, P> {

    fn apply_full_batch(vec_state: &mut [Vec<F>], last: bool) {

        // Use Montgomery simultaneous inversion
        let mut w: Vec<F> = Vec::new();
        let mut accum_prod = F::one();

        w.push(accum_prod);
        // Calculate the intermediate partial products
        for i in 0..vec_state.len() {
            for j in 0..P::T {
                accum_prod = accum_prod * &vec_state[i][j];
                w.push(accum_prod);
            }
        }

        // if the accum_prod is zero, it means that one of the S-Boxes is zero
        // in that case compute the inverses individually
        if accum_prod == F::zero() {
            for i in 0..vec_state.len() {
                for j in 0..P::T {
                    if vec_state[i][j] != F::zero() {
                        vec_state[i][j] = vec_state[i][j].inverse().unwrap();
                    }
                }
            }
        } else {

            // Calculate the inversion of the products
            // The inverse always exists in this case
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - 2;
            for i in (0..vec_state.len()).rev() {
                for j in (0..P::T).rev() {
                    let vec_1 = vec_state[i][j].clone();
                    vec_state[i][j] = w_bar * &w[idx as usize];
                    w_bar = w_bar * &vec_1;
                    idx -= 1;
                }
            }
        }

        if !last {
            // Perform the matrix mix
            for i in 0..vec_state.len() {
                Self::matrix_mix(&mut vec_state[i]);
            }
        }
    }

    fn apply_partial_batch(vec_state: &mut [Vec<F>]) {
        let mut w: Vec<F> = Vec::new();
        let mut accum_prod = F::one();

        w.push(accum_prod);
        // Calculate the intermediate partial products
        for i in 0..vec_state.len() {
            accum_prod = accum_prod * &vec_state[i][0];
            w.push(accum_prod);
        }

        // if the accum_prod is zero, it means that one of the S-Boxes is zero
        // in that case compute the inverses individually
        if accum_prod == F::zero() {
            for i in 0..(vec_state.len() - 1) {
                if vec_state[i][0] != F::zero() {
                    vec_state[i][0] = vec_state[i][0].inverse().unwrap();
                }
            }
        } else {

            // Calculate the inversion of the products
            // Use Montgomery simultaneous inversion
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - 2;
            for i in (0..vec_state.len()).rev() {
                let vec_1 = vec_state[i][0].clone();
                vec_state[i][0] = w_bar * &w[idx as usize];
                w_bar = w_bar * &vec_1;
                idx -= 1;
            }
        }

        // Perform the matrix mix
        for i in 0..vec_state.len() {
            Self::matrix_mix(&mut vec_state[i]);
        }
    }
}

impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonInverseParameters<Fr = F>> PoseidonBatchSBox<P>
for PoseidonInverseSBox<F, P> {}

/// S-Box: S(x) = x^5
#[derive(Debug)]
pub struct PoseidonQuinticSBox<F: PrimeField, P: PoseidonParameters<Fr = F>> {
    _field: PhantomData<F>,
    _parameters: PhantomData<P>,
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> PoseidonQuinticSBox<F, P> {
    #[inline]
    fn exp_by_5(elem: &mut F) {
        let w1 = *elem * &(*elem);
        let w2 = w1 * &w1;
        *elem = w2 * &(*elem);
    }
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> PoseidonSBox<P> for PoseidonQuinticSBox<F, P> {
    // It uses Montgomery multiplication
    // Constants are defined such that the result is x * t * 2^n mod M,
    // that is the Montgomery representation of the operand x * t mod M, and t is the 64-bit constant
    #[inline]
    fn scalar_mul(res: &mut F, state: &mut [F], mut start_idx_cst: usize) {
        state.iter().for_each(|x| {
            let elem = x.mul(&P::MDS_CST[start_idx_cst]);
            start_idx_cst += 1;
            *res += &elem;
        });
    }
}


impl<F: PrimeField, P: PoseidonParameters<Fr = F>> SBox for PoseidonQuinticSBox<F, P> {
    type Field = F;
    type Parameters = P;

    #[inline]
    fn apply_full(state: &mut Vec<F>, last: bool) {
        // Apply the S-BOX to each of the elements of the state vector
        for i in 0..P::T {
            Self::exp_by_5(&mut state[i]);
        }

        if !last {
            // Perform the matrix mix
            Self::matrix_mix(state);
        }
    }

    #[inline]
    fn apply_partial(state: &mut Vec<F>) {
        // Apply S-BOX only to the first element of the state vector
        Self::exp_by_5(&mut state[0]);
        // Apply the matrix mix
        Self::matrix_mix(state);
    }
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> BatchSBox for PoseidonQuinticSBox<F, P> {}
impl<F: PrimeField, P: PoseidonParameters<Fr = F>> PoseidonBatchSBox<P> for PoseidonQuinticSBox<F, P> {}