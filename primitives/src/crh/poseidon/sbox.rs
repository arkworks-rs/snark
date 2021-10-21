use crate::{BatchSBox, PoseidonParameters, SBox};
use algebra::PrimeField;
use std::marker::PhantomData;

/// S-Box: S(x) = x^-1
#[derive(Debug)]
pub struct PoseidonInverseSBox<F: PrimeField, P: PoseidonParameters<Fr = F>> {
    _field: PhantomData<F>,
    _parameters: PhantomData<P>,
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> SBox for PoseidonInverseSBox<F, P> {
    type Field = F;
    type Parameters = P;

    // Uses batch inversion on the state.
    #[inline]
    fn apply_full(state: &mut Vec<F>) {
        // Apply the S-BOX to each of the elements of the state vector
        // Use batch inversion
        let mut w: Vec<F> = Vec::new();
        let mut accum_prod = F::one();

        w.push(accum_prod);

        // Calculate the intermediate partial products
        for d in state.iter() {
            accum_prod = accum_prod * d;
            w.push(accum_prod);
        }

        if accum_prod == F::zero() {
            // At least one of the S-Boxes is zero
            // Calculate inverses individually
            for d in state.iter_mut() {
                // The S-BOX is an inversion function
                if *d != F::zero() {
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

    #[inline]
    fn apply_partial(state: &mut Vec<F>) {
        if state[0] != F::zero() {
            state[0] = state[0].inverse().unwrap();
        }
    }
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> BatchSBox for PoseidonInverseSBox<F, P> {
    // Uses batch inversion across all instances in the batch.
    fn apply_full_batch(vec_state: &mut [Vec<F>]) {
        // Apply the S-BOX to each of the elements of the state vector
        // Use batch inversion
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
            let mut idx: i64 = w.len() as i64 - P::R as i64;
            for i in (0..vec_state.len()).rev() {
                for j in (0..P::T).rev() {
                    let vec_1 = vec_state[i][j].clone();
                    vec_state[i][j] = w_bar * &w[idx as usize];
                    w_bar = w_bar * &vec_1;
                    idx -= 1;
                }
            }
        }
    }

    // Uses batch inversion across all instances in the batch.
    fn apply_partial_batch(vec_state: &mut [Vec<F>]) {
        // Apply the S-BOX to the first elements of each of the state vector
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
            // Use batch inversion
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - P::R as i64;
            for i in (0..vec_state.len()).rev() {
                let vec_1 = vec_state[i][0].clone();
                vec_state[i][0] = w_bar * &w[idx as usize];
                w_bar = w_bar * &vec_1;
                idx -= 1;
            }
        }
    }
}

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

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> SBox for PoseidonQuinticSBox<F, P> {
    type Field = F;
    type Parameters = P;

    #[inline]
    fn apply_full(state: &mut Vec<F>) {
        // Apply the S-BOX to each of the elements of the state vector
        for i in 0..P::T {
            Self::exp_by_5(&mut state[i]);
        }
    }

    #[inline]
    fn apply_partial(state: &mut Vec<F>) {
        // Apply S-BOX only to the first element of the state vector
        Self::exp_by_5(&mut state[0]);
    }
}

impl<F: PrimeField, P: PoseidonParameters<Fr = F>> BatchSBox for PoseidonQuinticSBox<F, P> {}
