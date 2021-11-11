use crate::FieldBasedHashParameters;
use algebra::Field;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

pub trait SBox {
    type Field: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Field>;

    // Apply this SBox to the state, if performing a full round
    fn apply_full(state: &mut Vec<Self::Field>);

    // Apply this SBox to the state, if performing a partial round
    fn apply_partial(state: &mut Vec<Self::Field>);
}

pub trait BatchSBox: SBox {
    fn apply_full_batch(vec_state: &mut [Vec<Self::Field>]) {
        vec_state.par_iter_mut().for_each(|s| Self::apply_full(s));
    }

    fn apply_partial_batch(vec_state: &mut [Vec<Self::Field>]) {
        vec_state
            .par_iter_mut()
            .for_each(|s| Self::apply_partial(s));
    }
}
