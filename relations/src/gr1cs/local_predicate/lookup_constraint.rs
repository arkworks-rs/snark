use core::marker::PhantomData;

use ark_ff::Field;
use ark_std::vec::Vec;

use super::Evaluatable;

// TODO
#[derive(Debug, Clone)]
pub struct LookupConstraint<F: Field> {
    _marker: PhantomData<F>,
}

impl<F: Field> Evaluatable<F> for LookupConstraint<F> {
    fn evaluate(&self, _point: Vec<F>) -> F {
        unimplemented!()
    }
}
