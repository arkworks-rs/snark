pub mod delegable_dpc;
pub mod plain_dpc;

pub trait Assignment<T> {
    fn get(&self) -> Result<&T, r1cs_core::SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(&self) -> Result<&T, r1cs_core::SynthesisError> {
        match *self {
            Some(ref v) => Ok(v),
            None => Err(r1cs_core::SynthesisError::AssignmentMissing),
        }
    }
}
