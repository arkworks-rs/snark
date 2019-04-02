use snark::SynthesisError;
////////////////////////////////////////////////////

///////////////////// CRYPTO PRIMITIVES /////////////////////////

pub mod commitment;
pub use self::commitment::*;

pub mod crh;
pub use self::crh::*;

pub mod prf;
pub use self::prf::*;

pub mod signature;
pub use self::signature::*;

pub mod verifier;
pub use self::verifier::*;

pub mod mht;
pub use self::mht::*;

////////////////////////////////////////////////////

///////////////////// DPC GADGETS /////////////////////////
pub mod dpc;

////////////////////////////////////////////////////

pub trait Assignment<T> {
    fn get(&self) -> Result<&T, SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(&self) -> Result<&T, SynthesisError> {
        match *self {
            Some(ref v) => Ok(v),
            None => Err(SynthesisError::AssignmentMissing),
        }
    }
}
