#![allow(
    clippy::upper_case_acronyms,
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::try_err,
    clippy::map_collect_result_unit,
    clippy::not_unsafe_ptr_arg_deref,
    clippy::suspicious_op_assign_impl,
    clippy::suspicious_arithmetic_impl,
    clippy::assertions_on_constants
)]

#[macro_use]
extern crate derivative;

pub mod crh;
pub use self::crh::*;

#[cfg(feature = "commitment")]
pub mod commitment;
#[cfg(feature = "commitment")]
pub use self::commitment::*;

#[cfg(feature = "merkle_tree")]
pub mod merkle_tree;
#[cfg(feature = "merkle_tree")]
pub use self::merkle_tree::*;

#[cfg(feature = "prf")]
pub mod prf;
#[cfg(feature = "prf")]
pub use self::prf::*;

#[cfg(feature = "signature")]
pub mod signature;
#[cfg(feature = "signature")]
pub use self::signature::*;

#[cfg(feature = "vrf")]
pub mod vrf;
#[cfg(feature = "vrf")]
pub use self::vrf::*;

#[cfg(feature = "nizk")]
pub mod nizk;
#[cfg(feature = "nizk")]
pub use self::nizk::*;
