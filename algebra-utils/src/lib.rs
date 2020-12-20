#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public)]
#![cfg_attr(use_asm, feature(llvm_asm))]
#![cfg_attr(not(use_asm), forbid(unsafe_code))]
#![cfg_attr(use_asm, deny(unsafe_code))]

#[cfg(feature = "parallel")]
pub mod msm;
#[cfg(feature = "parallel")]
pub use self::msm::*;

#[cfg(feature = "fft")]
pub mod fft;
#[cfg(feature = "fft")]
pub use self::fft::*;

#[cfg(feature = "parallel")]
pub mod polycommit;
#[cfg(feature = "parallel")]
pub use self::polycommit::*;

pub type Error = Box<dyn std::error::Error>;