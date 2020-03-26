#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[cfg(any(feature = "groth16", feature = "gm17"))]
#[macro_use]
extern crate bench_utils;

#[cfg(feature = "groth16")]
pub mod groth16;

#[cfg(feature = "gm17")]
pub mod gm17;
