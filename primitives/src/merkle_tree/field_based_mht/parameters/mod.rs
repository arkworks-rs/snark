#[cfg(feature = "mnt4_753")]
pub mod mnt4753;
#[cfg(feature = "mnt4_753")]
pub use self::mnt4753::*;

#[cfg(feature = "mnt6_753")]
pub mod mnt6753;
#[cfg(feature = "mnt6_753")]
pub use self::mnt6753::*;

#[cfg(feature = "tweedle")]
pub mod tweedle_dee;
#[cfg(feature = "tweedle")]
pub use self::tweedle_dee::*;

#[cfg(feature = "tweedle")]
pub mod tweedle_dum;
#[cfg(feature = "tweedle")]
pub use self::tweedle_dum::*;

#[cfg(feature = "bn_382")]
pub mod bn382;
#[cfg(feature = "bn_382")]
pub use self::bn382::*;

#[cfg(feature = "bn_382")]
pub mod bn382_dual;
#[cfg(feature = "bn_382")]
pub use self::bn382_dual::*;

use crate::FieldBasedHash;
use algebra::{PrimeField, ToConstraintField};

#[allow(dead_code)]
pub(crate) fn generate_phantom_merkle_root_from_magic_string<
    F: PrimeField,
    H: FieldBasedHash<Data = F>,
>(
    magic_string: &str,
) -> F {
    let magic_string_as_fes = magic_string.as_bytes().to_field_elements().unwrap();

    let mut digest = H::init_constant_length(magic_string_as_fes.len(), None);
    magic_string_as_fes.into_iter().for_each(|fe| {
        digest.update(fe);
    });
    digest.finalize().unwrap()
}

#[allow(dead_code)]
pub(crate) fn generate_mht_empty_nodes<F: PrimeField, H: FieldBasedHash<Data = F>>(
    arity: usize,
    max_height: usize,
    empty_leaf: F,
) -> Vec<F> {
    let mut empty_nodes = Vec::with_capacity(max_height);
    empty_nodes.push(empty_leaf.clone());

    let mut digest = H::init_constant_length(arity, None);
    let mut empty_node = empty_leaf;

    for _ in 1..max_height {
        for _ in 0..arity {
            digest.update(empty_node);
        }
        empty_node = digest.finalize().unwrap();
        //println!("{:?}", empty_node);
        empty_nodes.push(empty_node);
        digest.reset(None);
    }

    empty_nodes
}
