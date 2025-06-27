use ark_std::vec::Vec;

use ark_ff::Field;

use crate::utils::IndexMap;

/// Interned field element:
///
/// * bit 0   – tag: 0 = inline constant, 1 = pooled
/// * bits 1-23  (23 b) – signed small constant (two’s complement)
/// * bits 24-63 (40 b) – index into the coefficient pool
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub struct InternedField(u32);

/// An interner for field elements.
#[derive(Debug, Clone, Default)]
pub struct FieldInterner<F: Field> {
    /// The map from field elements to their interned representation.
    pub map: IndexMap<F, InternedField>,
    /// The vector of field elements in the order they were added.
    pub vec: Vec<F>,
}

impl<F: Field> FieldInterner<F> {
    /// Creates a new field interner.
    #[inline]
    pub(super) fn new() -> Self {
        let mut result = Self {
            map: IndexMap::default(),
            vec: Vec::new(),
        };
        result.intern(F::ONE);
        result.intern(-F::ONE);
        result
    }

    #[inline(always)]
    fn intern(&mut self, value: F) -> InternedField {
        let index = self.vec.len();
        let interned = InternedField(index as u32);
        self.map.insert(value, interned);
        self.vec.push(value);
        interned
    }

    /// Interns a field element, returning its unique identifier.
    #[inline]
    pub(crate) fn get_or_intern(&mut self, value: F) -> InternedField {
        if value == F::ONE {
            InternedField(0)
        } else if let Some(&index) = self.map.get(&value) {
            index
        } else {
            self.intern(value)
        }
    }

    /// Returns the field element corresponding to the given identifier.
    #[inline(always)]
    pub(crate) fn value(&self, id: InternedField) -> Option<F> {
        if id == InternedField(0) {
            Some(F::ONE)
        } else if id == InternedField(1) {
            return Some(-F::ONE);
        } else {
            self.vec.get(id.0 as usize).copied()
        }
    }
}
