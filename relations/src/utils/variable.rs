use core::cmp::Ordering;

use super::linear_combination::LcIndex;

/// Represents the different kinds of variables present in a constraint system.
#[derive(Copy, Clone, PartialEq, Debug, Eq)]
pub enum Variable {
    /// Represents the "zero" constant.
    Zero,
    /// Represents of the "one" constant.
    One,
    /// Represents a public instance variable.
    Instance(usize),
    /// Represents a private witness variable.
    Witness(usize),
    /// Represents of a linear combination.
    SymbolicLc(LcIndex),
}

impl Variable {
    /// Create a `new` Zero instance variable.
    pub fn zero() -> Variable {
        Variable::Zero
    }

    /// Create a `new` One instance variable.
    #[inline]
    pub fn one() -> Variable {
        Variable::One
    }

    /// Is `self` the zero variable?
    #[inline]
    pub fn is_zero(&self) -> bool {
        matches!(self, Variable::Zero)
    }

    /// Is `self` the one variable?
    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, Variable::One)
    }

    /// Is `self` an instance variable?
    #[inline]
    pub fn is_instance(&self) -> bool {
        matches!(self, Variable::Instance(_))
    }

    /// Is `self` a witness variable?
    #[inline]
    pub fn is_witness(&self) -> bool {
        matches!(self, Variable::Witness(_))
    }

    /// Is `self` a linear combination?
    #[inline]
    pub fn is_lc(&self) -> bool {
        matches!(self, Variable::SymbolicLc(_))
    }

    /// Get the `LcIndex` in `self` if `self.is_lc()`.
    #[inline]
    pub fn get_lc_index(&self) -> Option<LcIndex> {
        match self {
            Variable::SymbolicLc(index) => Some(*index),
            _ => None,
        }
    }

    /// Returns `Some(usize)` if `!self.is_lc()`, and `None` otherwise.
    #[inline]
    pub fn get_index_unchecked(&self, witness_offset: usize) -> Option<usize> {
        match self {
            // The one variable always has index 0
            Variable::One => Some(0),
            Variable::Instance(i) => Some(*i),
            Variable::Witness(i) => Some(witness_offset + *i),
            _ => None,
        }
    }
}

impl PartialOrd for Variable {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Variable {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        use Variable::*;
        match (self, other) {
            (Zero, Zero) => Ordering::Equal,
            (One, One) => Ordering::Equal,
            (Zero, _) => Ordering::Less,
            (One, _) => Ordering::Less,
            (_, Zero) => Ordering::Greater,
            (_, One) => Ordering::Greater,

            (Instance(i), Instance(j)) | (Witness(i), Witness(j)) => i.cmp(j),
            (Instance(_), Witness(_)) => Ordering::Less,
            (Witness(_), Instance(_)) => Ordering::Greater,

            (SymbolicLc(i), SymbolicLc(j)) => i.cmp(j),
            (_, SymbolicLc(_)) => Ordering::Less,
            (SymbolicLc(_), _) => Ordering::Greater,
        }

    }
}
