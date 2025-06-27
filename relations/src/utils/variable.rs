/// Variables in [`ConstraintSystem`]s
#[derive(Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
#[must_use]
pub struct Variable(u64);

impl Variable {
    // [ tag: payload: 61 bits | 3 bits ]
    const TAG_BITS: u64 = 3;

    /// Bit position where the tag field starts (61 = 64 − 3).
    const TAG_SHIFT: u64 = 64 - Self::TAG_BITS;

    /// Mask for the payload (low 61 bits) when the tag is in the top byte.
    const PAYLOAD_MASK: u64 = (1u64 << Self::TAG_SHIFT) - 1;

    /// The zero variable.
    #[allow(non_upper_case_globals)]
    pub const Zero: Variable = Variable::pack_unchecked(0, 0);

    /// The one variable.
    #[allow(non_upper_case_globals)]
    pub const One: Variable = Variable::pack_unchecked(1, 0);

    /// The zero variable.
    #[inline(always)]
    pub const fn zero() -> Self {
        Self::Zero
    }

    /// Is `self` the zero variable?
    #[inline(always)]
    #[must_use]
    pub const fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Is `self` the one variable?
    #[inline(always)]
    #[must_use]
    pub const fn is_one(&self) -> bool {
        self.0 == Self::One.0
    }

    /// The `one` variable.
    #[inline(always)]
    pub const fn one() -> Self {
        Self::One
    }

    /// Construct an instance variable.
    #[inline(always)]
    pub const fn instance(i: usize) -> Self {
        Self::pack_unchecked(0b010, i as u64)
    }

    /// Is `self` an instance variable?
    #[inline(always)]
    #[must_use]
    pub const fn is_instance(self) -> bool {
        self.tag() == VarKind::Instance as u8
    }

    /// Construct a new witness variable.
    #[inline(always)]
    pub const fn witness(i: usize) -> Self {
        Self::pack_unchecked(0b011, i as u64)
    }

    /// Is `self` a witness variable?
    #[inline(always)]
    #[must_use]
    pub const fn is_witness(self) -> bool {
        self.tag() == VarKind::Witness as u8
    }

    /// Construct a symbolic linear combination variable.
    #[inline(always)]
    pub const fn symbolic_lc(i: usize) -> Self {
        Self::pack_unchecked(0b100, i as u64)
    }

    /// Is `self` a symbolic linear combination variable?
    #[inline(always)]
    #[must_use]
    pub const fn is_lc(self) -> bool {
        self.tag() == VarKind::SymbolicLc as u8
    }

    /// Get the `usize` in `self` if `self.is_lc()`.
    #[inline(always)]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub const fn get_lc_index(&self) -> Option<usize> {
        if self.is_lc() {
            Some(self.payload() as usize)
        } else {
            None
        }
    }

    /// Returns `Some(usize)` if `!self.is_lc()`, and `None` otherwise.
    #[inline(always)]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub const fn get_variable_index(&self, witness_offset: usize) -> Option<usize> {
        match self.kind() {
            // The one variable always has index 0
            VarKind::One => Some(0),
            VarKind::Instance => Some(self.payload() as usize),
            VarKind::Witness => Some(self.payload() as usize + witness_offset),
            _ => None,
        }
    }

    /// Returns the tag of the variable.
    #[inline(always)]
    const fn tag(self) -> u8 {
        (self.0 >> Self::TAG_SHIFT) as u8
    }

    /// Unconditionally returns the payload of the variable.
    /// Note that when `self.tag() == 0` or `self.tag() == 1`, the data
    /// value is not meaningful.
    #[inline(always)]
    const fn payload(self) -> u64 {
        self.0 & Self::PAYLOAD_MASK
    }

    /// What kind of variable is this?
    #[inline(always)]
    #[allow(unsafe_code)]
    pub const fn kind(self) -> VarKind {
        match self.tag() {
            0 => VarKind::Zero,
            1 => VarKind::One,
            2 => VarKind::Instance,
            3 => VarKind::Witness,
            4 => VarKind::SymbolicLc,
            _ => unsafe { core::hint::unreachable_unchecked() },
        }
    }

    /// If `self` is an instance, witness, or symbolic linear combination,
    /// returns the index of that variable.
    #[inline(always)]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub const fn index(self) -> Option<usize> {
        match self.kind() {
            VarKind::Zero | VarKind::One => None,
            _ => Some(self.payload() as usize),
        }
    }

    /// Does not check that the tag and payload are valid.
    const fn pack_unchecked(tag: u64, payload: u64) -> Self {
        debug_assert!(payload <= Self::PAYLOAD_MASK);
        Variable((tag << Self::TAG_SHIFT) | payload & Self::PAYLOAD_MASK)
    }

    #[cfg(test)]
    const fn new(kind: VarKind, index: usize) -> Self {
        match kind {
            VarKind::Zero => Self::Zero,
            VarKind::One => Self::One,
            VarKind::Instance => Self::instance(index),
            VarKind::Witness => Self::witness(index),
            VarKind::SymbolicLc => Self::symbolic_lc(index),
        }
    }
}

/// The kinds of variables that can be used in a constraint system.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
#[allow(missing_docs)]
#[must_use]
pub enum VarKind {
    Zero = 0,
    One = 1,
    Instance = 2,
    Witness = 3,
    SymbolicLc = 4,
}

impl core::fmt::Debug for Variable {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match (self.kind(), self.index()) {
            (VarKind::Zero, _) => f.write_str("Zero"),
            (VarKind::One, _) => f.write_str("One"),
            (k, Some(i)) => f.debug_tuple(&format!("{k:?}")).field(&i).finish(),
            _ => unreachable!(),
        }
    }
}

// Compile-time proof it really is 8 B.
const _: () = assert!(core::mem::size_of::<Variable>() == 8);

#[cfg(test)]
mod tests {
    // test PartialOrd and Ord vs Eq and PartialEq
    use super::*;

    use ark_std::rand::Rng;

    #[test]
    fn test_variable_ordering() {
        use core::cmp::Ordering::*;
        use VarKind::*;
        let mut rng = ark_std::test_rng();
        let kinds = [Zero, One, Instance, Witness, SymbolicLc];
        for this_kind in kinds {
            let this_payload: u32 = rng.gen();
            let this = Variable::new(this_kind, this_payload as usize);
            for other_kind in kinds {
                let other_1 = Variable::new(other_kind, this_payload as usize);

                let other_payload: u32 = rng.gen();
                let other_2 = Variable::new(other_kind, other_payload as usize);

                let eq_case_with_payload = || {
                    assert_eq!(this, other_1, "{this:?} != {other_1:?}");
                    if this_payload < other_payload {
                        assert!(this < other_2, "{this:?} >= {other_2:?}");
                    } else if this_payload > other_payload {
                        assert!(this > other_2, "{this:?} <= {other_2:?}");
                    } else {
                        assert_eq!(this, other_2, "{this:?} != {other_2:?}");
                    }
                    assert_eq!(this.cmp(&other_1), Equal);
                };
                let eq_case = || {
                    assert_eq!(this, other_1, "{this:?} != {other_1:?}");
                    assert_eq!(this, other_2, "{this:?} != {other_2:?}");
                    assert_eq!(this.cmp(&other_1), Equal);
                };
                let lt_case = || {
                    assert!(this < other_1, "{this:?} >= {other_1:?}");
                    assert!(this < other_2, "{this:?} >= {other_2:?}");
                };
                let gt_case = || {
                    assert!(this > other_1, "{this:?} <= {other_1:?}");
                    assert!(this > other_2, "{this:?} <= {other_2:?}");
                };
                match (this_kind, other_kind) {
                    (Zero, Zero) => eq_case(),
                    (One, One) => eq_case(),
                    (Instance, Instance) => eq_case_with_payload(),
                    (Witness, Witness) => eq_case_with_payload(),
                    (SymbolicLc, SymbolicLc) => eq_case_with_payload(),

                    (Zero, _) => lt_case(),
                    (_, Zero) => gt_case(),

                    (One, _) => lt_case(),
                    (_, One) => gt_case(),

                    (Instance, _) => lt_case(),
                    (_, Instance) => gt_case(),

                    (Witness, _) => lt_case(),
                    (_, Witness) => gt_case(),
                }
            }
        }
    }
}
