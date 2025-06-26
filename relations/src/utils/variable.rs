use core::cmp::Ordering;

/// Variables in [`ConstraintSystem`]s
#[derive(Copy, Clone, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct Variable(u64);

impl Variable {
    // [ tag: 3 bits | payload: 61 bits ]
    const TAG_BITS: u64 = 3;
    const TAG_MASK: u64 = (1 << Self::TAG_BITS) - 1;
    const PAYLOAD_SHIFT: u64 = Self::TAG_BITS;

    /// The zero variable.
    #[allow(non_upper_case_globals)]
    pub const Zero: Variable = Variable(0);

    /// The one variable.
    #[allow(non_upper_case_globals)]
    pub const One: Variable = Variable(1);

    /// The zero variable.
    #[inline(always)]
    pub const fn zero() -> Self {
        Self(0)
    }

    /// Is `self` the zero variable?
    #[inline(always)]
    pub const fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Is `self` the one variable?
    #[inline(always)]
    pub const fn is_one(&self) -> bool {
        self.0 == 1
    }

    /// The `one` variable.
    #[inline(always)]
    pub const fn one() -> Self {
        Self(1)
    }

    /// Construct an instance variable.
    #[inline(always)]
    pub const fn Instance(i: usize) -> Self {
        Self::pack(0b010, i as u64)
    }

    /// Is `self` an instance variable?
    #[inline(always)]
    pub const fn is_instance(self) -> bool {
        matches!(self.kind(), VarKind::Instance)
    }

    /// Construct a new witness variable.
    #[inline(always)]
    pub const fn Witness(i: usize) -> Self {
        Self::pack(0b011, i as u64)
    }

    /// Is `self` a witness variable?
    #[inline(always)]
    pub const fn is_witness(self) -> bool {
        matches!(self.kind(), VarKind::Witness)
    }

    /// Construct a symbolic linear combination variable.
    #[inline(always)]
    pub const fn SymbolicLc(i: usize) -> Self {
        Self::pack(0b100, i as u64)
    }

    /// Is `self` a symbolic linear combination variable?
    #[inline(always)]
    pub const fn is_lc(self) -> bool {
        matches!(self.kind(), VarKind::SymbolicLc)
    }

    /// Get the `usize` in `self` if `self.is_lc()`.
    #[inline(always)]
    pub const fn get_lc_index(&self) -> Option<usize> {
        if self.is_lc() {
            Some(self.index().unwrap() as usize)
        } else {
            None
        }
    }

    /// Returns `Some(usize)` if `!self.is_lc()`, and `None` otherwise.
    #[inline(always)]
    pub const fn get_variable_index(&self, witness_offset: usize) -> Option<usize> {
        match self.kind() {
            // The one variable always has index 0
            VarKind::One => Some(0),
            VarKind::Instance => self.index(),
            VarKind::Witness => Some(self.index().unwrap() + witness_offset),
            _ => None,
        }
    }

    /// What kind of variable is this?
    #[inline(always)]
    pub const fn kind(self) -> VarKind {
        match self.0 & Self::TAG_MASK {
            0 => VarKind::Zero,
            1 => VarKind::One,
            2 => VarKind::Instance,
            3 => VarKind::Witness,
            4 => VarKind::SymbolicLc,
            _ => unreachable!(),
        }
    }

    /// If `self` is an instance, witness, or symbolic linear combination,
    /// returns the index of that variable.
    #[inline(always)]
    pub const fn index(self) -> Option<usize> {
        match self.kind() {
            VarKind::Zero | VarKind::One => None,
            _ => Some((self.0 >> Self::PAYLOAD_SHIFT) as usize),
        }
    }

    const fn pack(tag: u64, payload: u64) -> Self {
        debug_assert!(
            payload >> (64 - Self::PAYLOAD_SHIFT) == 0,
            "payload too large"
        );
        Self((payload << Self::PAYLOAD_SHIFT) | tag)
    }
}

/// The kinds of variables that can be used in a constraint system.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
#[allow(missing_docs)]
pub enum VarKind {
    Zero,
    One,
    Instance,
    Witness,
    SymbolicLc,
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

impl PartialOrd for Variable {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Variable {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        use VarKind::*;
        match (self.kind(), other.kind()) {
            (Zero, Zero) => Ordering::Equal,
            (One, One) => Ordering::Equal,
            (Zero, _) => Ordering::Less,
            (One, _) => Ordering::Less,
            (_, Zero) => Ordering::Greater,
            (_, One) => Ordering::Greater,

            (Instance, Instance) | (Witness, Witness) | (SymbolicLc, SymbolicLc) => {
                self.index().cmp(&other.index())
            },

            (Instance, Witness) => Ordering::Less,
            (Witness, Instance) => Ordering::Greater,
            (_, SymbolicLc) => Ordering::Less,
            (SymbolicLc, _) => Ordering::Greater,
        }
    }
}
