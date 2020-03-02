pub trait Flags: Default + Clone + Copy + Sized {
    fn u8_bitmask(&self) -> u8;
    fn from_u8(value: u8) -> Self;
    fn from_u8_remove_flags(value: &mut u8) -> Self;
    fn len() -> usize;
}

/// Flags to be encoded into the serialization.
/// The default flags (empty) should not change the binary representation.
#[derive(Default, Clone, Copy)]
pub struct EmptyFlags;

impl Flags for EmptyFlags {
    fn u8_bitmask(&self) -> u8 {
        0
    }

    fn from_u8(_value: u8) -> Self {
        EmptyFlags
    }

    fn from_u8_remove_flags(_value: &mut u8) -> Self {
        EmptyFlags
    }

    fn len() -> usize {
        0
    }
}

/// Flags to be encoded into the serialization.
/// The default flags (empty) should not change the binary representation.
#[derive(Default, Clone, Copy)]
pub struct SWFlags {
    pub y_sign:      bool,
    pub is_infinity: bool,
}

impl SWFlags {
    pub fn infinity() -> Self {
        SWFlags {
            y_sign:      false,
            is_infinity: true,
        }
    }

    pub fn y_sign(sign: bool) -> Self {
        SWFlags {
            y_sign:      sign,
            is_infinity: false,
        }
    }
}

impl Flags for SWFlags {
    fn u8_bitmask(&self) -> u8 {
        let mut mask = 0;
        if self.y_sign {
            mask |= 1 << 7;
        }
        if self.is_infinity {
            mask |= 1 << 6;
        }
        mask
    }

    fn from_u8(value: u8) -> Self {
        let x_sign = (value >> 7) & 1 == 1;
        let is_infinity = (value >> 6) & 1 == 1;
        SWFlags {
            y_sign: x_sign,
            is_infinity,
        }
    }

    fn from_u8_remove_flags(value: &mut u8) -> Self {
        let flags = Self::from_u8(*value);
        *value &= 0x3F;
        flags
    }

    /// Number of bits required for these flags.
    fn len() -> usize {
        2
    }
}

/// Flags to be encoded into the serialization.
/// The default flags (empty) should not change the binary representation.
#[derive(Default, Clone, Copy)]
pub struct EdwardsFlags {
    pub y_sign: bool,
}

impl EdwardsFlags {
    pub fn y_sign(sign: bool) -> Self {
        EdwardsFlags { y_sign: sign }
    }
}

impl Flags for EdwardsFlags {
    fn u8_bitmask(&self) -> u8 {
        let mut mask = 0;
        if self.y_sign {
            mask |= 1 << 7;
        }
        mask
    }

    fn from_u8(value: u8) -> Self {
        let x_sign = (value >> 7) & 1 == 1;
        EdwardsFlags { y_sign: x_sign }
    }

    fn from_u8_remove_flags(value: &mut u8) -> Self {
        let flags = Self::from_u8(*value);
        *value &= 0x7F;
        flags
    }

    /// Number of bits required for these flags.
    fn len() -> usize {
        1
    }
}
