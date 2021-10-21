use crate::Error;

pub trait ToBits {
    /// Serialize `self` into a bit vector using a BigEndian bit order representation.
    fn write_bits(&self) -> Vec<bool>;
}

pub trait FromBits: Sized {
    /// Reads `self` from `bits`, where `bits` are expected to be
    /// in a BigEndian bit order representation.
    fn read_bits(bits: Vec<bool>) -> Result<Self, Error>;
}

pub trait ToCompressedBits {
    fn compress(&self) -> Vec<bool>;
}

pub trait FromCompressedBits: Sized {
    fn decompress(compressed: Vec<bool>) -> Result<Self, Error>;
}

#[derive(Debug)]
pub enum BitSerializationError {
    InvalidFieldElement(String),
    UndefinedSqrt,
    NotPrimeOrder,
    NotOnCurve,
    NotInCorrectSubgroup,
    InvalidFlags,
}

impl std::fmt::Display for BitSerializationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            BitSerializationError::InvalidFieldElement(s) => s.to_owned(),
            BitSerializationError::UndefinedSqrt => "square root doesn't exist in field".to_owned(),
            BitSerializationError::NotPrimeOrder => {
                "point is not in the prime order subgroup".to_owned()
            }
            BitSerializationError::NotOnCurve => "point is not on curve".to_owned(),
            BitSerializationError::NotInCorrectSubgroup => {
                "point is not in the correct subgroup".to_owned()
            }
            BitSerializationError::InvalidFlags => "illegal flags combination".to_owned(),
        };
        write!(f, "{}", msg)
    }
}

impl std::error::Error for BitSerializationError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}
