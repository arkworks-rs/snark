use crate::io;
use core::fmt;

/// This is an error that could occur during serialization
#[derive(Debug)]
pub enum SerializationError {
    /// During serialization, the output buffer was of the wrong size.
    BufferWrongSize,
    /// During serialization, we didn't have enough space to write extra info.
    NotEnoughSpace,
    /// During serialization, the data was invalid.
    InvalidData,
    /// During serialization, extra info was of the wrong size.
    ExtraInfoWrongSize,
    /// During serialization, we countered an I/O error.
    IoError(io::Error),
}

#[cfg(feature = "std")]
impl std::error::Error for SerializationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl From<io::Error> for SerializationError {
    fn from(e: io::Error) -> SerializationError {
        SerializationError::IoError(e)
    }
}

impl fmt::Display for SerializationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        match self {
            SerializationError::BufferWrongSize => {
                write!(f, "the output buffer is of the wrong size")
            },
            SerializationError::NotEnoughSpace => write!(
                f,
                "the last byte does not have enough space to encode the extra info bits"
            ),
            SerializationError::InvalidData => write!(f, "the input buffer contained invalid data"),
            SerializationError::ExtraInfoWrongSize => write!(f, "extra info is of the wrong size"),
            SerializationError::IoError(err) => write!(f, "I/O error: {:?}", err),
        }
    }
}
