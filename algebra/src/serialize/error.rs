use core::fmt;
use crate::fake_io;

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
    IoError(fake_io::Error),
}

impl From<fake_io::Error> for SerializationError {
    fn from(e: fake_io::Error) -> SerializationError {
        SerializationError::IoError(e)
    }
}

impl fmt::Display for SerializationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        if let SerializationError::IoError(_) = self {
            write!(f, "I/O error")
        } else {
            let description = match self {
                SerializationError::BufferWrongSize => {
                    "the output buffer is of the wrong size"
                },
                SerializationError::NotEnoughSpace => {
                    "the last byte does not have enough space to encode the extra info bits"
                },
                SerializationError::InvalidData => {
                    "the input buffer contained invalid data"
                },
                SerializationError::ExtraInfoWrongSize => {
                    "extra info is of the wrong size"
                },
                SerializationError::IoError(_) => "encountered an I/O error",
            };
            write!(f, "{}", description)
        }
    }
}
