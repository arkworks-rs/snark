use crate::darlin::pcd::error::PCDError;
use marlin::Error as MarlinError;
use poly_commit::Error as PCError;

#[derive(Debug)]
pub enum FinalDarlinError {
    MarlinError(MarlinError<PCError>),
    PCDError(PCDError),
    Other(String),
}

impl std::fmt::Display for FinalDarlinError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FinalDarlinError::MarlinError(err) => write!(f, "{}", err),
            FinalDarlinError::PCDError(err) => write!(f, "{}", err),
            FinalDarlinError::Other(err) => write!(f, "{}", err),
        }
    }
}

impl From<MarlinError<PCError>> for FinalDarlinError {
    fn from(err: MarlinError<PCError>) -> Self {
        FinalDarlinError::MarlinError(err)
    }
}

impl From<PCDError> for FinalDarlinError {
    fn from(err: PCDError) -> Self {
        FinalDarlinError::PCDError(err)
    }
}

impl std::error::Error for FinalDarlinError {}
