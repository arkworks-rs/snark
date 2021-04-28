use marlin::Error as MarlinError;
use poly_commit::Error as PCError;
use crate::darlin::pcd::error::PCDError;

#[derive(Debug)]
pub enum FinalDarlinError {
    MarlinError(MarlinError<PCError>),
    PCDError(PCDError),
    Other(String),
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