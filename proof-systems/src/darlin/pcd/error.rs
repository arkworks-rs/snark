#[derive(Debug)]
pub enum PCDError {
    FailedSuccinctVerification(String),
    FailedHardVerification(String),
    MissingSystemInputs(String),
    MissingUserInputs(String),
}

impl std::fmt::Display for PCDError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PCDError::FailedSuccinctVerification(err) => {
                write!(f, "Succinct check failed: {}", err)
            }
            PCDError::FailedHardVerification(err) => write!(f, "Hard check failed: {}", err),
            PCDError::MissingSystemInputs(missing_field) => {
                write!(f, "Unable to retrieve system input: {}", missing_field)
            }
            PCDError::MissingUserInputs(missing_field) => {
                write!(f, "Unable to retrieve user input: {}", missing_field)
            }
        }
    }
}

impl std::error::Error for PCDError {}
