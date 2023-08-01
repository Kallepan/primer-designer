use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub struct GenericError {
    pub message: String,
}

#[allow(dead_code)]
impl GenericError {
    pub fn new(message: &str) -> GenericError {
        GenericError {
            message: message.to_string(),
        }
    }
}

impl fmt::Display for GenericError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for GenericError {
    fn description(&self) -> &str {
        &self.message
    }
}