use ocelot::Error as OcelotError;
use std::error;
use std::fmt;

#[derive(Debug)]
pub enum OleError {
    TODOError,
    OTError(OcelotError),
}

impl fmt::Display for OleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "something went wrong")
    }
}

impl error::Error for OleError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

impl From<OcelotError> for OleError {
    fn from(e: OcelotError) -> OleError {
        OleError::OTError(e)
    }
}
