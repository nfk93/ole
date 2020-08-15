use ocelot::Error as OcelotError;
use std::error;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum OleError {
    TODOError,
    IOError(std::io::Error),
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

impl From<std::io::Error> for OleError {
    fn from(e: std::io::Error) -> OleError {
        OleError::IOError(e)
    }
}
