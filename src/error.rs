use ocelot::Error as OcelotError;
use std::error;
use std::fmt;

#[derive(Debug)]
pub enum OleError {
    SenderError(&'static str),
    ReceiverError(&'static str),
    IOError(std::io::Error),
    OTError(OcelotError),
}

impl fmt::Display for OleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            OleError::IOError(e) => write!(f, "{}", e),
            OleError::SenderError(s) => write!(f, "{}", s),
            OleError::ReceiverError(s) => write!(f, "{}", s),
            OleError::OTError(e) => write!(f, "{}", e),
        }

    }
}

impl error::Error for OleError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            OleError::IOError(e) => Some(e),
            // OleError::OTError(e) => Some(e), doesn't work, as OcelotError doesn't implement Error
            _ => None,
        }
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
