use core::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KeyGenError {
    RngFailure,
    InvalidSeed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SignError {
    /// All FIS slots rejected — probability < 2^{-24} for PRISM-128
    AllSlotsRejected,
    InvalidSecretKey,
    InvalidContext,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerifyError {
    InvalidSignature,
    InvalidPublicKey,
    Forgery,
}

impl fmt::Display for SignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SignError::AllSlotsRejected => write!(f, "all FIS slots rejected (P < 2^-24)"),
            SignError::InvalidSecretKey => write!(f, "malformed secret key"),
            SignError::InvalidContext => write!(f, "context too long"),
        }
    }
}

impl fmt::Display for VerifyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerifyError::InvalidSignature => write!(f, "malformed signature"),
            VerifyError::InvalidPublicKey => write!(f, "malformed public key"),
            VerifyError::Forgery => write!(f, "signature verification failed"),
        }
    }
}
