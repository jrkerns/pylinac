use ndarray::ShapeError;
use pyo3::PyErr;
use pyo3::exceptions::PyValueError;
use std::fmt::Display;
use std::string::FromUtf8Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    InvalidPixels,
    InvalidProperty,
    InvalidHistogram,
    InvalidCompressionIndicator,
    FailedDecompression,
    InvalidWidth,
    InvalidHeight,
    InvalidPixelBufferSize,
    InvalidLookupTableSize,
    InvalidOther(String),
}

impl core::error::Error for Error {}

impl Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidCompressionIndicator => {
                todo!()
            }
            Error::InvalidWidth => todo!(),
            Error::InvalidHeight => todo!(),
            Error::InvalidPixelBufferSize => {
                todo!()
            }
            Error::InvalidOther(val) => write!(f, "Failed: {}", val),
            Error::InvalidPixels => todo!(),
            Error::InvalidHistogram => todo!(),
            Error::InvalidProperty => todo!(),
            Error::FailedDecompression => todo!(),
            _ => todo!(),
        }
    }
}

impl From<ShapeError> for Error {
    fn from(value: ShapeError) -> Self {
        Self::InvalidOther(value.to_string())
    }
}

impl From<std::io::Error> for Error {
    fn from(value: std::io::Error) -> Self {
        Self::InvalidOther(value.to_string())
    }
}

impl From<FromUtf8Error> for Error {
    fn from(value: FromUtf8Error) -> Self {
        Self::InvalidOther(value.to_string())
    }
}

impl From<Error> for PyErr {
    fn from(value: Error) -> Self {
        match value {
            Error::InvalidPixels => PyValueError::new_err("Failed to read pixels for XIM image."),
            Error::InvalidProperty => {
                PyValueError::new_err("Failed to read property for XIM image.")
            }
            Error::InvalidHistogram => {
                PyValueError::new_err("Failed to read histogram for XIM image.")
            }
            Error::InvalidCompressionIndicator => {
                PyValueError::new_err("Unable to tell if the image is compressed or not.")
            }
            Error::FailedDecompression => PyValueError::new_err("Failed to decompress image."),
            Error::InvalidWidth => PyValueError::new_err("Invalid width for image."),
            Error::InvalidHeight => PyValueError::new_err("Invalid height for image"),
            Error::InvalidPixelBufferSize => PyValueError::new_err("Invalid pixel buffer size."),
            Error::InvalidLookupTableSize => PyValueError::new_err("Invalid lookup table size"),
            Error::InvalidOther(_) => PyValueError::new_err("Unknown error"),
        }
    }
}
