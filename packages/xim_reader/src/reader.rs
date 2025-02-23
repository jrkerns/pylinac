use byteorder::ByteOrder;
use numpy::PyArray2;
use numpy::PyArrayMethods;
use pyo3::{IntoPyObject, PyResult, prelude::Bound, pyclass, pymethods};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use pyo3_stub_gen::impl_stub_type;
use std::collections::HashMap;
use std::fmt::Debug;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use crate::error::{Error, Result};
use byteorder::{LE, ReadBytesExt};
use ndarray::ArrayViewMut2;

fn read_i32_as_usize(mut reader: impl Read, err: Error) -> Result<usize> {
    let val = reader.read_i32::<LE>()?;
    usize::try_from(val).map_err(|_err| err)
}

fn read_i32_into_buf<const N: usize>(mut reader: impl Read) -> Result<[i32; N]> {
    let mut values = [0i32; N];
    let _ = reader.read_i32_into::<LE>(&mut values)?;
    Ok(values)
}

fn read_i32_into_vec(mut reader: impl Read, n: usize) -> Result<Vec<i32>> {
    let mut values = vec![0i32; n];
    let _ = reader.read_i32_into::<LE>(&mut values)?;
    Ok(values)
}

trait ReadFromReader<I> {
    fn read_type_into<B: ByteOrder>(&mut self, dst: &mut [I]) -> Result<()>;
}

impl<R: ReadBytesExt> ReadFromReader<i8> for R {
    fn read_type_into<B: ByteOrder>(&mut self, dst: &mut [i8]) -> Result<()> {
        Ok(self.read_i8_into(dst)?)
    }
}

impl<R: ReadBytesExt> ReadFromReader<i16> for R {
    fn read_type_into<B: ByteOrder>(&mut self, dst: &mut [i16]) -> Result<()> {
        Ok(self.read_i16_into::<B>(dst)?)
    }
}

impl<R: ReadBytesExt> ReadFromReader<i32> for R {
    fn read_type_into<B: ByteOrder>(&mut self, dst: &mut [i32]) -> Result<()> {
        Ok(self.read_i32_into::<B>(dst)?)
    }
}

impl<R: ReadBytesExt> ReadFromReader<i64> for R {
    fn read_type_into<B: ByteOrder>(&mut self, dst: &mut [i64]) -> Result<()> {
        Ok(self.read_i64_into::<B>(dst)?)
    }
}

#[derive(IntoPyObject)]
pub enum XIMArray<'py> {
    #[pyo3(transparent)]
    Int8(Bound<'py, PyArray2<i8>>),
    #[pyo3(transparent)]
    Int16(Bound<'py, PyArray2<i16>>),
    #[pyo3(transparent)]
    Int32(Bound<'py, PyArray2<i32>>),
    #[pyo3(transparent)]
    Int64(Bound<'py, PyArray2<i64>>),
}

impl_stub_type!(XIMArray<'_> = PyArray2<i8> | PyArray2<i16> |PyArray2<i32> |PyArray2<i64>);

#[gen_stub_pyclass]
#[pyclass]
pub struct XIMImage {
    #[pyo3(get)]
    pub header: XIMHeader,
    pixel_data: PixelDataSupported,
    pub histogram: XIMHistogram,
    pub properties: XIMProperties,
}

#[gen_stub_pymethods]
#[pymethods]
impl XIMImage {
    #[new]
    pub fn new(image_path: PathBuf) -> PyResult<Self> {
        let file = File::open(image_path)?;
        let mut reader = BufReader::new(file);
        let header = XIMHeader::from_reader(&mut reader)?;
        let pixel_data = if header.is_compressed {
            PixelDataSupported::from_compressed(&mut reader, &header)?
        } else {
            PixelDataSupported::from_uncompressed(&mut reader, &header)?
        };
        let histogram = XIMHistogram::from_reader(&mut reader)?;
        let properties = XIMProperties::from_reader(&mut reader)?;
        Ok(Self {
            header,
            pixel_data,
            histogram,
            properties,
        })
    }

    #[getter]
    pub fn numpy<'py>(this: Bound<'py, Self>) -> XIMArray<'py> {
        match &this.borrow().pixel_data {
            PixelDataSupported::Int8(pixel_data) => {
                let array = &pixel_data.0;
                unsafe {
                    let pyarray = PyArray2::borrow_from_array(array, this.into_any());
                    pyarray.readwrite().make_nonwriteable();
                    XIMArray::Int8(pyarray)
                }
            }
            PixelDataSupported::Int16(pixel_data) => {
                let array = &pixel_data.0;
                unsafe {
                    let pyarray = PyArray2::borrow_from_array(array, this.into_any());
                    pyarray.readwrite().make_nonwriteable();
                    XIMArray::Int16(pyarray)
                }
            }
            PixelDataSupported::Int32(pixel_data) => {
                let array = &pixel_data.0;
                unsafe {
                    let pyarray = PyArray2::borrow_from_array(array, this.into_any());
                    pyarray.readwrite().make_nonwriteable();
                    XIMArray::Int32(pyarray)
                }
            }
            PixelDataSupported::Int64(pixel_data) => {
                let array = &pixel_data.0;
                unsafe {
                    let pyarray = PyArray2::borrow_from_array(array, this.into_any());
                    pyarray.readwrite().make_nonwriteable();
                    XIMArray::Int64(pyarray)
                }
            }
        }
    }

    #[getter]
    pub fn histogram(&self) -> Vec<i32> {
        self.histogram.histogram.clone()
    }

    #[getter]
    pub fn properties(&self) -> HashMap<String, PropertyValue> {
        self.properties.properties.clone()
    }
}

#[derive(Debug, Clone)]
#[gen_stub_pyclass]
#[pyclass]
pub struct XIMHeader {
    #[pyo3(get)]
    pub identifier: String,
    #[pyo3(get)]
    pub version: i32,
    #[pyo3(get)]
    pub width: i32,
    #[pyo3(get)]
    pub height: i32,
    #[pyo3(get)]
    pub bits_per_pixel: i32,
    #[pyo3(get)]
    pub bytes_per_pixel: i32,
    #[pyo3(get)]
    pub is_compressed: bool,
}

#[derive(Debug, Clone)]
pub struct PixelData<I>(ndarray::Array2<I>);

#[derive(Debug, Clone)]
pub enum PixelDataSupported {
    Int8(PixelData<i8>),
    Int16(PixelData<i16>),
    Int32(PixelData<i32>),
    Int64(PixelData<i64>),
}

#[derive(Debug, Clone)]
#[gen_stub_pyclass]
#[pyclass]
pub struct XIMHistogram {
    #[pyo3(get)]
    pub histogram: Vec<i32>,
}

#[derive(Debug, Clone)]
enum PropertyType {
    Integer,
    Double,
    String,
    DoubleArray,
    IntegerArray,
}

#[derive(Debug, Clone, IntoPyObject)]
pub enum PropertyValue {
    #[pyo3(transparent)]
    Integer(i32),
    #[pyo3(transparent)]
    Double(f64),
    #[pyo3(transparent)]
    String(String),
    #[pyo3(transparent)]
    DoubleArray(Vec<f64>),
    #[pyo3(transparent)]
    IntegerArray(Vec<i32>),
}

impl_stub_type!(PropertyValue = i32 | f64 | String | Vec<f64> | Vec<i32>);

#[derive(Debug, Clone)]
pub struct Property {
    pub property_name: String,
    pub property_value: PropertyValue,
}

#[derive(Debug, Clone)]
#[gen_stub_pyclass]
#[pyclass]
pub struct XIMProperties {
    pub properties: HashMap<String, PropertyValue>,
}

impl XIMHeader {
    pub fn from_reader<R: Read>(reader: &mut R) -> Result<Self> {
        let identifier = {
            let mut buf = [0u8; 8];
            reader.read_exact(&mut buf)?;
            String::from_utf8(buf.to_vec())?
        };

        let [
            version,
            width,
            height,
            bits_per_pixel,
            bytes_per_pixel,
            compression,
        ] = read_i32_into_buf(reader)?;

        let is_compressed = match compression {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(Error::InvalidCompressionIndicator),
        }?;

        Ok(Self {
            identifier,
            version,
            width,
            height,
            bits_per_pixel,
            bytes_per_pixel,
            is_compressed,
        })
    }

    pub fn width(&self) -> Result<usize> {
        usize::try_from(self.width).map_err(|_err| Error::InvalidWidth)
    }

    pub fn height(&self) -> Result<usize> {
        usize::try_from(self.height).map_err(|_err| Error::InvalidHeight)
    }
}

impl<I> PixelData<I> {
    pub fn new(array: ndarray::Array2<I>) -> Self {
        Self(array)
    }
}

impl PixelDataSupported {
    fn read_to_arr<I, R>(mut reader: R, width: usize, height: usize) -> Result<PixelData<I>>
    where
        I: num_traits::ConstZero + Clone + Copy,
        R: Read + ReadFromReader<I>,
    {
        let array = {
            let mut data: Vec<I> = vec![I::ZERO; width * height];
            let _ = ReadFromReader::<I>::read_type_into::<LE>(&mut reader, &mut data)
                .map_err(|_err| Error::InvalidPixels)?;
            ndarray::Array2::from_shape_vec((width, height), data)?
        };
        Ok(PixelData::new(array))
    }

    pub fn from_uncompressed(mut reader: impl Read, header: &XIMHeader) -> Result<Self> {
        let num_bytes = header.bytes_per_pixel;

        let width = header.width()?;
        let height = header.height()?;

        let _pixel_buffer_size = read_i32_as_usize(&mut reader, Error::InvalidPixelBufferSize)?;

        match num_bytes {
            1 => Self::read_to_arr(&mut reader, width, height).map(Self::Int8),
            2 => Self::read_to_arr(&mut reader, width, height).map(Self::Int16),
            4 => Self::read_to_arr(&mut reader, width, height).map(Self::Int32),
            8 => Self::read_to_arr(&mut reader, width, height).map(Self::Int64),
            _ => todo!(),
        }
    }

    pub fn parse_lookup(lookup_table: Vec<u8>) -> Vec<u8> {
        let num_bytes_table = lookup_table
            .into_iter()
            .flat_map(|vals| {
                vec![
                    (vals & 0b00000011),
                    (vals & 0b00001100) >> 2,
                    (vals & 0b00110000) >> 4,
                    (vals & 0b11000000) >> 6,
                ]
            })
            .map(|val| 1u8 << val)
            .collect::<Vec<u8>>();
        num_bytes_table
    }

    fn decompress_array<I>(
        compressed_pixel_buffer: Vec<u8>,
        num_bytes_table: impl Iterator<Item = u8>,
        width: usize,
        height: usize,
    ) -> Result<PixelData<I>>
    where
        I: num_traits::ConstZero
            + TryFrom<i32>
            + TryFrom<i8>
            + TryFrom<i16>
            + Clone
            + Copy
            + num_traits::WrappingAdd
            + num_traits::WrappingSub,
    {
        let (mut uncompressed_buffer, mut compressed_diffs) =
            compressed_pixel_buffer.split_at((width + 1) * 4);

        let initial_uncompressed = {
            let initial_uncompressed = read_i32_into_vec(&mut uncompressed_buffer, width + 1)?;
            initial_uncompressed
                .into_iter()
                .map(|val| I::try_from(val).map_err(|_err| Error::InvalidPixels))
                .collect::<Result<Vec<I>>>()?
        };

        let differences = num_bytes_table
            .map(|num_bytes| match num_bytes {
                1 => compressed_diffs
                    .read_i8()
                    .map_err(|_err| Error::InvalidPixels)
                    .and_then(|x| I::try_from(x).map_err(|_err| Error::InvalidPixels)),
                2 => compressed_diffs
                    .read_i16::<LE>()
                    .map_err(|_err| Error::InvalidPixels)
                    .and_then(|x| I::try_from(x).map_err(|_err| Error::InvalidPixels)),
                4 => compressed_diffs
                    .read_i32::<LE>()
                    .map_err(|_err| Error::InvalidPixels)
                    .and_then(|x| I::try_from(x).map_err(|_err| Error::InvalidPixels)),
                _ => todo!(),
            })
            .collect::<Result<Vec<I>>>()?;

        let array = {
            let uncompressed_data = [initial_uncompressed, differences].concat();
            let mut array = ndarray::Array2::from_shape_vec((height, width), uncompressed_data)?;
            Self::decompress_diffs(array.view_mut())?;
            array
        };

        Ok(PixelData::new(array))
    }

    pub fn from_compressed(mut reader: impl Read, header: &XIMHeader) -> Result<Self> {
        let width = header.width()?;
        let height = header.height()?;

        let lookup_table_size = read_i32_as_usize(&mut reader, Error::InvalidLookupTableSize)?;

        let lookup_table = {
            let lookup_table: Vec<u8> = {
                let mut buf = vec![0u8; lookup_table_size];
                let _ = reader.read_exact(&mut buf)?;
                buf
            };
            let num_bytes_table = Self::parse_lookup(lookup_table);
            let full_len = num_bytes_table.len();
            num_bytes_table
                .into_iter()
                .take(full_len - (width * (height - 1)) % 4 - 1)
        };

        let compressed_pixel_buffer_size =
            read_i32_as_usize(&mut reader, Error::InvalidPixelBufferSize)?;

        let compressed_pixel_buffer = {
            let mut buf = vec![0; compressed_pixel_buffer_size];
            let _ = reader.read_exact(&mut buf);
            buf
        };

        let pixel_data = match header.bytes_per_pixel {
            1 => {
                let arr =
                    Self::decompress_array(compressed_pixel_buffer, lookup_table, width, height)?;
                PixelDataSupported::Int8(arr)
            }
            2 => {
                let arr =
                    Self::decompress_array(compressed_pixel_buffer, lookup_table, width, height)?;
                PixelDataSupported::Int16(arr)
            }
            4 => {
                let arr =
                    Self::decompress_array(compressed_pixel_buffer, lookup_table, width, height)?;
                PixelDataSupported::Int32(arr)
            }
            _ => todo!(),
        };

        let _uncompressed_buffer_size = reader.read_i32::<LE>()?;
        Ok(pixel_data)
    }

    pub fn decompress_diffs<I>(mut compressed_arr: ArrayViewMut2<I>) -> Result<ArrayViewMut2<I>>
    where
        I: num_traits::WrappingAdd + num_traits::WrappingSub + Copy,
    {
        let width = compressed_arr.ncols();

        let arr = compressed_arr
            .as_slice_mut()
            .ok_or(Error::FailedDecompression)?;
        let first_index = width + 1;
        for i in first_index..arr.len() {
            let [left, above, upper_left] = (|| {
                Some([
                    *arr.get(i - 1)?,
                    *arr.get(i - width)?,
                    *arr.get(i - width - 1)?,
                ])
            })()
            .ok_or(Error::FailedDecompression)?;

            let diff = arr.get_mut(i).ok_or(Error::FailedDecompression)?;
            *diff = diff
                .wrapping_add(&left)
                .wrapping_add(&above)
                .wrapping_sub(&upper_left);
        }
        Ok(compressed_arr)
    }
}

impl XIMHistogram {
    pub fn from_reader<R: Read>(mut reader: R) -> Result<Self> {
        let number_of_bins = read_i32_as_usize(&mut reader, Error::InvalidHistogram)?;

        let histogram = read_i32_into_vec(&mut reader, number_of_bins)?;
        Ok(Self { histogram })
    }
}

impl Property {
    pub fn from_reader<R: Read>(mut reader: R) -> Result<Self> {
        let property_name_length = read_i32_as_usize(&mut reader, Error::InvalidProperty)?;

        let property_name = {
            let mut buf = vec![0u8; property_name_length];
            let _ = reader.read_exact(&mut buf)?;
            String::from_utf8(buf)?
        };

        let property_type = {
            let val = reader.read_i32::<LE>()?;
            match val {
                0 => PropertyType::Integer,
                1 => PropertyType::Double,
                2 => PropertyType::String,
                4 => PropertyType::DoubleArray,
                5 => PropertyType::IntegerArray,
                _ => todo!(),
            }
        };

        let property_value = match property_type {
            PropertyType::Integer => {
                let value = reader.read_i32::<LE>()?;
                PropertyValue::Integer(value)
            }
            PropertyType::Double => {
                let value = reader.read_f64::<LE>()?;
                PropertyValue::Double(value)
            }
            PropertyType::String => {
                let value_len = read_i32_as_usize(&mut reader, Error::InvalidProperty)?;

                let mut value = vec![0u8; value_len];
                reader.read_exact(&mut value)?;
                let value = String::from_utf8(value)?;
                PropertyValue::String(value)
            }
            PropertyType::DoubleArray => {
                let value_len = {
                    let len_bytes = read_i32_as_usize(&mut reader, Error::InvalidProperty)?;
                    len_bytes.checked_div(8).ok_or(Error::InvalidProperty)?
                };

                let mut value = vec![0f64; value_len];
                let _ = reader.read_f64_into::<LE>(&mut value)?;
                PropertyValue::DoubleArray(value)
            }
            PropertyType::IntegerArray => {
                let value_len = {
                    let len_bytes = read_i32_as_usize(&mut reader, Error::InvalidProperty)?;
                    len_bytes.checked_div(4).ok_or(Error::InvalidProperty)?
                };

                let mut value = vec![0i32; value_len];
                let _ = reader.read_i32_into::<LE>(&mut value)?;

                PropertyValue::IntegerArray(value)
            }
        };

        Ok(Self {
            property_name,
            property_value,
        })
    }
}

impl XIMProperties {
    pub fn from_reader<R: Read>(mut reader: R) -> Result<Self> {
        let num_properties = read_i32_as_usize(&mut reader, Error::InvalidProperty)?;

        let mut properties = HashMap::with_capacity(num_properties);

        for _ in 0..num_properties {
            let property = Property::from_reader(&mut reader)?;
            properties.insert(property.property_name, property.property_value);
        }
        Ok(Self { properties })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decompression() {
        let input: [i16; 8] = [4, 3, 10, 1, 10, 30, 20, 40];
        let mut input_array = ndarray::Array2::from_shape_vec((4, 2), input.to_vec()).unwrap();
        let calculated_output = PixelDataSupported::decompress_diffs(input_array.view_mut())
            .expect("Failed to decompress diffs");
        let output =
            ndarray::Array2::from_shape_vec((4, 2), vec![4, 3, 10, 10, 27, 57, 94, 164]).unwrap();
        assert_eq!(calculated_output, output);
    }
    #[test]
    fn test_parse_lookup() {
        let test: Vec<u8> = vec![1, 10, 30, 20, 40];
        println!("{:#010b}", test.get(1).unwrap());
        let test = test
            .into_iter()
            .map(|val| val.to_le_bytes().to_vec())
            .collect::<Vec<_>>()
            .concat();
        let calculated_output = PixelDataSupported::parse_lookup(test);
        let output: Vec<u8> = vec![2, 1, 1, 1, 4, 4, 1, 1, 4, 8, 2, 1, 1, 2, 2, 1, 1, 4, 4, 1];
        assert_eq!(output, calculated_output);
    }
}
