use pyo3::prelude::{Bound, PyModule, PyModuleMethods, PyResult, pymodule};
mod error;
mod reader;
use pyo3_stub_gen::define_stub_info_gatherer;

#[pymodule]
fn xim_reader(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let _ = m.add_class::<reader::XIMImage>();
    let _ = m.add_class::<reader::XIMHeader>();
    let _ = m.add_class::<reader::XIMHistogram>();
    let _ = m.add_class::<reader::XIMProperties>();
    Ok(())
}

define_stub_info_gatherer!(stub_info);
