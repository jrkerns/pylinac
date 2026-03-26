# AGENTS.md: Pylinac Test Suite Guide

## Purpose

Primary test suite for pylinac. All feature-level and core module tests live here. Tests use `unittest.TestCase` as the base class with `pytest` as the runner.

Refer to `../AGENTS.md` for the full codebase guide.

## Running Tests

```bash
# All tests
uv run pytest

# Single analysis module
uv run pytest tests_basic/test_cbct.py

# Core module tests only
uv run pytest tests_basic/core

# Contrib tests only
uv run pytest tests_basic/contrib

# With explicit parallel workers (--dist loadscope is the default)
uv run pytest tests_basic/test_cbct.py -n 3

# Nox compatibility matrix (core tests across NumPy/SciPy versions)
uv run nox -s run_basic_test_suite_312
```

## Directory Layout

| Path | Contents |
| --- | --- |
| `test_<module>.py` | One test file per analysis module (test_cbct, test_acr, test_starshot, test_winstonlutz, etc.) |
| `core/` | Tests for `pylinac/core/` modules (test_image, test_geometry, test_profile, test_gamma, test_mtf, etc.) |
| `contrib/` | Tests for `pylinac/contrib/` modules (test_quasar, test_orthogonality) |
| `test_files/` | Local test data cache -- GCP downloads are stored here and reused across runs |
| `utils.py` | Shared GCP download helpers and test mixins |
| `__init__.py` | Module-level config: `TEST_FILES_DIR`, `TEST_BANK_DIR`, `HIDE_PLOTS`, `DELETE_FILES` |

## Test Data Management

### GCP Bucket

Large test datasets live in the GCP bucket `pylinac_test_files`. Files are downloaded on first use and cached locally in `tests_basic/test_files/`.

- **Local credentials**: `GCP_creds.json` at the repo root.
- **CI credentials**: `GOOGLE_CREDENTIALS` environment variable (base64-encoded JSON).

### Download Functions (from `utils.py`)

**`get_file_from_cloud_test_repo(path: list[str], force: bool = False) -> str`**
Downloads a single file from GCP. Returns the local path. Skips download if the file already exists locally (unless `force=True`). Logs the MD5 hash on download.

**`get_folder_from_cloud_repo(folder: list[str], skip_exists: bool = True) -> str`**
Downloads an entire folder from GCP. Skips the network call entirely if the local folder exists and is non-empty (when `skip_exists=True`).

### Cache Behavior

- Downloaded files persist in `test_files/` between runs.
- CI workflows cache specific `test_files/` subdirectories per feature to avoid redundant downloads.
- Set `DELETE_FILES=true` (env var) to clean up downloaded files after each test class.

### GCP Path Convention

GCP paths mirror the local `test_files/` structure. For example:
- GCP: `CBCT/CatPhan_604/CBCTCatPhan604.zip`
- Local: `tests_basic/test_files/CBCT/CatPhan_604/CBCTCatPhan604.zip`

## Test Mixins

All mixins are defined in `utils.py`. Combine them with `unittest.TestCase` to get automatic constructor and loading tests.

### CloudFileMixin

Provides GCP-backed file resolution. Set `dir_path` and `file_name`, then call `get_filename()` to get the local path (downloading if needed).

```python
class MyTest(CloudFileMixin, TestCase):
    dir_path = ["Starshot"]
    file_name = "Starshot-1.tif"
```

### InitTesterMixin

Tests that the analysis class can be instantiated from a file or folder. Provides `test_init()`.

Required attributes: `klass`, `init_file` (list of path segments), `init_kwargs` (dict), `is_folder` (bool).

### FromZipTesterMixin

Tests `klass.from_zip()`. Provides `test_from_zip()`.

Required attributes: `klass`, `zip` (list of path segments), `zip_kwargs` (dict).

### FromDemoImageTesterMixin

Tests `klass.from_demo_image()` (or a custom method name). Provides `test_from_demo()`.

Required attributes: `klass`. Optional: `demo_load_method` (default `"from_demo_image"`).

### FromURLTesterMixin

Tests `klass.from_url()` against the demo files bucket. Provides `test_from_url()`.

Required attributes: `klass`, `url` (relative path within the demo bucket).

### PlotlyTestMixin

Tests `plotly_analyzed_images()` output. Provides `test_plotly_render()` and `test_plotly_options()`.

Required attributes: `instance` (analyzed object), `num_figs` (expected figure count). Optional: `fig_data` (dict mapping figure index to expected title, trace count, axis labels).

### DataBankMixin

Bulk-processes an entire directory of files. Used primarily in `tests_bank/`, not `tests_basic/`. Requires a local `pylinac test files` directory.

## Test Class Conventions

### Naming

- **Files**: `test_<module>.py` matching the source module name.
- **Classes**: `Test<Feature>` (e.g. `TestACRCT`, `TestStarshotLoading`, `CatPhan604Test`).
- **Methods**: `test_<behavior>` (e.g. `test_from_multiples`, `test_results_data`, `test_hu_values`).

### Structure

```python
class TestACRCT(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = ACRCT
    zip = [*TEST_DIR_CT, "Philips.zip"]
    init_file = [*TEST_DIR_CT, "Philips_CT"]
    is_folder = True

    def test_phantom_center(self):
        self.assertAlmostEqual(self.instance.phantom_center.x, 256, delta=3)
```

### Assertions

- **unittest**: `self.assertEqual`, `self.assertAlmostEqual`, `self.assertIsInstance`, `self.assertRaises`
- **NumPy**: `numpy.testing.assert_array_almost_equal`
- **Parameterized**: `@parameterized.expand` from the `parameterized` package

### Module-Level Data

Fetch shared test files at module level so they download once per test session:

```python
dcm_path = get_file_from_cloud_test_repo(["VMAT", "DRGSdmlc-105-example.dcm"])


class TestSomething(TestCase):
    def test_load(self):
        img = image.load(dcm_path)
```

## Markers

| Marker | Purpose |
| --- | --- |
| `@pytest.mark.catphan604` | CatPhan 604 tests. Used in CI for cache partitioning. |

## conftest.py (repo root)

The root `conftest.py` provides:

- **Active test tracking**: Writes currently-running test node IDs to `active_tests.json` via the `pytest_runtest_protocol` hook. Used for memory monitoring in CI.
- **Session cleanup**: Removes `active_tests.json` after the test session completes.

## Adding New Tests

1. Create `test_<module>.py` in `tests_basic/` (or the appropriate subdirectory for `core/` or `contrib/` modules).
2. Use the standard mixins for constructor and loading tests (`InitTesterMixin`, `FromZipTesterMixin`, `FromDemoImageTesterMixin`, `FromURLTesterMixin`).
3. Upload test data to the GCP bucket `pylinac_test_files` using the same directory structure that `test_files/` expects.
4. Add CI cache entries in `.github/workflows/ci.yml` if the test data is large enough to warrant caching.
5. Add a `@pytest.mark.<name>` marker if the tests need special CI handling, and register it in `pyproject.toml`.

## Maintenance Protocol

Update this file whenever:

- New test mixins or base classes are added to `utils.py`.
- The test data strategy changes (new buckets, new download patterns).
- New markers are introduced.
- The pytest configuration in `pyproject.toml` changes.
