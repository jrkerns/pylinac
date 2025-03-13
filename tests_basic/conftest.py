import pytest
from tests_basic.utils import get_file_from_cloud_test_repo, get_folder_from_cloud_repo


@pytest.fixture(scope="session")
def anonymous_source_folder():
    return get_folder_from_cloud_repo(["mlc_logs", "_anonbase"])


@pytest.fixture(scope="session")
def anonymous_dest_folder():
    return get_folder_from_cloud_repo(["mlc_logs", "anonymous"])


@pytest.fixture(scope="session")
def test_files():
    return get_folder_from_cloud_repo(["test_files"])


@pytest.fixture(scope="session")
def bad_tif_path():
    return get_file_from_cloud_test_repo(["Winston-Lutz", "AQA_A_03082023.tif"])


@pytest.fixture(scope="session")
def tif_path():
    return get_file_from_cloud_test_repo(["Starshot", "Starshot-1.tif"])


@pytest.fixture(scope="session")
def png_path():
    return get_file_from_cloud_test_repo(["Starshot", "Starshot-1.png"])


@pytest.fixture(scope="session")
def dcm_path():
    return get_file_from_cloud_test_repo(["VMAT", "DRGSdmlc-105-example.dcm"])


@pytest.fixture(scope="session")
def as500_path():
    return get_file_from_cloud_test_repo(["picket_fence", "AS500#5.dcm"])


@pytest.fixture(scope="session")
def xim_path():
    return get_file_from_cloud_test_repo(["ximdcmtest.xim"])


@pytest.fixture(scope="session")
def xim_dcm_path():
    return get_file_from_cloud_test_repo(["ximdcmtest.dcm"])

@pytest.fixture(scope="session")
def rt_plan_file():
    return get_file_from_cloud_test_repo(["plan_generator", "Murray-plan.dcm"])

@pytest.fixture(scope="session")
def halcyon_plan_file():
    return get_file_from_cloud_test_repo(["plan_generator", "Halcyon Prox.dcm"])
