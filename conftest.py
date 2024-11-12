import json
import os
import threading

import pytest

ACTIVE_TESTS_FILE = "active_tests.json"

active_tests_lock = threading.Lock()
active_tests = set()


def update_active_tests_file():
    with active_tests_lock:
        data = list(active_tests)
    with open(ACTIVE_TESTS_FILE, "w") as f:
        json.dump(data, f)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    with active_tests_lock:
        active_tests.add(item.nodeid)
    update_active_tests_file()

    # outcome = yield

    with active_tests_lock:
        active_tests.discard(item.nodeid)
    update_active_tests_file()


@pytest.fixture(scope="session", autouse=True)
def ensure_active_tests_file_cleanup():
    yield
    if os.path.exists(ACTIVE_TESTS_FILE):
        os.remove(ACTIVE_TESTS_FILE)
