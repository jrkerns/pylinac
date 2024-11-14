# conftest.py

import json
import os
import threading

import pytest

# Define the path for the active tests file
ACTIVE_TESTS_FILE = "active_tests.json"

# Initialize a thread-safe set to store active tests
active_tests_lock = threading.Lock()
active_tests = set()


def update_active_tests_file():
    """Writes the current active tests to a JSON file."""
    with active_tests_lock:
        data = list(active_tests)
    with open(ACTIVE_TESTS_FILE, "w") as f:
        json.dump(data, f)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    """Hook to track test start and end, and trace memory allocations."""
    # Before the test runs
    with active_tests_lock:
        active_tests.add(item.nodeid)
    update_active_tests_file()

    try:
        # Run the actual test
        yield
    finally:
        # After the test runs
        with active_tests_lock:
            active_tests.discard(item.nodeid)
        update_active_tests_file()


@pytest.fixture(scope="session", autouse=True)
def ensure_active_tests_file_cleanup():
    """Ensure that the active tests file is removed after the test session."""
    yield
    if os.path.exists(ACTIVE_TESTS_FILE):
        os.remove(ACTIVE_TESTS_FILE)
