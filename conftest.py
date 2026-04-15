# conftest.py

import json
import os
import threading
import time

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

    temp_file = f"{ACTIVE_TESTS_FILE}.{os.getpid()}.tmp"
    # This file is best-effort telemetry for memory monitoring; test execution
    # should never fail due to lock contention from xdist workers.
    for attempt in range(5):
        try:
            with open(temp_file, "w") as f:
                json.dump(data, f)
            os.replace(temp_file, ACTIVE_TESTS_FILE)
            return
        except PermissionError:
            time.sleep(0.01 * (2**attempt))
        finally:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except OSError:
                    pass


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
        try:
            os.remove(ACTIVE_TESTS_FILE)
        except OSError:
            pass
