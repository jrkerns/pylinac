import os
from pathlib import Path

import nox


@nox.session(python=False)
@nox.parametrize("py_ver", ["3.9", "3.10", "3.11", "3.12"])
def run_tests(session, py_ver):
    py_str = f".venv{py_ver}"
    session.run("uv", "python", "pin", py_ver)
    session.run(
        "uv",
        "sync",
        "--quiet",
        env={"UV_PROJECT_ENVIRONMENT": py_str, "VIRTUAL_ENV": py_str},
    )
    session.run(
        "uv",
        "run",
        "pytest",
        env={"UV_PROJECT_ENVIRONMENT": py_str, "VIRTUAL_ENV": py_str},
    )


@nox.session(python=False)
def serve_docs(session):
    session.run(
        "sphinx-autobuild",
        "docs/source",
        "docs/build",
        "--port",
        "8777",
        "--open-browser",
    )


@nox.session(python=False)
def build_docs(session):
    """Build the docs; used in CI pipelines to test the build. Will always rebuild and will always fail if there are any warnings"""
    session.run(
        "sphinx-build",
        "docs/source",
        "docs/build",
        "-W",
        "--keep-going",
        "-a",
        "-E",
        "-b",
        "html",
        "-q",
    )


@nox.session(python=False)
def update_dev_kraken(session):
    """Run the Kraken build to update it with new pylinac changes"""
    if Path("GCP_creds.json").exists():
        os.environ["GCP_BUILD_CREDS"] = Path("gcp_build_creds.json").open().read()
    key_info = os.environ["GCP_BUILD_CREDS"]
    with open("service_key.json", "w") as key_file:
        key_file.write(key_info)
    session.run(
        "gcloud", "auth", "activate-service-account", "--key-file", "service_key.json"
    )
    session.run("gcloud", "config", "set", "project", "radmachine")
    session.run("gcloud", "builds", "triggers", "run", "Kraken", "--branch=master")
