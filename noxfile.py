from __future__ import annotations

import os
from pathlib import Path

import nox
from nox import Session


def set_uv_project_env(session) -> dict[str, str]:
    """Set uv environment vars used for uv sync/run against the session virtualenv."""
    return {
        "UV_PROJECT_ENVIRONMENT": session.virtualenv.location,
        "UV_PYTHON": session.virtualenv.location,
        "UV_NO_SYNC": "true",
        "UV_NO_DEFAULT_GROUPS": "true",
    }


def sync_uv_groups(session, groups: str | list[str]) -> dict[str, str]:
    """Run ``uv sync`` for one or more dependency groups using session-local uv env vars."""
    uv_env = set_uv_project_env(session)
    group_list = [groups] if isinstance(groups, str) else list(groups)
    cmd = ["uv", "sync"]
    for group in group_list:
        cmd.extend(["--group", group])
    session.run_install(*cmd, env=uv_env)
    return uv_env


@nox.parametrize(
    "numpy", ["1.23.5", "1.24.4", "1.25.2", "1.26.4", "2.0.0", "2.1.0", "2.2.0"]
)
@nox.parametrize("scipy", ["1.13.1", "1.14.1", "1.15.0"])
@nox.session(
    python="3.10",
    reuse_venv=True,
    venv_backend="uv",
)
def run_basic_test_suite_310(session: Session, numpy: str, scipy: str):
    sync_uv_groups(session, "test")
    session.install(f"numpy~={numpy}", f"scipy~={scipy}")
    session.install("pip")
    session.run("pip", "freeze")
    session.run(
        "pytest",
        "-n5",
        "tests_basic/core",
    )


@nox.parametrize(
    "numpy",
    ["1.23.5", "1.24.4", "1.25.2", "1.26.4", "2.0.0", "2.1.0", "2.2.0"],
)
@nox.parametrize("scipy", ["1.13.1", "1.14.1", "1.15.0"])
@nox.session(
    python="3.11",
    reuse_venv=True,
    venv_backend="uv",
)
def run_basic_test_suite_311(session: Session, numpy: str, scipy: str):
    sync_uv_groups(session, "test")
    session.install(f"numpy~={numpy}", f"scipy~={scipy}")
    session.install("pip")
    session.run("pip", "freeze")
    session.run(
        "pytest",
        "-n5",
        "tests_basic/core",
    )


@nox.parametrize("numpy", ["1.26.4", "2.0.0", "2.1.0", "2.2.0"])
@nox.parametrize("scipy", ["1.13.1", "1.14.1", "1.15.0"])
@nox.session(
    python="3.12",
    reuse_venv=True,
    venv_backend="uv",
)
def run_basic_test_suite_312(session: Session, numpy: str, scipy: str):
    sync_uv_groups(session, "test")
    session.install(f"numpy~={numpy}", f"scipy~={scipy}")
    session.install("pip")
    session.run("pip", "freeze")
    session.run(
        "pytest",
        "-n5",
        "tests_basic/core",
    )


@nox.parametrize("numpy", ["2.1.0", "2.2.0"])
@nox.parametrize("scipy", ["1.14.1", "1.15.0"])
@nox.session(
    python="3.13",
    reuse_venv=True,
    venv_backend="uv",
)
def run_basic_test_suite_313(session: Session, numpy: str, scipy: str):
    sync_uv_groups(session, "test")
    session.install(f"numpy~={numpy}", f"scipy~={scipy}")
    session.install("pip")
    session.run("pip", "freeze")
    session.run(
        "pytest",
        "-n5",
        "tests_basic/core",
    )


@nox.session(reuse_venv=True, venv_backend="uv|virtualenv")
def serve_docs(session: Session):
    sync_uv_groups(session, "docs")
    session.run(
        "sphinx-autobuild",
        "docs/source",
        "docs/build",
        "--port",
        "8777",
        "--open-browser",
    )


@nox.session(reuse_venv=True, venv_backend="uv|virtualenv")
def build_docs(session: Session):
    """Build the docs; used in CI pipelines to test the build. Will always rebuild and will always fail if there are any warnings"""
    sync_uv_groups(session, "docs")
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


@nox.session(reuse_venv=True, venv_backend="uv|virtualenv")
def build_docs_pdf(session: Session):
    """Build the docs as PDF. Sphinx has PDF capabilities but it doesn't like SVGs, GIFs, latex, etc.
    Instead, we build .rst -> HTML via Sphinx and then HTML -> PDF via plutoprint."""
    sync_uv_groups(session, "docs")
    session.install("plutoprint")
    session.run("sphinx-build", "docs/source", "docs/build-pdf", "-b", "singlehtml")
    session.log("Single HTML built with Sphinx. Building PDF with plutoprint...")
    output_file = "pylinac-docs.pdf"
    session.run(
        "plutoprint",
        "docs/build-pdf/index.html",
        output_file,
        "--size=A4",
        "--title=Pylinac",
    )
    session.log(f"PDF built successfully at {output_file}")


@nox.session(reuse_venv=True, venv_backend="uv")
def build_wheel(session: Session):
    """Build the wheel and sdist"""
    sync_uv_groups(session, "dev")
    session.run("uv", "build")


@nox.session(python=False)
def update_dev_kraken(session: Session):
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
