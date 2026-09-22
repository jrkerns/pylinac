"""Print the project version from the source configured in pyproject.toml."""

import runpy
from pathlib import Path

import tomllib

project_root = Path(__file__).resolve().parents[1]
with (project_root / "pyproject.toml").open("rb") as project_file:
    version_path = tomllib.load(project_file)["tool"]["hatch"]["version"]["path"]

print(runpy.run_path(str(project_root / version_path))["__version__"])
