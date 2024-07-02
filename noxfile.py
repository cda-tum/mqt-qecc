"""Nox sessions."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from typing import TYPE_CHECKING

import nox

if TYPE_CHECKING:
    from collections.abc import Sequence

nox.needs_version = ">=2024.3.2"
nox.options.default_venv_backend = "uv|virtualenv"

nox.options.sessions = ["lint", "tests"]

PYTHON_ALL_VERSIONS = ["3.8", "3.9", "3.10", "3.11", "3.12"]

# The following lists all the build requirements for building the package.
# Note that this includes transitive build dependencies of package dependencies,
# since we use `--no-build-isolation` to install the package in editable mode
# and get better caching performance. This only concerns dependencies that are
# not available via wheels on PyPI (i.e., only as source distributions).
BUILD_REQUIREMENTS = [
    "scikit-build-core[pyproject]>=0.8.1",
    "setuptools_scm>=7",
    "pybind11>=2.12",
    "wheel>=0.40",  # transitive dependency of pytest on Windows
    "Cython>=3; python_version > '3.11'",  # required to build ldpc on Python 3.12+
    "numpy>=1.26,<2; python_version > '3.11'",  # required to build ldpc on Python 3.12+
]

if os.environ.get("CI", None):
    nox.options.error_on_missing_interpreters = True


@nox.session(reuse_venv=True)
def lint(session: nox.Session) -> None:
    """Lint the Python part of the codebase using pre-commit.

    Simply execute `nox -rs lint` to run all configured hooks.
    """
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files", *session.posargs)


def _run_tests(
    session: nox.Session,
    *,
    install_args: Sequence[str] = (),
    run_args: Sequence[str] = (),
    extras: Sequence[str] = (),
) -> None:
    posargs = list(session.posargs)
    env = {"PIP_DISABLE_PIP_VERSION_CHECK": "1"}

    if os.environ.get("CI", None) and sys.platform == "win32":
        env["SKBUILD_CMAKE_ARGS"] = "-T ClangCL"

    if shutil.which("cmake") is None and shutil.which("cmake3") is None:
        session.install("cmake")
    if shutil.which("ninja") is None:
        session.install("ninja")

    _extras = ["test", *extras]
    if "--cov" in posargs:
        _extras.append("coverage")
        posargs.append("--cov-config=pyproject.toml")

    session.install(*BUILD_REQUIREMENTS, *install_args, env=env)
    install_arg = f"-ve.[{','.join(_extras)}]"
    session.install("--no-build-isolation", install_arg, *install_args, env=env)
    session.run("pytest", *run_args, *posargs, env=env)


@nox.session(reuse_venv=True, python=PYTHON_ALL_VERSIONS)
def tests(session: nox.Session) -> None:
    """Run the test suite."""
    _run_tests(session)


@nox.session(reuse_venv=True, venv_backend="uv", python=PYTHON_ALL_VERSIONS)
def minimums(session: nox.Session) -> None:
    """Test the minimum versions of dependencies."""
    _run_tests(
        session,
        install_args=["--resolution=lowest-direct"],
        run_args=["-Wdefault"],
    )
    session.run("uv", "pip", "list")


@nox.session(reuse_venv=True)
def docs(session: nox.Session) -> None:
    """Build the docs. Use "--non-interactive" to avoid serving. Pass "-b linkcheck" to check links."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", dest="builder", default="html", help="Build target (default: html)")
    args, posargs = parser.parse_known_args(session.posargs)

    serve = args.builder == "html" and session.interactive
    extra_installs = ["sphinx-autobuild"] if serve else []
    session.install(*BUILD_REQUIREMENTS, *extra_installs)
    session.install("--no-build-isolation", "-ve.[docs]")
    session.chdir("docs")

    if args.builder == "linkcheck":
        session.run("sphinx-build", "-b", "linkcheck", "source", "_build/linkcheck", *posargs)
        return

    shared_args = (
        "-n",  # nitpicky mode
        "-T",  # full tracebacks
        f"-b={args.builder}",
        "source",
        f"_build/{args.builder}",
        *posargs,
    )

    if serve:
        session.run("sphinx-autobuild", *shared_args)
    else:
        session.run("sphinx-build", "--keep-going", *shared_args)
