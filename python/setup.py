"""Builds the Python bindings against the SAME C++ core the R package uses.

Nothing here is a reimplementation: ../src/core/mash.cpp is compiled directly,
so the two languages cannot drift apart. The core includes no R headers, which
is what makes this possible.
"""
NOTE = """macOS: build with Apple's clang, not conda's. An old ld64 produces a
shared object that fails at import with 'mis-aligned LINKEDIT string pool',
and the build itself reports success. See README.md."""

from pathlib import Path

import os
import re
import subprocess

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


def eigen_include():
    """Where Eigen lives.

    The core's solver uses Eigen, which the R package gets from RcppEigen.
    A Python build needs it separately: EIGEN_INCLUDE_DIR if set, otherwise
    pkg-config, otherwise the usual places.
    """
    if os.environ.get("EIGEN_INCLUDE_DIR"):
        return os.environ["EIGEN_INCLUDE_DIR"]
    try:
        out = subprocess.run(["pkg-config", "--cflags-only-I", "eigen3"],
                             capture_output=True, text=True, check=True).stdout
        for token in out.split():
            if token.startswith("-I"):
                return token[2:]
    except Exception:
        pass
    for candidate in ("/opt/homebrew/include/eigen3", "/usr/local/include/eigen3",
                      "/usr/include/eigen3", "/opt/anaconda3/include/eigen3"):
        if os.path.isdir(candidate):
            return candidate
    raise SystemExit(
        "Eigen headers not found. Install it (brew install eigen, apt install "
        "libeigen3-dev, conda install eigen) or set EIGEN_INCLUDE_DIR.")


def check_eigen_version(include_dir):
    """Warn if Eigen is not 3.4.

    The version is not cosmetic: a fit built against Eigen 5.0.1 differs from
    the R package's by about 3e-4 in the random effects, because the majors
    vectorise and reduce in different orders and thirteen IRLS iterations
    amplify it. Against 3.4 the same fit is bit-identical. Set
    EIGEN_INCLUDE_DIR to a 3.4 tree if the default is something else.
    """
    # Eigen 3.4 keeps the version in Core/util/Macros.h; later trees moved it,
    # so try both rather than silently skipping the check.
    text = None
    for relative in (("Eigen", "src", "Core", "util", "Macros.h"),
                     ("Eigen", "Version"),
                     ("Eigen", "Core")):
        try:
            with open(os.path.join(include_dir, *relative)) as handle:
                candidate = handle.read()
        except OSError:
            continue
        if "EIGEN_WORLD_VERSION" in candidate:
            text = candidate
            break
    if text is None:
        return
    found = {}
    for part in ("WORLD", "MAJOR", "MINOR"):
        match = re.search(r"#define EIGEN_%s_VERSION\s+(\d+)" % part, text)
        if match:
            found[part] = int(match.group(1))
    if found.get("WORLD") == 3 and found.get("MAJOR") == 4:
        return
    version = "%s.%s.%s" % (found.get("WORLD", "?"), found.get("MAJOR", "?"),
                            found.get("MINOR", "?"))
    print("WARNING: building against Eigen %s. PACE reproduces the R package's "
          "fits bit-for-bit only with Eigen 3.4; expect differences around 3e-4 "
          "in the random effects. Set EIGEN_INCLUDE_DIR to an Eigen 3.4 tree."
          % version)

# setuptools rejects absolute paths in `sources`, so the core is referenced
# relative to this file. It still compiles the one shared copy; nothing is
# duplicated into the Python tree.
CORE = Path("..") / "src" / "core"
EIGEN = eigen_include()
check_eigen_version(EIGEN)

ext_modules = [
    Pybind11Extension(
        "pace._core",
        sources=[
            "src/pace_py.cpp",
            str(CORE / "mash.cpp"),
            str(CORE / "decomposition.cpp"),
            str(CORE / "thread_pool.cpp"),
            str(CORE / "stage_timer.cpp"),
            str(CORE / "neighbourhood.cpp"),
            str(CORE / "random_design.cpp"),
            str(CORE / "special.cpp"),
            str(CORE / "preprocess.cpp"),
            str(CORE / "statistics.cpp"),
            str(CORE / "sparse_product.cpp"),
            str(CORE / "count_stats.cpp"),
            str(CORE / "linear_predictor.cpp"),
            str(CORE / "hyperparameters.cpp"),
            str(CORE / "dispersion.cpp"),
            str(CORE / "gene_solve.cpp"),
            str(CORE / "irls_chunk.cpp"),
            str(CORE / "irls_driver.cpp"),
        ],
        include_dirs=[str(CORE), EIGEN],
        # The R package passes this through configure and the core depends on
        # it: fused multiply-add changes the arithmetic, and a fit must not
        # depend on which language called it. Without it the random effects
        # drift from R's by about 1e-4.
        extra_compile_args=["-ffp-contract=off"],
        cxx_std=17,
    )
]

setup(cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)
