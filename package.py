# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Phyex(CMakePackage, PythonExtension):
    """PHYEX: Physics parameterization package with Python bindings.

    PHYEX provides physics parameterizations for atmospheric models,
    including microphysics (ICE3/ICE4, LIMA), shallow convection,
    and turbulence schemes. This package includes Python bindings
    built with scikit-build-core and Cython.
    """

    homepage = "https://github.com/maurinl26/PHYEX"
    git = "https://github.com/maurinl26/PHYEX.git"
    url = "https://github.com/maurinl26/PHYEX/archive/refs/tags/v0.1.0.tar.gz"

    maintainers("maurinl26")

    license("Apache-2.0")

    version("main", branch="python-bindings")
    version("0.1.0", sha256="0000000000000000000000000000000000000000000000000000000000000000")

    # Build system variants
    variant("python", default=True, description="Build Python bindings")
    variant("double", default=True, description="Build with double precision")
    variant("single", default=False, description="Build with single precision")
    variant("progs", default=False, description="Build PHYEX program executables")

    # CMake and Fortran dependencies
    depends_on("cmake@3.15:", type="build")
    depends_on("ecbuild@3.11:", type="build")
    depends_on("fortran", type=("build", "link", "run"))

    # fiat dependency (fetched automatically via FetchContent)
    # Note: fiat is fetched from https://github.com/ecmwf-ifs/fiat during build

    # Python build dependencies (when +python)
    extends("python", when="+python")
    depends_on("python@3.11:", when="+python", type=("build", "link", "run"))
    depends_on("py-scikit-build-core@0.10:", when="+python", type="build")
    depends_on("py-cython", when="+python", type="build")
    depends_on("py-numpy@:1.999", when="+python", type=("build", "run"))

    # Python runtime dependencies
    depends_on("py-pip", when="+python", type="build")
    depends_on("py-setuptools", when="+python", type="build")

    # Optional development dependencies
    depends_on("py-pytest@7:", when="+python", type="test")

    # Conflicts
    conflicts("+double +single", msg="Cannot build both single and double precision simultaneously")
    conflicts("~double ~single", msg="Must enable either double or single precision")

    def cmake_args(self):
        """Configure CMake arguments for PHYEX build."""
        args = []

        # Precision options
        args.append(self.define_from_variant("ENABLE_DOUBLE_PRECISION", "double"))
        args.append(self.define_from_variant("ENABLE_SINGLE_PRECISION", "single"))
        args.append(self.define_from_variant("ENABLE_PHYEX_BUILD_PROGS", "progs"))

        # Shared libraries (required for Python bindings)
        if "+python" in self.spec:
            args.append(self.define("BUILD_SHARED_LIBS", True))

        # Fortran compiler
        args.append(self.define("CMAKE_Fortran_COMPILER", self.compiler.fc))

        return args

    @when("+python")
    def install(self, spec, prefix):
        """Install Python package using pip."""
        # First install the Fortran library with CMake
        super().install(spec, prefix)

        # Then install Python bindings
        with working_dir(self.stage.source_path):
            pip = which("pip")
            pip(
                "install",
                "--no-deps",
                "--prefix=" + prefix,
                "--no-build-isolation",
                "-v",
                ".",
            )

    @when("~python")
    def install(self, spec, prefix):
        """Install Fortran library only (no Python bindings)."""
        super().install(spec, prefix)

    def setup_build_environment(self, env):
        """Set up build environment."""
        # Ensure Fortran compiler is available
        env.set("FC", self.compiler.fc)
        env.set("F90", self.compiler.fc)

        if "+python" in self.spec:
            # scikit-build-core configuration
            env.set("SKBUILD_CMAKE_ARGS", ";".join(self.cmake_args()))

            # Set build type
            env.set("SKBUILD_BUILD_TYPE", "Release")

    def setup_run_environment(self, env):
        """Set up runtime environment."""
        if "+python" in self.spec:
            # Add Python package to PYTHONPATH
            python_version = self.spec["python"].version.up_to(2)
            env.prepend_path(
                "PYTHONPATH",
                join_path(
                    prefix.lib,
                    f"python{python_version}",
                    "site-packages"
                )
            )

    @run_after("install")
    @on_package_attributes(python=True)
    def install_test(self):
        """Run basic import test after installation."""
        python = which("python")
        python("-c", "import phyex; print(phyex.__version__)")
