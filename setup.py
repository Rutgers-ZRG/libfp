import os
import sys
import sysconfig
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import numpy
from ctypes.util import find_library

class CustomBuildExtCommand(build_ext):
    def run(self):
        self.include_dirs.append(numpy.get_include())
        build_ext.run(self)

def find_lapack():
    # Check for Intel MKL first
    mklroot = os.environ.get('MKLROOT')
    if mklroot:
        print(f"Found MKLROOT at {mklroot}")
        return ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core']

    # Fallback to other LAPACK implementations
    lapack_libs = ['openblas', 'lapack', 'mkl_rt']
    for lib in lapack_libs:
        if find_library(lib):
            return [lib]

    print("Warning: No LAPACK library found. Trying to build without explicit LAPACK linkage.")
    return []

# Workaround Python issue 21121
config_var = sysconfig.get_config_var("CFLAGS")
if (config_var is not None and
    "-Werror=declaration-after-statement" in config_var):
    os.environ['CFLAGS'] = config_var.replace(
        "-Werror=declaration-after-statement", "")

sources = ['libfp.c', 'fingerprint.c', 'rcov.c', 'apc.c']
source_dir = "src" if os.path.exists('src') else "../src"
include_dirs = [source_dir, numpy.get_include()]
sources = [f"{source_dir}/{s}" for s in sources]

# Detect LAPACK library
lapack_lib = find_lapack()

# Try to find LAPACK library directory
lapack_dir = []
mklroot = os.environ.get('MKLROOT')
if mklroot:
    lapack_dir = [f"{mklroot}/lib/intel64"]
    include_dirs.append(f"{mklroot}/include")
elif sys.platform == "darwin":
    lapack_dir = ['/usr/local/opt/openblas/lib', '/usr/local/opt/lapack/lib']
elif sys.platform.startswith("linux"):
    lapack_dir = ['/usr/lib', '/usr/local/lib']

extra_link_args = []
if mklroot:
    extra_link_args = [f'-Wl,-rpath,{mklroot}/lib/intel64']
elif sys.platform == "darwin":
    extra_link_args = ['-Wl,-rpath,/usr/local/opt/openblas/lib', '-Wl,-rpath,/usr/local/opt/lapack/lib']
elif sys.platform.startswith("linux"):
    extra_link_args = ['-Wl,-rpath,/usr/lib', '-Wl,-rpath,/usr/local/lib']

extension = Extension('libfp._libfp',
                      include_dirs=include_dirs,
                      library_dirs=lapack_dir,
                      libraries=lapack_lib,
                      sources=['fppy/_libfp.c'] + sources,
                      extra_compile_args=[],
                      extra_link_args=extra_link_args,
                      define_macros=[])

setup(
    packages=find_packages(where='fppy'),
    package_dir={'': 'fppy'},
    ext_modules=[extension],
    cmdclass={'build_ext': CustomBuildExtCommand},
)

if not lapack_lib:
    print("No LAPACK library was found. The installation might fail or the library might not work correctly.")
    print("Please ensure you have LAPACK, OpenBLAS, or MKL installed on your system.")
else:
    print(f"Using LAPACK libraries: {lapack_lib}")
    print(f"Library directories: {lapack_dir}")
