import os
import sys
import sysconfig
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import numpy

class CustomBuildExtCommand(build_ext):
    def run(self):
        self.include_dirs.append(numpy.get_include())
        build_ext.run(self)

# Workaround Python issue 21121
config_var = sysconfig.get_config_var("CFLAGS")
if (config_var is not None and
    "-Werror=declaration-after-statement" in config_var):
    os.environ['CFLAGS'] = config_var.replace(
        "-Werror=declaration-after-statement", "")

sources = ['fplib.c', 'fingerprint.c', 'rcov.c', 'apc.c']
source_dir = "src" if os.path.exists('src') else "../src"
include_dirs = [source_dir, ]
sources = [f"{source_dir}/{s}" for s in sources]

lapack_dir = []
lapack_lib = ['lapack']
extra_link_args = []

extension = Extension('fplib._fplib',
                      include_dirs=include_dirs,
                      library_dirs=lapack_dir,
                      libraries=lapack_lib,
                      sources=['fppy/_fplib.c'] + sources,
                      extra_compile_args=[],
                      extra_link_args=extra_link_args,
                      define_macros=[])

setup(
    packages=find_packages(where='fppy'),
    package_dir={'': 'fppy'},
    ext_modules=[extension],
    cmdclass={'build_ext': CustomBuildExtCommand},
)
