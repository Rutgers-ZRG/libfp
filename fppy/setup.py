import os
import sys
import sysconfig


setup_type = sys.argv[1]

lapack_dir=['/opt/local/lib']
lapack_lib=['openblas']

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext
    use_setuptools = True
    print("setuptools is used.")

    class CustomBuildExtCommand(build_ext):
        def run(self):
            import numpy
            self.include_dirs.append(numpy.get_include())
            build_ext.run(self)

except ImportError:
    from distutils.core import setup, Extension
    use_setuptools = False
    print("distutils is used.")

    try:
        from numpy.distutils.misc_util import get_numpy_include_dirs
    except ImportError:
        print("Please install numpy first before installing fplib...")
        sys.exit(1)

# Workaround Python issue 21121
config_var = sysconfig.get_config_var("CFLAGS")
if (config_var is not None and
    "-Werror=declaration-after-statement" in config_var):
    os.environ['CFLAGS'] = config_var.replace(
        "-Werror=declaration-after-statement", "")

sources = ['fplib.c',
           'fingerprint.c',
           'rcov.c',
           'apc.c']

if os.path.exists('src'):
    source_dir = "src"
else:
    source_dir = "../src"

include_dirs = [source_dir, ]
if not use_setuptools:
    include_dirs += get_numpy_include_dirs()

for i, s in enumerate(sources):
    sources[i] = "%s/%s" % (source_dir, s)

extra_compile_args = []
extra_link_args = []
define_macros = []

extension = Extension('fplib._fplib',
                      include_dirs=include_dirs,
                      library_dirs=lapack_dir,
                      libraries=lapack_lib,
                      sources=['_fplib.c'] + sources,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args,
                      define_macros=define_macros)

version_nums = [None, None, None]
with open("%s/version.h" % source_dir) as w:
    for line in w:
        for i, chars in enumerate(("MAJOR", "MINOR", "MICRO")):
            if chars in line:
                version_nums[i] = int(line.split()[2])

# To deploy to pypi by travis-CI
nanoversion = 0
if os.path.isfile("__nanoversion__.txt"):
    with open('__nanoversion__.txt') as nv:
        try:
            for line in nv:
                nanoversion = int(line.strip())
                break
        except ValueError:
            nanoversion = 0
if nanoversion != 0:
    version_nums.append(nanoversion)

if None in version_nums:
    print("Failed to get version number in setup.py.")
    raise

version = ".".join(["%d" % n for n in version_nums[:3]])
if len(version_nums) > 3:
    version += "-%d" % version_nums[3]
if use_setuptools:
    setup(name='fplib',
          version=version,
          cmdclass={'build_ext': CustomBuildExtCommand},
          setup_requires=['numpy', 'setuptools>=18.0'],
          license='MIT',
          description='This is the fplib module.',
          author='Li Zhu',
          author_email='z@zhuli.name',
          url='https://github.com/zhuligs/fplib',
          packages=['fplib'],
          install_requires=['numpy', ],
          provides=['fplib'],
          platforms=['all'],
          ext_modules=[extension])
else:
    setup(name='fplib',
          version=version,
          license='MIT',
          description='This is the fplib module.',
          author='Li Zhu',
          author_email='z@zhuli.name',
          url='https://github.com/zhuligs/fplib',
          packages=['fplib'],
          requires=['numpy'],
          provides=['fplib'],
          platforms=['all'],
          ext_modules=[extension])
