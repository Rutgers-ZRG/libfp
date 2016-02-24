from distutils.core import setup, Extension
import numpy
include_dirs_numpy = [numpy.get_include()]

ext = Extension('fppy',
                include_dirs=['../src/'] + include_dirs_numpy,
                library_dirs=['/usr/lib/'],
                libraries=['lapack'],
                sources=['pyfp.c',
                         '../src/fplib.c',
                         '../src/fingerprint.c',
                         '../src/rcov.c',
                         '../src/apc.c'])


setup(name='fppy', version='0.1', description='This if fppy.',
      ext_modules=[ext])

