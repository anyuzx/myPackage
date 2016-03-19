from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

from sys import platform as _platform

if _platform == 'linux' or _platform == 'linux2':
	ext_modules = [Extension('core._pdist',['src/_pdist.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
				   Extension('core._cdist',['src/_cdist.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
	               Extension('analysis._contactmap',['src/_contactmap.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
	               Extension('analysis._contactprob',['src/_contactprob.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
	               Extension('analysis._subchaindist',['src/_subchaindist.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
	               Extension('analysis._hbond',['src/_hbond.pyx'],
							  include_dirs = [numpy.get_include()],
							  libraries=['m']),
	               Extension('create._lattice_chain',
	                          sources = ['src/_lattice_chain.pyx','src/c_lattice_chain.cpp'],
	                          extra_compile_args=['-std=c++11'],
	                          language='c++',
	                          libraries = ['m'],
	                          include_dirs = [numpy.get_include()])]
elif _platform == 'darwin':
	ext_modules = [Extension('core._pdist',['src/_pdist.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('core._cdist',['src/_cdist.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('analysis._contactmap',['src/_contactmap.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('analysis._contactprob',['src/_contactprob.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('analysis._subchaindist',['src/_subchaindist.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('analysis._hbond',['src/_hbond.pyx'],
							  include_dirs = [numpy.get_include()]),
	               Extension('create._lattice_chain',
	               	          sources = ['src/_lattice_chain.pyx','src/c_lattice_chain.cpp'],
	               	          extra_compile_args =['-std=c++11','-stdlib=libc++'],
	               	          language='c++',
	               	          include_dirs = [numpy.get_include()])]


setup(
	name = 'myPackage',
	version = '0.0.6',
	author = 'Guang Shi',
	author_email = 'stefanshi1988@gmail.com',
	description = ('A collection of tools I used in the research'),
	packages = ['core','analysis','create'],
	ext_modules = cythonize(ext_modules),
)
