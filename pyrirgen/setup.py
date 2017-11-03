from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os.path
import numpy as np

# python3 setup.py build_ext --inplace

extensions = [
	Extension("*", 
        sources=["pyrirgen.pyx"],
		include_dirs=[os.path.abspath(".")],
		language='c++'
	),
]

setup(
	name = 'pyrirgen',
    author = "ty274", 
    description = "Cython-based image method implementation for room impulse response generation",
    packages = setuptools.find_packages(), 

	ext_modules=cythonize(extensions, compiler_directives = {
		'language_level': 3, # Build for Python 3
		'embedsignature': True,
		'c_string_encoding': 'default',
		'c_string_type': 'str',
	}),
)
