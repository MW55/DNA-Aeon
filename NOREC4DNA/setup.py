# -*- coding: utf-8 -*-
import numpy
from setuptools import find_packages
from distutils.core import setup, Extension

try:
    long_description = open("README.md").read()
except IOError:
    long_description = ""

setup(
    name="norec4dna",
    version="0.1.1",
    description="NOREC4DNA - a Fountain Code based approach to DNA-Storage",
    author="Michael Schwarz",
    packages=find_packages(),
    install_requires=['reedsolo', 'pytest', 'numpy', 'matplotlib'],
    zip_safe=False,
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    ext_modules=[Extension('cdnarules', ['cdnarules.c'], include_dirs=[numpy.get_include()]), ],
    include_dirs=[numpy.get_include()]
)
