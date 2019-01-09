# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
from distutils.extension import Extension

setup(
    name="msvlm",
    version="0.1",
    author="Francis Brochu",
    author_email="francis.brochu.2@ulaval.ca",
    description='Python and c++ implementations of the Virtual Lock Mass correction and alignment algorithms',
    license="MIT",
    url="https://github.com/francisbrochu/msvlm",
    packages=find_packages(), 
    requires=['h5py', 'numpy', 'scipy', 'bisect', 'heappq'],
    ext_modules= [Extension(
        'msvlm/msAlign/msAlign',
        libraries = ['boost_python'],
        extra_compile_args=['-std=c++11'],
        sources = ['msvlm/msAlign/msAlignForPy.cpp', 
            'msvlm/msAlign/ActiveSequence.cpp', 
            'msvlm/msAlign/Heap.cpp'])]
)
