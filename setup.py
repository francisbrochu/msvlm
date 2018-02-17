# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name="msvlm",
    version="0.1",
    author="Francis Brochu",
    author_email="francis.brochu.2@ulaval.ca",
    description='Python and c++ implementations of the Virtual Lock Mass correction and alignment algorithms',
    license="MIT",
    url="",
    packages=find_packages(), requires=['h5py', 'numpy', 'scipy', 'bisect', 'heappq']
)
