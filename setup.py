#! /usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="apogeebh",
    version='0.1.dev',
    author="Adrian Price-Whelan",
    author_email="adrn@princeton.edu",
    packages=["apogeebh", "apogeebh.db"],
    url="https://github.com/adrn/apogeebh",
    license="MIT",
    description="Search for black holes in APOGEE",
    long_description=""
)
