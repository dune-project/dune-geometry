import sys, os
from setuptools import find_packages
from skbuild import setup

from dune.common.dunepackaging import metaData

data, cmake_flags = metaData()
requires=(data.install_requires+data.dune_dependencies).__str__()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name=data.name,
    version=data.version,
    author=data.author,
    author_email=data.author_email,
    description=data.description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=data.url if data.url is not None else '',
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    install_requires=requires,
    classifiers=[
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
    ],
    python_requires='>=3.4',
    cmake_args=cmake_flags
)
