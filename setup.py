"""
setup script for seq_deposit
"""
# !/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

with open("README.md", "r", encoding="utf-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

setup(
    name="seq_deposit",
    version="0.1.0",
    description="a set of commands for ensuring that all sequences ordered are properly organized and documented",
    long_description=readme,
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/seq_deposit",
    packages=[
        "seq_deposit",
    ],
    package_dir={"seq_deposit": "seq_deposit"},
    py_modules=["seq_deposit / seq_deposit"],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="seq_deposit",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={
        "console_scripts": [
            "seq-deposit = seq_deposit.cli:cli",
        ]
    },
)
