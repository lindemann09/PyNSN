#!/usr/bin/env python
"""
Installer
"""

import codecs
import os
from sys import version_info as _vi

from setuptools import find_packages, setup

PACKAGE_NAME = "pynsn"

if _vi.major < 3 or _vi.minor < 8:
    raise RuntimeError(f"{PACKAGE_NAME} requires Python 3.8 or larger.")

install_requires = ["numpy>=1.20",
                    "scipy>=1.5",
                    "Pillow>=8.4"]

extras_require = {
    'svg':                ["svgwrite>=1.4"],
    'pygame':             ["pygame>=1.9"],
    'expyriment':         ["expyriment>=0.9"],
    'matplotlib':         ["matplotlib>=3.2"]
}

entry_points = {}

packages = find_packages(".")


def readme():
    directory = os.path.dirname(os.path.join(
        os.getcwd(), __file__, ))
    with codecs.open(
        os.path.join(directory, "README.md"),
        encoding="utf8",
        mode="r",
        errors="replace",
    ) as file:
        return file.read()


def get_version(package):
    """Get version number"""

    with open(os.path.join(package, "__init__.py"), encoding="utf-8") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("'")[1]
    return "None"


if __name__ == '__main__':
    setup(
        name=PACKAGE_NAME,
        version=get_version(PACKAGE_NAME),
        description='Creating Non-Symbolic Number Displays',
        author='Oliver Lindemann',
        author_email='lindemann@cognitive-psychology.eu',
        license='GNU GPLv3',
        url='https://github.com/lindemann09/PyNSN',
        packages=packages,
        include_package_data=True,
        setup_requires=[],
        install_requires=install_requires,
        entry_points=entry_points,
        extras_require=extras_require,
        keywords="",
        classifiers=[
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering"
        ],
        long_description=readme(),
        long_description_content_type='text/markdown'
    )
