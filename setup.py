#!/usr/bin/env python3
"""
Installer
"""

from setuptools import setup, find_packages
import codecs
import os
from sys import version_info as _vi

package_name = "pynsn"

if _vi.major < 3 or _vi.minor < 6:
    raise RuntimeError("{0} requires Python 3.6 or larger.".format(package_name))

install_requires = ["numpy>=1.19",
                    "scipy>=1.5",
                    "Pillow>=8.4",
                    "svgwrite>=1.4"]

extras_require = {
    'gui':                ["PyQt5>=5.14"],
    'pygame':             ["pygame>=1.9"],
    'expyriment':         ["expyriment>=0.9"],
    'matplotlib':         ["matplotlib>=3.2"]
}

entry_points = {'console_scripts': ['pynsn=pynsn.gui:start']}

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

    with open(os.path.join(package, "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("'")[1]
    return "None"


if __name__ == '__main__':
    setup(
        name = package_name,
        version=get_version(package_name),
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
        keywords = "",
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
