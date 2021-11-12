#!/usr/bin/env python3
"""
Installer
"""

from setuptools import setup
import codecs
import os
from sys import version_info as _vi

package_name = "pynsn"

install_requires = ["numpy>=1.6",
                    "scipy>=1.0",
                    "Pillow>=5.0"]

extras_require = {
    'gui':                ['PyQt5>=5.14'],
    'pygame':             ["pygame>=1.9"],
    'expyriment':         ["expyriment>=0.9"],
    'dataframe':          ["pandas>=1.3"],
    'svg':                ["svgwrite>=1.4"]
}

entry_points = {'console_scripts': ['pynsn-gui=pynsn.gui:start']}

packages = [package_name]
for subp in ["dot_array","image", "lib", "qt", "database", "sequence"]:
    packages.append("{}.{}".format(package_name, subp))

if _vi.major< 1:
    raise RuntimeError("{0} requires Python 3 or larger.".format(package_name))

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
        keywords = "", #ToDo
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

    try:
        import PyQt5
    except:
        print("Pleae note:\n"
          "Running the PyNSN GUI requires the installation 'PyQt5' ")
