# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import pynsn
from pynsn import distributions

print("PyNSN Version: " + pynsn.__version__)

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'PyNSN'
copyright = '2022, Oliver Lindemann'
author = 'Oliver Lindemann'

# The full version, including alpha/beta/rc tags
release = pynsn.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions =["sphinx.ext.autodoc",  # Core Sphinx library for auto html doc generation from docstrings
             "sphinx.ext.autosummary",  # Create neat summary tables for modules/classes/methods etc
             #"sphinx.ext.intersphinx",  # Link to other project's documentation (see mapping below)
             "sphinx.ext.napoleon",
             "sphinx_autodoc_typehints",  # Automatically document param types (less noise in class signature)
             ]

# generate autosummary even if no references
autosummary_generate = False
autosummary_imported_members = True

autodoc_inherit_docstrings = True
autodoc_member_order = 'groupwise'
always_document_param_types =True
#typehints_fully_qualified=False
#autodoc_default_flags = ['members', 'undoc-members' ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'

#html_theme = 'sphinx_rtd_theme'
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "github_url": "https://github.com/lindemann09/PyNSN",
}

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

autodoc_mock_imports = ['bs4', 'requests']
