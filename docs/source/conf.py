from __future__ import annotations

import os
import sys
from importlib.metadata import distribution

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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "Arrakis"
copyright = "2021, Alec Thomson"
author = "Alec Thomson"

# The full version, including alpha/beta/rc tags
release = distribution("arrakis").version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "numpydoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "myst_parser",
    "autoapi.extension",
    "sphinx.ext.autosectionlabel",
    "sphinxcontrib.mermaid",
]

source_suffix = [".rst"]

napoleon_google_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

autoapi_type = "python"
autoapi_dirs = ["../../arrakis", "../../scripts"]
autoapi_member_order = "groupwise"
autoapi_keep_files = False
autoapi_root = "autoapi"
autoapi_add_toctree_entry = True
# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
html_title = "arrakis"
html_logo = "SPICE-RACS_circ.png"
logo_only = True
html_favicon = "favicon.ico"

# autodoc_member_order = 'bysource'

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("../"))

mermaid_params = [
    "--theme",
    "forest",
    "--width",
    "600",
    "--backgroundColor",
    "transparent",
]
