import os
import sys

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

project = 'SPICE-RACS'
copyright = '2021, Alec Thomson'
author = 'Alec Thomson'

# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram',
    # 'myst_parser',
    'autoapi.extension',
    'm2r2',
]

source_suffix = ['.rst']

napoleon_google_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

autoapi_type = 'python'
autoapi_dirs = ['../../spiceracs', '../../scripts']
# autoapi_dirs = ['../../spiceracs']
autoapi_member_order = 'groupwise'
autoapi_keep_files = False
# autoapi_root = 'api'
autoapi_template_dir = '_autoapi_templates'
autoapi_add_toctree_entry = True
# autoapi_generate_api_docs = True
autoapi_generate_api_docs = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_theme_options = {

    "repository_url": "https://bitbucket.csiro.au/projects/SPICE/repos/spiceracs",
    "use_repository_button": True,
    # "use_issues_button": True,
    # "use_edit_page_button": True,
    "logo_only": True,
    "show_navbar_depth":1,

}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
html_title = "SPICE-RACS"
html_logo = "SPICE-RACS_circ.png"
logo_only=True
html_favicon="favicon.ico"

# autodoc_member_order = 'bysource'

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../spiceracs'))
sys.path.insert(0, os.path.abspath('../askap_surveys'))
sys.path.insert(0, os.path.abspath('../askap_surveys/racs'))
sys.path.insert(0, os.path.abspath('../rmtable'))