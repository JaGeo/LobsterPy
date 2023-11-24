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

sys.path.insert(0, os.path.abspath("../../"))


# -- Project information -----------------------------------------------------


project = "LobsterPy"
copyright = "2022-2023, LobsterPy Development Team"
author = "Janine George"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinxarg.ext",
    "myst_nb",
    "sphinx_design",
    "sphinx_copybutton",
]

# napoleon_include_private_with_doc = True
# napoleon_include_special_with_doc = True
# autoclass_content = "both"


source_suffix = [".rst", ".md", ".ipynb"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

exclude_patterns = [
    "test*.py",
    "test",
    "Thumbs.db",
    ".DS_Store",
]

myst_heading_anchors = 2  # enable headings as link targets
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "dollarmath",
    "html_admonition",
    "html_image",
]

nb_execution_timeout = 500


# use type hints
autodoc_typehints = "description"
autoclass_content = "class"
autodoc_member_order = "bysource"

# better napoleon support
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_logo = "_static/logo.png"

html_theme_options = {
    "repository_provider": "github",
    "repository_url": "https://github.com/JaGeo/LobsterPy",
    "use_repository_button": True,
    "use_issues_button": True,
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = ['custom.css']

# hide sphinx footer
html_show_sphinx = False
html_show_sourcelink = False
html_title = "lobsterpy"
