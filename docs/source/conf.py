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

sys.path.insert(0, os.path.abspath("../"))


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
    "sphinx.ext.coverage",
    "sphinxarg.ext",
    "m2r2",
]

napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
autoclass_content = "both"

source_suffix = [".rst", ".md"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

exclude_patterns = [
    "../../lobsterpy/test",
    "../../lobsterpy/cohp/test",
    "../../lobsterpy/plotting/test",
    "../../lobsterpy/featurize/test",
    "../../lobsterpy/structuregraph/test",
    "../../lobsterpy/TestData",
    "Thumbs.db",
    ".DS_Store",
]


def run_apidoc(_):
    import subprocess
    import glob

    output_path = os.path.abspath(os.path.dirname(__file__))
    excludes = glob.glob(os.path.join(output_path, "../../lobsterpy/cohp/test"))
    excludes1 = glob.glob(os.path.join(output_path, "../../lobsterpy/test"))
    excludes2 = glob.glob(os.path.join(output_path, "../../lobsterpy/plotting/test"))
    excludes3 = glob.glob(os.path.join(output_path, "../../lobsterpy/featurize/test"))
    excludes4 = glob.glob(
        os.path.join(output_path, "../../lobsterpy/structuregraph/test")
    )
    module = os.path.join(output_path, "../../lobsterpy")
    cmd_path = "sphinx-apidoc"
    command = [
        cmd_path,
        "-e",
        "-o",
        output_path,
        module,
        " ".join(excludes),
        " ".join(excludes1),
        " ".join(excludes2),
        " ".join(excludes3),
        " ".join(excludes4),
        "--force",
    ]
    subprocess.check_call(command)


def setup(app):
    app.connect("builder-inited", run_apidoc)


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
