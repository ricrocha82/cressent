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

project = "CRESSENT"
copyright = "2025, Ricardo R. Pavan"
author = "Ricardo R. Pavan"
html_title = "CRESSENT"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_design",
    "sphinx_copybutton",
    "myst_parser",
    'sphinx_rtd_theme',
]
source_suffix = [".rst", ".md"]
myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_favicon = "_static/figures/fig_cressent_new.png"
pygments_style = "tango"
pygments_dark_style = "tango"

html_theme_options = {
    # 📌 Sidebar layout
    "collapse_navigation": True,   # Collapse navigation bar entries
    "navigation_depth": 4,         # How many nested levels to show
    "titles_only": False,          # Only show page titles (no subsections)

    # 📌 Version control / menu
    # "display_version": True,       # Show version in the sidebar
    "sticky_navigation": True,     # Keeps sidebar scroll position
    "includehidden": True,         # Include hidden toctree entries

    # 📌 Branding
    "logo_only": False,            # If True, only show logo image (no project name)

    # 📌 Style tweaks (RTD theme uses CSS variables)
    # As of v1.3+, you can use theme-specific colors via CSS, but not many built‑in variables.
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = [
    'js/custom.js',
]