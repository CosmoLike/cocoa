# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CoCoA'
copyright = '2025, Tim Eifler, Elisabeth Krause, Vivian Miranda and the Cosmolike Org'
author = 'Tim Eifler, Elisabeth Krause, Vivian Miranda and the Cosmolike Org'
release = '2020'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']


extensions = [
    "myst_parser",
]

# Tell Sphinx what file to start with
master_doc = "index"

# Optional: make .md recognized
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}