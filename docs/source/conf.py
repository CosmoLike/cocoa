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

extensions = [
    "myst_parser",
    "sphinxemoji.sphinxemoji"
]

templates_path = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

language = "en"


# -- Options for HTML output -------------------------------------------------
html_theme = 'alabaster'
html_static_path = ['_static']

show_authors = True

# Output file base name for HTML help builder.
htmlhelp_basename = "cocoadoc"


# GPT code for rendering GitHub [TIP] and [Note]

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.transforms import SphinxTransform

class GithubAdmonitionTransform(SphinxTransform):
  default_priority = 210  # run after parsing

  def apply(self):
    for node in self.document.traverse(nodes.literal_block):
      continue  # skip code
    for node in self.document.traverse(nodes.paragraph):
      text = node.astext().strip()
      if text.startswith("[!TIP]"):
        node.replace_self(nodes.admonition("", nodes.paragraph(text=text[6:].strip()), classes=["tip"]))
      elif text.startswith("[!NOTE]"):
        node.replace_self(nodes.admonition("", nodes.paragraph(text=text[7:].strip()), classes=["note"]))

def setup(app):
  app.add_transform(GithubAdmonitionTransform)