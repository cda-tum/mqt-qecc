import sys

if sys.version_info < (3, 10, 0):
    import importlib_metadata as metadata
else:
    from importlib import metadata

import pybtex.plugin
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.template import field, href

# -- Project information -----------------------------------------------------
project = "QECC"
author = "Lucas Berent"

release = metadata.version("mqt.qecc")
version = ".".join(release.split(".")[:3])
language = "en"
copyright = "Chair for Design Automation, Technical University of Munich"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "hoverxref.extension",
    "nbsphinx",
    "sphinxext.opengraph",
    "sphinx_rtd_dark_mode",
]

nbsphinx_execute = "auto"  # auto, never

highlight_language = "python3"

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'png2x'}",
    "--InlineBackend.rc=figure.dpi=96",
]

nbsphinx_kernel_name = "python3"

autosectionlabel_prefix_document = True

hoverxref_auto_ref = True
hoverxref_domains = ["cite", "py"]
hoverxref_roles = []
hoverxref_mathjax = True
hoverxref_role_types = {
    "ref": "tooltip",
    "p": "tooltip",
    "labelpar": "tooltip",
    "class": "tooltip",
    "meth": "tooltip",
    "func": "tooltip",
    "attr": "tooltip",
    "property": "tooltip",
}
exclude_patterns = ["_build", "build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store", ".env"]


class CDAStyle(UnsrtStyle):
    def format_url(self, e):
        url = field("url", raw=True)
        return href()[url, "[PDF]"]


pybtex.plugin.register_plugin("pybtex.style.formatting", "cda_style", CDAStyle)

bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "cda_style"

copybutton_prompt_text = r"(?:\(venv\) )?(?:\[.*\] )?\$ "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"

autosummary_generate = True

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_baseurl = "https://qecc.readthedocs.io/en/latest/"
html_logo = "_static/mqt_light.png"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
autodoc_member_order = "groupwise"
