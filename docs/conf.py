# Copyright (c) 2025 Chair for Design Automation, TUM
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Sphinx configuration file."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import pybtex.plugin
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.template import field, href

ROOT = Path(__file__).parent.parent.resolve()


try:
    # Python 3.8+
    from importlib.metadata import version as get_version
except ImportError:
    # For older Python versions
    from importlib_metadata import version as get_version

try:
    version = get_version("mqt.qecc")
except ModuleNotFoundError:
    msg = (
        "Package should be installed to produce documentation! "
        "Assuming a modern git archive was used for version discovery."
    )
    warnings.warn(msg, stacklevel=1)


# Filter git details from version
release = version.split("+")[0]

if TYPE_CHECKING:
    from pybtex.database import Entry
    from pybtex.richtext import HRef

# -- Project information -----------------------------------------------------
project = "QECC"
author = "Chair for Design Automation, Technical University of Munich"
language = "en"
project_copyright = "2025, Chair for Design Automation, Technical University of Munich"

master_doc = "index"

templates_path = ["_templates"]
html_css_files = ["custom.css"]


# -- General configuration ---------------------------------------------------
extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxext.opengraph",
    "sphinx.ext.viewcode",
    "sphinxcontrib.inkscapeconverter",
    "sphinxcontrib.bibtex",
    # "sphinx.ext.napoleon",
    # "sphinx.ext.autodoc",
    # "sphinx.ext.autosummary",
    # "sphinx.ext.mathjax",
    # "sphinx.ext.intersphinx",
    # "sphinx.ext.autosectionlabel",
    # "sphinx.ext.viewcode",
    # "sphinx.ext.githubpages",
    # "sphinxcontrib.bibtex",
    # "sphinx_copybutton",
    # "nbsphinx",
    # "sphinxext.opengraph",
    # "sphinx_autodoc_typehints",
]

source_suffix = [".rst", ".md"]

exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "**.jupyter_cache",
    "**jupyter_execute",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

pygments_style = "colorful"

add_module_names = False

modindex_common_prefix = ["mqt.qecc."]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "typing_extensions": ("https://typing-extensions.readthedocs.io/en/latest", None),
    "qiskit": ("https://qiskit.org/documentation", None),
    "mqt": ("https://mqt.readthedocs.io/en/latest", None),
    "core": ("https://mqt.readthedocs.io/projects/core/en/latest", None),
    "ddsim": ("https://mqt.readthedocs.io/projects/ddsim/en/latest", None),
    "qcec": ("https://mqt.readthedocs.io/projects/qcec/en/latest", None),
    "qmap": ("https://mqt.readthedocs.io/projects/qmap/en/latest", None),
    "syrec": ("https://mqt.readthedocs.io/projects/syrec/en/latest", None),
}

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "substitution",
    "deflist",
    "dollarmath",
]
myst_substitutions = {
    "version": version,
}
myst_heading_anchors = 3


# nbsphinx_execute = "auto"
# highlight_language = "python3"
# nbsphinx_execute_arguments = [
#     "--InlineBackend.figure_formats={'svg', 'pdf'}",
#     "--InlineBackend.rc=figure.dpi=200",
# ]
# nbsphinx_kernel_name = "python3"

# autosectionlabel_prefix_document = True

# exclude_patterns = [
#     "_build",
#     "build",
#     "**.ipynb_checkpoints",
#     "Thumbs.db",
#     ".DS_Store",
#     ".env",
# ]

# -- Options for {MyST}NB ----------------------------------------------------

nb_execution_mode = "cache"
nb_mime_priority_overrides = [
    # builder name, mime type, priority
    ("latex", "image/svg+xml", 15),
]


class CDAStyle(UnsrtStyle):
    """Custom style for including PDF links."""

    def format_url(self, _e: Entry) -> HRef:  # noqa: PLR6301
        """Format URL field as a link to the PDF."""
        url = field("url", raw=True)
        return href()[url, "[PDF]"]


pybtex.plugin.register_plugin("pybtex.style.formatting", "cda_style", CDAStyle)

bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "cda_style"

copybutton_prompt_text = r"(?:\(venv\) )?\$ "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"

modindex_common_prefix = ["mqt.qecc."]

autoapi_dirs = ["../src/mqt"]
autoapi_python_use_implicit_namespaces = True
autoapi_root = "api"
autoapi_add_toctree_entry = False
autoapi_ignore = [
    "*/**/_version.py",
]
autoapi_options = [
    "members",
    "show-inheritance",
    "show-module-summary",
]
autoapi_keep_files = True
add_module_names = False
toc_object_entries_show_parents = "hide"
python_use_unqualified_type_names = True
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# -- Options for HTML output -------------------------------------------------
html_theme = "furo"
html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "mqt_dark.png",
    "dark_logo": "mqt_light.png",
    "source_repository": "https://github.com/cda-tum/mqt-qecc/",
    "source_branch": "main",
    "source_directory": "docs/source",
    "navigation_with_keys": True,
}
