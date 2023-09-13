# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pyMOE'
copyright = '2023, J. Cunha, D. E. Aguiam'
author = 'J. Cunha, D. E. Aguiam'
release = 'v1.3'
import os               # line 13
import sys              # line 14
sys.path.insert(0, os.path.abspath('../pyMOE/'))
sys.path.insert(0, os.path.abspath('..')) 
sys.path.insert(0, os.path.abspath('../..'))  
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'rst2pdf.pdfbuilder']

#pdf_documents = [('index', u'rst2pdf', u'Sample rst2pdf doc', u'Your Name'),]

templates_path = ['_templates']
exclude_patterns = []

autoclass_content = 'both'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'   
html_static_path = ['_static']
