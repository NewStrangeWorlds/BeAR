# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BeAR'
copyright = '2020 - 2026, Daniel Kitzmann'
author = 'Daniel Kitzmann'
release = '2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = []
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_title = 'BeAR documentation'
html_static_path = ['static']
html_css_files = ['custom.css']
html_logo = 'static/bear_small.png'
html_favicon = 'static/bear_small.png'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation: https://pradyunsg.me/furo/customisation/
#
html_theme_options = {
    'sidebar_hide_name': False,
    'navigation_with_keys': True,
    'top_of_page_buttons': ['view'],
    'source_repository': 'https://github.com/daniel-kitzmann/beardevdoc/',
    'source_branch': 'main',
    'source_directory': 'docs_src/',
    'light_css_variables': {
        'color-brand-primary': '#8a5a2b',
        'color-brand-content': '#8a5a2b',
        'font-stack--monospace': 'SFMono-Regular, Menlo, Consolas, "Liberation Mono", monospace',
    },
    'dark_css_variables': {
        'color-brand-primary': '#d9a06b',
        'color-brand-content': '#d9a06b',
    },
    'footer_icons': [
        {
            'name': 'GitHub',
            'url': 'https://github.com/newstrangeworlds/bear',
            'html': '''
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            ''',
            'class': '',
        },
    ],
}
