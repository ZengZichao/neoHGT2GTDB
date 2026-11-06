"""
////////////////////////////////////////////////////////////
│	
│	Filename:	setup.py
│	Comment:	Installation script for the neoHGT2GTDB distribution
│			(unofficial GTDB adaptation of neoHGT; provides the `neoHGT` command)
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the BSD 3-Clause License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

import re
from setuptools import setup, find_packages
from glob import glob

import neoHGT.__init__ as init


REPO_URL = 'https://github.com/ZengZichao/neoHGT2GTDB'
BLOB_URL = f'{REPO_URL}/blob/main'


def long_description():
    """README text with relative links rewritten for the PyPI page.

    PyPI renders the README outside the repository, so a link such as
    `docs/GTDB-数据库.md` would resolve against pypi.org. Rewriting the few
    repository-relative targets to absolute GitHub URLs keeps the package page
    navigable; the repository copy keeps using relative links.
    """
    with open('README.md', encoding='utf-8') as fh:
        text = fh.read()
    targets = (
        'README.zh-CN.md', 'CHANGELOG.md', 'CONTRIBUTING.md', 'LICENSE',
    )
    for name in targets:
        text = text.replace(f']({name})', f']({BLOB_URL}/{name})')
    text = re.sub(r'\]\((docs/[^)]+)\)', rf']({BLOB_URL}/\1)', text)
    return text


classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS
    Operating System :: Microsoft :: Windows
"""

params = {
    # the distribution carries this repository's name (neoHGT2GTDB); the import
    # package and the command it installs stay `neoHGT`
    'name':             'neoHGT2GTDB',
    'version':          init.__version__,
    'license':          init.__license__,
    'description':      init.__description__,
    'long_description': long_description(),
    'long_description_content_type': 'text/markdown',
    'author':           init.__author__,
    'maintainer':       '曾子超 (Zeng Zichao)',
    'maintainer_email': init.__email__,
    'url':              REPO_URL,
    'project_urls':     {
        'Source': REPO_URL,
        'Documentation': f'{REPO_URL}#documentation',
        'Changelog': f'{BLOB_URL}/CHANGELOG.md',
        'Upstream project (official neoHGT2)': init.__url__,
        'Maintainer ORCID': 'https://orcid.org/0000-0001-6553-970X',
    },
    'install_requires': ['pyyaml',
                         'numpy',
                         'scipy',
                         'pandas',
                         'scikit-learn',
                         'matplotlib'],
    'classifiers':      classes.strip().split('\n    '),
    'python_requires':  '>=3.6',
    'scripts':          glob(f'scripts/{init.__name__}'),
    'package_data':     {init.__name__: ['config.yml']},
    'include_package_data': True
}

setup(**params, packages=find_packages(exclude=['neoHGT.tests*']))
