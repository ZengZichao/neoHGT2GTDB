"""
////////////////////////////////////////////////////////////
│	
│	Filename:	setup.py
│	Comment:	Installation script for the neoHGT2GTDB distribution
│			(unofficial GTDB adaptation of neoHGT; provides the `neoHGT` command)
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the Modified BSD License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

from setuptools import setup, find_packages
from glob import glob

import neoHGT.__init__ as init


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
    # package and the command it installs stay `neoHGT` -- see the README's
    # "Naming, in one place" paragraph
    'name':             'neoHGT2GTDB',
    'version':          init.__version__,
    'license':          init.__license__,
    'description':      init.__description__,
    'long_description': open('README.md').read(),
    'long_description_content_type': 'text/markdown',
    'author':           init.__author__,
    # maintainer of this adaptation; `author` credits the upstream code
    # authorship, whose personal address is deliberately not published
    'maintainer':       '曾子超 (Zeng Zichao)',
    'maintainer_email': init.__email__,
    # unofficial derivative: no homepage of its own is claimed; the upstream
    # project is linked for attribution only
    'project_urls':     {'Upstream project (official neoHGT2)': init.__url__,
                         'Maintainer ORCID': 'https://orcid.org/0000-0001-6553-970X'},
    'install_requires': ['pyyaml',
                         'numpy',
                         'scipy',
                         'pandas',
                         'scikit-learn',
                         'matplotlib'],
    'classifiers':      classes.strip().split('\n    '),
    'python_requires':  '>=3.6',
    # entry_points used to be the orphan string '[console_scripts]' (upstream
    # inheritance), which setuptools ignored; the command comes from `scripts=`
    'scripts':          glob(f'scripts/{init.__name__}'),
    'package_data':     {init.__name__: ['config.yml']},
    'include_package_data': True
}

# the test tree ships fixture-independent code but depends on data files the
# wheel does not carry: installed dead tests only confuse users (B 12.3)
setup(**params, packages=find_packages(exclude=['neoHGT.tests*']))
