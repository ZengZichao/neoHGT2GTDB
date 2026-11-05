"""
////////////////////////////////////////////////////////////
│	
│	Filename:	__init__.py
│	Comment:	Installation script
│	
│	=======================================================
│	Authorship:	Zhihui Li (Qiyun Zhu and Katharina Dittmar)
│	Copyright:	Distributed under the terms of the Modified BSD License.
│	
////////////////////////////////////////////////////////////
"""
#!/usr/bin/env python3

# the import package and the CLI entry point keep the upstream spelling
# `neoHGT` (README, "Naming, in one place"); the repository and the pip
# distribution are named neoHGT2GTDB (see setup.py).
__name__ = 'neoHGT'
# PEP 440: upstream base 2.0b3 + a local version segment marking this fork
# ('2.0b3.gtdb1' was rejected by modern setuptools and broke `pip install .`)
__version__ = '2.0b3+gtdb1'
__description__ = (
    '🧬 Unofficial GTDB adaptation of neoHGT: improved genome-wide detection of horizontal gene transfer (HGT) events based on sequence homology search hit distribution statistics')
__license__ = 'BSD-3-Clause'
__author__ = 'Zhihui Li and Qiyun Zhu'
# maintainer contact for THIS adaptation (曾子超 / Zeng Zichao, ORCID
# 0000-0001-6553-970X); the upstream authors' personal address is not
# published for an unofficial derivative
__email__ = 'zengzichao@sjtu.edu.cn'
__url__ = 'https://github.com/cgneo/neoHGT'
# homepage of the upstream project, kept for attribution via setup.py
# project_urls only -- this derivative claims no homepage of its own yet
__url__ = 'https://github.com/cgneo/neoHGT'
