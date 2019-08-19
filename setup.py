#!/usr/bin/env python
from setuptools import setup,find_packages
import sys
from pathlib import Path

version="0.1.0"

if sys.version_info.major != 3:
    raise RuntimeError('MotifRaptor requires Python 3')


setup(name='MotifRaptor',
    version=version,
    description='TF binding effects on cell type specific genomic varations for GWAS',
    long_description=Path('README.md').read_text('utf-8'),
    url='https://github.com/pinellolab/MotifRaptor',
    author='Qiuming Yao',
    author_email='qyao1 AT mgh DOT harvard DOT edu',
    license='AGPL-3',
    packages=find_packages(include=['MotifRaptor','MotifRaptor.*']),
    package_dir={'MotifRaptor':'MotifRaptor'},
    package_data={'MotifRaptor': ['MotifRaptor/Database/*']},
    include_package_data = True,
    entry_points = {'console_scripts': ['MotifRaptor=MotifRaptor:main']},
          
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Cython',
        ],
    install_requires=[
        'argparse',
        'matplotlib',
        'numpy',
        'pandas',
        'pybedtools',
    	'seaborn==0.9.0',
        'scipy',
        'twobitreader'
        ]


)
