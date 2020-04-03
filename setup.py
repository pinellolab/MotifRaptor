#!/usr/bin/env python
from setuptools import setup,find_packages,Extension
import sys
from pathlib import Path
import os

version="0.3.0"

if sys.version_info.major != 3:
    raise RuntimeError('MotifRaptor requires Python 3')


from numpy import get_include as numpy_get_include
numpy_include_dir = [numpy_get_include()]

try:
    import Cython.Distutils
    has_cython = True
except:
    has_cython = False

ext = '.pyx' if has_cython else '.c'
if os.path.exists("MotifRaptor/SNPScanner/motif_matching_lcp.c"):
    ext = '.c'
ext_modules = [Extension("MotifRaptor.SNPScanner.motif_matching_lcp",["MotifRaptor/SNPScanner/motif_matching_lcp"+ext],include_dirs=numpy_include_dir,extra_compile_args=['-w']),]
if has_cython:
    from Cython.Build import cythonize
    ext_modules = cythonize(ext_modules, language_level=2)


try:
    __import__('imp').find_module('pyBigWig')
    # Make things with supposed existing module
except ImportError:
    os.system("pip install pybigwig")
    pass

print("---------Install Thrid Party ----------")
os.system("bash install.sh")


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
    entry_points = {'console_scripts': ['MotifRaptor=MotifRaptor.MotifRaptor:main']},
          
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
        'statsmodels',
        'twobitreader'
        ],
   
    ext_modules = ext_modules,
    cmdclass = {},

)
