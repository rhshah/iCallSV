from __future__ import print_function
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
import io
import codecs
import os
import sys
import iCallSV

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst')

setup(
    name='iCallSV',
    version=iCallSV.__version__,
    description='The module helps to call structural variants using NGS data set on human.',
    long_description=long_description,
    include_package_data=True,
    url='https://github.com/rhshah/iCallSV',
    download_url='https://github.com/rhshah/iCallSV/tarball/0.0.7',
    author=iCallSV.__author__,
    author_email='rons.shah@gmail.com',
    license=iCallSV.__license__,
    platforms='any',
    packages=['iCallSV'],
    install_requires=[
                      'pandas==0.16.2',
                      'pysam==0.8.4',
                      'pyvcf==0.6.7',
                      'biopython==1.65',
                      'coloredlogs==5.2'
                      ],
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    zip_safe=False,
    classifiers=(
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 5 - Production/Stable'
    ),
)
