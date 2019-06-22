#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools.command.build_ext import build_ext as _build_ext

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    
# see https://stackoverflow.com/a/21621689/1862861 for why this is here
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

        

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open("requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

test_requirements = [
    "py",
    "pytest",
    "coverage"
]

setup_requirements = [
    'numpy',
    'setuptools_scm'
]

setup(
    name='gravpy',
    use_scm_version=True,
    description="A sandbox for gravitational wave astronomy.",
    long_description=readme + '\n\n' + history,
    author="Daniel Williams",
    author_email='daniel.williams@glasgow.ac.uk',
    url='https://github.com/transientlunatic/grasshopper',
    packages=[
        'gravpy',
    ],
    package_dir={'gravpy':
                 'gravpy'},
    include_package_data=True,
    setup_requires = setup_requirements,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='gravpy',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
