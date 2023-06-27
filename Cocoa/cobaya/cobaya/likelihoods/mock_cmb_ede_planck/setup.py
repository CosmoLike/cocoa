# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
from itertools import chain
import re


def find_version():
    init_file = open(path.join(path.dirname(__file__), 'cobaya_mock_cmb', '__init__.py')).read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", init_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def get_long_description():
    """Get the long description from the README file."""
    with open(path.join(path.abspath(path.dirname(__file__)), 'README.md'),
              encoding='utf-8') as f:
        return f.read()


setup(
    name='cobaya_mock_cmb',
    version=find_version(),
    description='Mock CMB likelihood for Cobaya',
    long_description=get_long_description(),
    url="https://github.com/misharash/cobaya_mock_cmb",
    project_urls={
        'Source': 'https://github.com/misharash/cobaya_mock_cmb',
        'Tracker': 'https://github.com/misharash/cobaya_mock_cmb/issues'
    },
    author="Michael Rashkovetskyi, Julian B. Muñoz, Daniel J. Eisenstein and Cora Dvorkin",
    zip_safe=False,
    classifiers=[
        'Operating System :: OS Independent',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    python_requires='>=3.6.1',
    keywords='montecarlo sampling MCMC cosmology',
    packages=find_packages(),
    install_requires=['cobaya>=3.0.4', 'numpy>=1.7.0'],
    extras_require={},
    package_data={
        'cobaya_mock_cmb': list(chain(['*/*.yaml', '*/*.bibtex']))},
    entry_points={},
)
