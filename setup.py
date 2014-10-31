from pip.req import parse_requirements
from setuptools import setup

from pylinac import __version__

# PyPI checklist:
# 1. Ensure version has incremented.
# 2. In PyCharm, go to Tools/Run setup.py task.../sdist with command line option "upload" (to push to PyPI)

setup(
    name='pylinac',
    version=__version__,
    packages=['pylinac/starshot', 'pylinac/vmatqa'],
    include_package_data=True,
    zip_safe=False,  # allows users to view data in egg/distribution
    url='https://github.com/jrkerns/pylinac',
    keywords='medical physics AAPM TG142 quality assurance starshot',
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    install_requires=[str(ir.req) for ir in parse_requirements('requirements.txt')],
    test_suite='tools.test_all',
    classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Developers",
          "Intended Audience :: End Users/Desktop",
          "Intended Audience :: Healthcare Industry",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Natural Language :: English",
          "Programming Language :: Python",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.3",
          "Programming Language :: Python :: 3.4",
          "Topic :: Scientific/Engineering :: Medical Science Apps.",
          "Topic :: Scientific/Engineering :: Image Recognition",
          "Topic :: Scientific/Engineering :: Physics",
          "Topic :: Software Development :: Libraries",
    ],

)

