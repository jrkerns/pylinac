from setuptools import setup, find_packages

# PyPI checklist:
# 1. Ensure version has incremented.
# 2. To make life easy, in PyCharm, go to Tools/Run setup.py task.../sdist with command line option "upload" to push to PyPI

__version__ = '0.1.3'
__version_info__ = (0, 1, 3)


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    # include_package_data=True,
    package_data={'pylinac': ['starshot/demo_files/*', 'vmatqa/demo_files/*', 'cbctqa/demo_files/*']},
    zip_safe=False,  # allows users to view files in egg/distribution
    url='https://github.com/jrkerns/pylinac',
    keywords='medical physics AAPM TG142 quality assurance starshot',
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    install_requires=["numpy >= 1.8",
                      "scipy >= 0.13",
                      "pydicom >= 0.9.8",
                      "matplotlib >= 1.3.1",
                      "future >= 0.13",
                      "Pillow >= 2.5",
                      "setuptools >= 5.8"],
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
