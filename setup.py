from setuptools import setup, find_packages

# PyPI checklist:
# 1. Ensure version has incremented.
# 2. To make life easy, in PyCharm, go to Tools/Run setup.py task.../sdist with command line option "upload" to push to PyPI

# 2 alt. To make a wheel, run $ python setup.py sdist bdist_wheel [upload]

__version__ = '1.2.2'
__version_info__ = (1, 2, 2)


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    # include_package_data=True,
    package_data={'pylinac':
                      ['demo_files/cbct/*',
                       'demo_files/log_reader/*',
                       'demo_files/starshot/*',
                       'demo_files/picket_fence/*',
                       'demo_files/vmat/*',
                       'demo_files/flatsym/*',
                       'demo_files/winston_lutz/*',
                       'demo_files/planar_imaging/*']}, # http://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    zip_safe=False,  # allows users to view files in egg/distribution
    url='https://github.com/jrkerns/pylinac',
    keywords='medical physics AAPM TG142 quality assurance starshot cbct vmat dynalog trajectory log',
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    install_requires=["numpy >= 1.9",
                      "scipy >= 0.15",
                      "pydicom >= 0.9.9",
                      "matplotlib >= 1.3.1",
                      "Pillow >= 2.5"],
    test_suite='tests._test_all',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Image Recognition",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries"]
)
