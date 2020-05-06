from setuptools import setup, find_packages

__version__ = '2.3.1'


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    package_data={'pylinac': ['watcher_config.yml', 'files/*.png']},
    url='https://github.com/jrkerns/pylinac',
    keywords="""medical physics AAPM TG142 quality assurance starshot cbct vmat dynalog starshot linac Varian Elekta
             trajectory log kv MV planar Leeds Las Vegas Standard Imaging PipsPro TG51""",
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    install_requires=["numpy >= 1.14",
                      "scipy >= 1.0",
                      "pydicom >= 1.0",
                      "matplotlib >= 1.4",
                      "scikit-image >= 0.13",
                      "scikit-learn >= 0.18",
                      "Pillow >= 4.0",
                      "tqdm >= 3.8",
                      "pyyaml >= 3.10",
                      "yagmail",
                      "argue",
                      "reportlab >= 3.3"],
    license='MIT',
    test_suite='test_basic',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Image Recognition",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries"]
)
