from setuptools import setup, find_packages

__version__ = '2.0.2'


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    package_data={'pylinac': ['watcher_config.yml']},
    zip_safe=False,  # allows users to view files in egg/distribution
    url='https://github.com/jrkerns/pylinac',
    keywords="""medical physics AAPM TG142 quality assurance starshot cbct vmat dynalog starshot linac Varian Elekta
             trajectory log kv MV planar Leeds Las Vegas Standard Imaging PipsPro TG51""",
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    install_requires=["numpy >= 1.11",
                      "scipy >= 0.17",
                      "pydicom >= 0.9.9",
                      "matplotlib >= 1.4",
                      "scikit-image >= 0.12",
                      "scikit-learn >= 0.18",
                      "Pillow >= 4.0",
                      "tqdm == 3.8",
                      "click",
                      "pyyaml >= 3.10",
                      "yagmail",
                      "mpld3",
                      "reportlab >= 3.3",
                      "trf2csv"],
    entry_points={
        'console_scripts':
            ['pylinac=pylinac.scripts:cli']
    },
    license='MIT',
    test_suite='tests._test_all',
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
