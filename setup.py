from setuptools import setup, find_packages

__version__ = '1.5.6'


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    package_data={'pylinac':
                      ['watcher_config.yml',]}, # http://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
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
                      "scikit-image >= 0.11",
                      "Pillow >= 2.5",
                      "tqdm == 3.8",],
    extras_require={
        'console': ["click", "pyyaml", "watchdog", "yagmail"],
        'interactive': ["mpld3"]
    },
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
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Image Recognition",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries"]
)
