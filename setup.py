from setuptools import setup, find_packages

# anti-patterns be damned; Python should come up with a better system
with open('requirements.txt') as f:
    required = f.read().splitlines()

from pylinac import __version__


setup(
    name='pylinac',
    version=__version__,
    packages=find_packages(),
    package_data={'pylinac': ['files/*.png',]},
    url='https://github.com/jrkerns/pylinac',
    keywords="""medical physics AAPM TG142 quality assurance starshot cbct vmat dynalog starshot linac Varian Elekta
             trajectory log kv MV planar Leeds Las Vegas Standard Imaging PipsPro TG51""",
    author='James Kerns',
    author_email='jkerns100@gmail.com',
    description='A toolkit for performing TG-142 QA-related tasks on a linear accelerator',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    install_requires=required,
    python_requires='>=3.6',
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
