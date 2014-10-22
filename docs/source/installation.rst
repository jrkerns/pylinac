.. _installation:

============
Installation
============

Installing pylinac is easy no matter your skill! However, no matter your skill level,
please take a look at the Note About Dependencies section to avoid some initial frustration.
Determine where you're at and then read the relevant section:

I know Python already
---------------------

Great! To get started, simply install pylinac through ``pip install pylinac`` via PyPI.
Alternatively, you can install from the github repository directly::

    pip install https://github.com/jrkerns/pylinac.git

Finally, you can download the repo (from github) and run ``setup.py install`` in the root directory

I'm new to Python
-----------------

That's okay! If you're not a programmer at all you'll have a few things to do to get up and running,
but never fear. Using pylinac requires not just the base language Python, but a few dependencies as well (see below). If you're new to
programming or come from, say MATLAB, you probably want to use an Integrated Development Environment (IDE). Unlike MATLAB,
you have a choice in IDEs when using Python; in fact the choices are numerous! Two that are worth considering are PyCharm and Spyder. An
easy way to get all of this in one lump sum is to download a distribution "stack", which includes Python and a number of other packages
(see note below for options).

Dependencies
------------

Pylinac, as a scientific package, has fairly standard scientific dependencies (>= means at least that version or newer):

* numpy >= 1.8
* scipy >= 0.13
* matplotlib >= 1.3.1
* pydicom >= 0.9.8
* Pillow >= 2.5
* future >= 0.13

A Note About Dependencies
-------------------------

So, if you're new to Python or you don't already have the packages listed, particularly numpy/scipy,
it is best to get these from a source that has already compiled these packages. Installing them via ``pip`` won't work unless
you have the underlying libraries installed, which is often not the case. There are a number of places you can get them
pre-compiled. Check out `Scipy's installation options <http://www.scipy.org/install.html>`_ for more on scientific distribution stacks. For
people new to Python, the easiest solution is likely `Python(x,y) <https://code.google.com/p/pythonxy/>`_ or
`Anaconda <http://continuum.io/downloads>`_. Both include the `Spyder IDE <https://bitbucket.org/spyder-ide/spyderlib/overview>`_,
which appears similar to MATLAB's. Installing one of those will have you up and running in no time.


