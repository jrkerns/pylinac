.. _installation:

============
Installation
============

Installing pylinac is easy no matter your skill!
Determine where you're at and then read the relevant section:

I know Python already
---------------------

Great! To get started install via pip:

.. code-block:: bash

    $ pip install pylinac

.. note::

    If you plan on developing or contributing to pylinac, you'll want to set up a development environment. See :ref:`setup_dev_env`.

I'm new to Python
-----------------

That's okay! If you're not a programmer at all you'll have a few things to do to get up and running,
but never fear. Using pylinac requires not just the base language Python, but a few dependencies as well.
Since most physicists don't program, or if they do it's in MATLAB, this section will help jumpstart your use of not
just pylinac but Python in general and all its wonderful goodness! Getting started with Python takes some work to
get set up and running, but it's well worth the effort.

.. _distro_stack:

Get a Distribution Stack
^^^^^^^^^^^^^^^^^^^^^^^^

Scientific computing with Python requires some specialized packages which require some specialized computing libraries.
While it's possible you have those libraries (for some odd reason), it's not likely. Thus, it's often best to install
the libraries *pre-compiled*. There are several options out there; I'll list just a few. Be sure to download the 3.x version,
preferably the newest:

* `Anaconda <http://continuum.io/downloads#py34>`__ - Continuum Analytics provides this one-stop-shop for tons of
  scientific libraries in an easy to install format. Just download and run the installer. If you don't want to install
  all 200+ packages, a slimmer option exists: `Miniconda <http://conda.pydata.org/miniconda.html>`__, which only installs
  ``conda`` and python installation tools. You can then use `conda <http://conda.pydata.org/index.html>`__ to install packages individually.
  Here's the Anaconda `quick start guide <https://store.continuum.io/static/img/Anaconda-Quickstart.pdf>`__.

  .. note:: Unlike the other options, individual packages can be upgraded on demand using the ``conda`` tool.

* `WinPython <https://winpython.github.io/>`_ - (Windows only) This grassroots project functions similarly to Anaconda, where all
  packages are precompiled and run out of the box. There are no corporate sponsors for this project, so support is not
  guaranteed.

See `Scipy's Installation Options <http://www.scipy.org/install.html>`__ for more options.

.. warning:: Python(x,y) is not yet available for Python 3, so don't choose this to try running pylinac.

.. note::
   If this is the first/only Python distribution you'll be using it'd be a good idea to activate it when the
   installer prompts you.

.. note:: You can install multiple Python stacks/versions, but only one is "active" at any given time.


Get an IDE (optional)
^^^^^^^^^^^^^^^^^^^^^

If you come from MATLAB, it's helpful to realize that MATLAB is both a language and an Integrated Development Environment (IDE).
Most languages don't have an official IDE, and some people may tell you IDEs are a crutch. If being a cyborg with superpowers is a crutch, then
call me a cripple because I find them extremely useful. As with all power, it must be wielded carefully though. The option of getting an IDE
is completely up to you. If you want one, here are some options:

* `PyCharm <https://www.jetbrains.com/pycharm/>`__ - A fully-featured, rich IDE. It's arguably king of the heavyweights and *free*. At least try it.
  Here's the PyCharm `quick start guide <https://www.jetbrains.com/pycharm/quickstart/>`__.

  .. image:: https://www.jetbrains.com/lp/pycharm-pro/static/5-refactoring-34f412bbb8d494a1fa4c8239eb7f0d5b.png
     :height: 400px
     :width: 600px

* `Spyder <https://www.spyder-ide.org/>`__ - A MATLAB-like IDE with similar layout, preferred by many working in the scientific realm.
  Here are the `Spyder docs <https://docs.spyder-ide.org/current/index.html>`__.

  .. note:: Spyder is part of the Anaconda distribution.

  .. image:: http://1.bp.blogspot.com/-KfAKKK_YN38/TkaV08KWgLI/AAAAAAAAB-s/TEDUviTJBeU/s1600/spyder_ipython012b.png
     :height: 400px
     :width: 600px
