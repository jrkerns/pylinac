.. _developer-guide:

=================
Developer's Guide
=================

.. note::

    In this context, "developer" refers to anyone who wants to set up ``pylinac`` locally,
    contribute to the code base, create a fork for personal use, or otherwise modify the code. A "contributor" is a looser definition and could refer to
    someone who reports bugs, requests features, etc. A contributor is not necessarily a developer.

This section will describe how to set up ``pylinac`` locally for development, either
for contribution to ``pylinac`` itself or for use as a customized version in your clinic.

.. important::

  This assumes a basic understanding of Python, git, the command line, and virtual environments.

  If you're new to git, see the `git guide <https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup>`__.

Prerequisites
-------------

You should have a basic understanding of and install the following:

* `Git <https://git-scm.com/>`__.
* An IDE or text editor of your choice (e.g., VSCode, PyCharm, etc.)
* `uv <https://docs.astral.sh/uv/>`__ - a Python tool for managing virtual environments and dependencies.

.. important::

    Install ``uv`` globally via the `Standalone installer <https://docs.astral.sh/uv/getting-started/installation/#standalone-installer>`__.
    Do not use the ``uv`` pypi package.
    This will let you use ``uv`` outside a single virtual environment, install other versions of Python, etc.

* ``nox`` - `A scripting and test automation tool <https://nox.thea.codes/en/stable/>`__. You do not need to explicitly install this
  as it is a developer dependency installed automatically when setting up your virtual environment, but it is helpful to know what it is as it is used often
  for various tasks in ``pylinac`` and is an all-around extremely useful tool for Python development.


Setting Up Your Development Environment
---------------------------------------

Fork
^^^^

Fork the pylinac repository and clone it to your local machine. You may use the command line or a GUI tool like GitHub Desktop, GitKraken, or your IDE's built-in git support.
I would recommend using a remote VCS like GitHub, Gitlab, etc which allows you to easily fork repositories
as well as keep your fork safe (vs local-only git repositories).

Python
^^^^^^

Python is an obvious prerequisite, but note that ``uv`` can install `Python itself <https://docs.astral.sh/uv/guides/install-python/>`__ and manage multiple versions of it!
This is extremely helpful to, say, test ``pylinac`` against multiple versions of Python.

.. code-block:: bash

    # install Python 3.11
    uv python install 3.11

    # or install older versions of Python
    uv python install 3.9 3.10

Virtual Environment
^^^^^^^^^^^^^^^^^^^

Set up a virtual environment and install the required dependencies using ``uv``:

.. code-block:: bash

  # in your cloned directory of pylinac
  uv venv
  uv pip install -e .[developer]

Documentation
^^^^^^^^^^^^^

To render documentation locally run:

.. tip::

    You don't need to install the documentation dependencies as this is handled for you by ``nox`` which manages
    its own isolated virtual environments.

.. code-block:: bash

    uv run nox -s serve_docs

This will build the documentation, start a local web server, and open the documentation in your default web browser.
Any changes you make to the documentation files will be automatically reflected in the browser momentarily.


Running Tests
^^^^^^^^^^^^^

.. warning::

   This section is still under construction. There is a PR open to allow for running public tests, but it is not yet merged.

To run tests you can directly use ``pytest``, ``nox``, or ``uv`` (they all use pytest under the hood, but can provide better isolation) to run the tests.


.. tab-set::

   .. tab-item:: pytest

      .. code-block:: bash

          # run all tests
          pytest

          # run a specific test file
          pytest tests_basic/test_starshot.py

   .. tab-item:: nox

      .. code-block:: bash

          # run all tests
          nox -s run_tests

   .. tab-item:: uv

    .. code-block:: bash

        # run all tests
        uv run pytest

Building Wheels
^^^^^^^^^^^^^^^

.. note::

    This section only applies to internal developers or developers who want to publish a fork to PyPI

To build a wheel:

.. code-block:: bash

    uv run nox -s build_wheel

To upload your wheel to PyPI, you can use the following command:

.. code-block:: bash

    uv run nox -s upload_wheel
