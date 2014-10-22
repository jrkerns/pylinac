Pylinac
=======

.. image:: https://pypip.in/version/pylinac/badge.svg?text=version
    :target: https://pypi.python.org/pypi/pylinac/
    :alt: Latest Version

.. image:: https://pypip.in/py_versions/pylinac/badge.svg
    :target: https://pypi.python.org/pypi/pylinac/
    :alt: Supported Python versions

.. image:: https://pypip.in/license/pylinac/badge.svg
    :target: https://pypi.python.org/pypi/pylinac/
    :alt: License

Pylinac provides TG-142 quality assurance (QA) tools to Python programmers as well as non-programmers in the field of 
therapy medical physics. The package comes in two flavors: source-level and GUI-level. The source-level
allows programmers and those familiar with Python to create custom tests with pylinac while the GUI-level will implement
pylinac into a GUI and executable that any user can use, without having to program software.

Below are the tools currently available; tools will be added one at a time as they are developed.

* `Starshot Analysis <http://pylinac.readthedocs.org/en/latest/starshot_docs.html>`_ -
    Tools for analyzing film or superimposed EPID images for gantry, collimator, or MLC star (aka spoke) shots. Can determine
    the minimum circle that touches all the radiation spokes (wobble).
* `VMAT QA <http://pylinac.readthedocs.org/en/latest/vmat_docs.html>`_ -
    A module for analyzing EPID images after performing the Jorgensen et al tests for VMAT QA, specifically the Dose Rate & Gantry Speed 
    (DRGS) and Dose Rate MLC (DRMLC) tests. Can load the open and MLC field images and calculate segment ratios as per the Jorgensen implementation.


Documentation
-------------
To get started, install, run the demos, view the API docs, and learn the module design, visit the `Full Documentation <http://pylinac.readthedocs.org/en/latest/index.html>`_.

Contributing
------------

Contributions to pylinac can be many. The most useful thing a non-programmer can contribute is images to analyze. If
you have VMAT images, starshot images, machine log files, CBCT DICOM files, or anything else you want analyzed, email me
 or share then via dropbox, drive, etc. Obviously, don't send patient files, just QA files.


Code is Python 2/3 compatible thanks to `future <http://python-future.org/>`_.


