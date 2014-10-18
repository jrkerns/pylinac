# **pylinac**
[![PyPi version](https://pypip.in/v/pylinac/badge.png)](https://crate.io/packages/pylinac/)
[![Supported Python versions](https://pypip.in/py_versions/pylinac/badge.svg)](https://pypi.python.org/pypi/pylinac/)
[![License](https://pypip.in/license/pylinac/badge.svg)](https://pypi.python.org/pypi/pylinac/)

Pylinac provides TG-142 quality assurance (QA) tools to Python programmers as well as non-programmers in the field of 
therapy medical physics. The package comes in two flavors: source-level and GUI-level. The source-level
allows programmers and those familiar with Python to create custom tests with pylinac while the GUI-level will implement
pylinac into a GUI and executable that any user can use, without having to program software.

Below are the tools currently available; tools will be added one at a time as they are developed.

+ **[Starshot Analysis][star]** -
    Tools for analyzing film or superimposed EPID images for gantry, collimator, or MLC star (or spoke) shots. Can determine
    the minimum circle that touches all the radiation spokes (wobble)
----------------

**To get started, visit the [Documentation][docs].**

Code is Python 2/3 compatible thanks to [future][fut]. 



[install]: https://github.com/Medical-Physics-Pythonistas/pylinac/wiki/Installation
[vmatqa]: https://github.com/Medical-Physics-Pythonistas/pylinac/wiki/Using-the-VMAT-module
[star]: https://github.com/Medical-Physics-Pythonistas/pylinac/wiki/Using-the-Starshot-module
[fut]: http://python-future.org/
