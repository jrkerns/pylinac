
==================
Directory Watching
==================

As of version 0.8, you can now set up a script to run that watches a directory for incoming files and analyzes them with
pylinac if certain keywords are in the file.

Setup
-----

.. note::
    To use the watcher script, `watchdog <http://pythonhosted.org/watchdog/>`_ must be installed.

To use the watcher script, locate the pylinac package on the computer. Then, open up a terminal and run the watcher.py file::

    > python "C:\path\to\pylinac\watcher.py"

This will watch the directory the watcher script is in. To watch a different directory, pass it as an argument::

    > python "C:\path\to\pylinac\watcher.py" "C:\path\to\watch"

A logger will notify when the script has started, when a file gets added, and what the analysis status is. If a file
gets analyzed successfully, a .png and .txt file with the same name as the original file will be generated in the directory.

Tweaking
--------

The watcher script runs using default values for keywords and tolerance. To adjust these values, change the ``keywords``
and ``args`` attributes. See the API docs.

The defaults are as follows:

.. note::
    Only .zip files are accepted for CBCT and VMAT analysis. Only .bin or .dlg files are
    acceptable for machine logs.

* **Starshots** - keywords: 'star'; radius: 0.8; tolerance: 1
* **Picket Fences** - keywords: 'pf', 'picket'; tolerance: 0.5mm; action tolerance: 0.3mm
* **CBCTs** - keywords: 'ct', 'cbct'; HU tolerance: 40HUs; scaling tolerance: 1mm
* **VMATs** - keywords: 'vmat', 'drgs', 'drmlc'; tolerance: 1.5%
* **Logs** - keywords: None; resolution: 0.1mm; distance to agreement: 1mm; dose to agreement: 1%; threshold: 10%

API Documentation
-----------------

.. autoclass:: pylinac.watcher.AnalyzeMixin
    :no-show-inheritance:

.. autoclass:: pylinac.watcher.AnalyzeStar
    :no-members:

.. autoclass:: pylinac.watcher.AnalyzePF
    :no-members:

.. autoclass:: pylinac.watcher.AnalyzeCBCT
    :no-members:

.. autoclass:: pylinac.watcher.AnalyzeVMAT
    :no-members:

.. autoclass:: pylinac.watcher.AnalyzeLog
    :no-members:
