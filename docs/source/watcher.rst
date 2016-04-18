
.. _dir_watching:

==================
Directory Watching
==================

The directory watcher service allows one to "set it and forget it" and let pylinac analyze files that
are moved to an appointed directory. Results can be emailed upon analysis.
The service allows for configurable analysis settings and email settings.

Getting Started
---------------

.. note::
    To use the watcher service a few dependencies must be installed that are not required for normal pylinac usage;
    see :ref:`dependencies` for more.

To use the watcher service, make sure you have v1.5+, then open up a terminal and run the console command:

.. code-block:: bash

    $ pylinac watch "dir/to/watch"

The ``watch`` argument initiates a thread that runs in the terminal. The directory to start watching is also
required. A logger will notify when the script has started, when a file gets added, and what the analysis status is. If a file
gets analyzed successfully, a .png and .txt file with the same name as the originals plus a suffix (default is ``_analysis``) will be generated in the directory.
You can also set up an email service when analysis runs, described below.

How it works
------------

The watcher service keeps an eye on a directory and watches for any files being added. If a file is moved
into the directory, the file is immediately checked to see if pylinac can analyze it. Because many files use the
same format (e.g. DICOM), keywords are used to filter which type of analysis should be done. When a file is
deemed analysis-worthy, pylinac will then run the analysis automatically and generate a .png and .txt file with
the analysis summary image and quantitative results, respectively. If the email service is setup, an email
can be sent either on any analysis done or only on failing analyses.

Configuration
-------------

The watcher service runs using default values for keywords and tolerance. These values are in a YAML
configuration file. Pylinac comes with a default file and settings. You can make your own YAML config file
and pass that into the service initialization call:

.. code-block:: bash

    $ pylinac watch "dir/to/watch" --config="my/config.yaml"

The YAML configuration file is the way to change keywords, change analysis settings, and set up email service.
You can use/copy the `pylinac default YAML <https://github.com/jrkerns/pylinac/blob/master/pylinac/watcher_config.yaml>`_
file as a starting template.

.. note::
    Only .zip files are accepted for CBCT, VMAT, and Winston-Lutz analyses. The ZIP archive
    filename must contain a keyword.

Setting up Email
----------------

The pylinac watcher service allows the user to set up an email trigger. The user must supply a
Gmail account (...@gmail.com). The gmail account name and password must be supplied in the YAML
configuration file.

.. warning::
    It is strongly recommended to create an *ad hoc* account for the watcher service. To
    use the email service requires that the account have lower-than-normal security by nature of
    the non-Gmail origin (i.e. you didn't log in and send it yourself).

To allow gmail to send the emails, log into the gmail account and go to account settings. Go to the
sign in & security section. At the very bottom in the section "Connected apps & sites" will be an
option to "Allow less secure apps". Turn this ON. This account will now allow the watcher service to
send emails.

.. warning::
    I'll say it again: don't use your personal account for this. Create a new account for the
    sole purpose of sending pylinac analysis emails.

In the YAML config file you can set emails to be sent after every analysis or only when an analysis
fails. The emails will contain a simple message, let you know when it was analyzed, and where to find
the results. All the emails also have the results attached, so no need to dig for the files.

Default YAML Configuration
--------------------------

The default configuration is reproduced here. All options are listed. You may remove or add keywords at will.
The analysis options must match the parameter names exaclty (e.g. ``hu_tolerance``, ``doseTA``).

.. literalinclude:: ../../pylinac/watcher_config.yml
