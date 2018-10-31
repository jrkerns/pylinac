
=========
Changelog
=========

V 2.2.0
-------

General
^^^^^^^

* `#131 <https://github.com/jrkerns/pylinac/issues/131>`_ Typing has been added to almost every function and class in pylinac.
* F-strings have been incorporated. This bumps the minimum version for Python to 3.6.
* The ``publish_pdf`` method of every module has had its signature changed. Before, not all the signatures matched
  and only included a few parameters like author and unit name. This has been changed to
  ``filename: str, notes: str, list of str, open_file: bool, metadata: dict``. Filename and open file are straightforward.
  notes is a string or list of strings that are placed at the bottom of the report (e.g. 'April monthly redo'). Metadata is a dictionary that will print
  both the key and value at the top of each page of the report (e.g. physicist and date of measurement)
* The TG-51 module has been placed under a new module: :ref:`calibration_module`. This is because:
* A TRS-398 calibration module has been created :ref:`trs398`.
* The default colormap for arrays is now Viridis, the matplotlib default.
* A contributer's guide has been added: :ref:`contributer_guide`.
* `#141 <https://github.com/jrkerns/pylinac/issues/141>`_ The Pylinac logo has been included in the package so that PDFs can be generated without needing www access.
* A new dependency has been added: `argue <https://pypi.org/project/argue/>`_ which handles input parameters.


Flatness & Symmetry
^^^^^^^^^^^^^^^^^^^

* `#130 <https://github.com/jrkerns/pylinac/issues/130>`_ The flatsym module has been completely rewritten.
  Documentation has also been updated and should be consulted given the number of changes: :ref:`flatsym_module`.

VMAT
^^^^

* The overall simplicity of use has been increased by automating & removing several parameters.
* `#128 <https://github.com/jrkerns/pylinac/issues/128>`_ The ``VMAT`` class has been split into two classes: :class:`~pylinac.vmat.DRGS` and :class:`~pylinac.vmat.DRMLC`. Although there are now two classes
  instead of one, the overall simplicity has been increased, such as the following:

  * The ``test`` parameter in ``analyze()`` is no longer required and has been removed.
  * The ``type`` is no longer required in ``.from_demo_images()``.
  * The demo method matches the other modules: ``.run_demo()``
  * All naming conventions have been deprecated.
* The ``x_offset`` parameter has been removed. The x-position is now based on the FWHM of the DMLC field itself.
  This means the x-position is dynamic and automatic.
* The ``delivery_types`` parameter has been removed. The delivery types of the images are now automatically determined.
* The methods for plotting and saving subimages (each image & the profiles) has been converted to a private method
  (``_plot_subimage()``, ...). There is little need for a public method to plot individually.

TG-51/Calibration
^^^^^^^^^^^^^^^^^

* `#129 <https://github.com/jrkerns/pylinac/issues/129>`_ The TG-51 module has been refactored to add a ``TG51ElectronLegacy`` and ``TG51ElectronModern`` calibration class.
  The Legacy class uses the classic TG-51 values that require a kecal value and a Pgradient measurement. The Modern
  class uses the equations from Muir & Rogers 2014 to calculate kQ that updates and incorporates the Pgradient and
  kecal values. While not strictly TG-51, these values are very likely to be incorporated into the next TG-51 addendum
  as the kQ values for photons already have.
* Certain parameters have been refactored: ``volt_high`` and ``volt_low`` have been refactored to ``voltage_reference``
  and ``voltage_reduced``, ``m_raw``, ``m_low``, and ``m_opp`` have been refactored to ``m_reference``, ``m_reduced``,
  and ``m_opposite``. These parameters are also the same for the TRS-398 classes
* The ``kq`` function has been separated into three functions: ``kq_photon_pdd10x``, ``kq_photon_tpr2010``, and
  ``kq_electron``.
* A PDD(20,10) to TPR(20,10) converter function has been added.
  This can be used in either TG-51 or TRS-398 to get TPR without actually needing to measure it.
* Defaults were removed from most functions to avoid possible miscalibration/miscalculation.
* Most parameters were changed to be keyword only. This will prevent accidental miscalculations from simple positional arguments.
* `#127 <https://github.com/jrkerns/pylinac/issues/127>`_ A TRS-398 module has been added. There are two main classes: ``TRS398Photon`` and ``TRS398Electron``.

Bug Fixes
^^^^^^^^^
* `#138 <https://github.com/jrkerns/pylinac/issues/138>`_/`#139 <https://github.com/jrkerns/pylinac/issues/139>`_: Too
  many arguments when plotting the leaf error subplot for picketfence.
* `#133 <https://github.com/jrkerns/pylinac/issues/133>`_: Trajectory log HDMLC status was reversed. This only affected
  fluence calculations using the ``equal_aspect`` argument.
* `#134 <https://github.com/jrkerns/pylinac/issues/134>`_: Trajectory log fluence array values were not in absolute MU.


V 2.1.0
-------

General
^^^^^^^

* After reflection, the package seems to have bloated in some respects.
  Certain behaviors are only helpful in very few circumstances and are hard to maintain w/ proper testing.
  They are described below or in their respective sections.
* The command line commands have been deprecated. All commands were simply shortcuts that are just as easy to place in
  a 1-2 line Python script. There was no good use case for it in the context of how typical physicists work.
* The interactive plotting using MPLD3 has been deprecated. Matplotlib figures and PDF reports should be sufficient.
  This was a testing nightmare and no use cases have been presented.
* The transition of the method ``return_results()`` to ``results()`` is complete. This was baked-in from the very
  beginning of the package. It is expected that results would return something, nor is there any other corresponding
  method prefixed with ``return_``.
* Pip is now the recommended way to install pylinac. Packaging for conda was somewhat cumbersome. Pylinac itself is just
  Python and was always installable via pip; it is the dependencies that are complicated.
  The wheels format seems to be changing that.
* Some dependency minimum versions have been bumped.

CatPhan
^^^^^^^

* The module was refactored to easily alter existing and add new catphan models.
* The CatPhan HU module classifier has been deprecated. Its accuracy was not as high as the original brute force method.
  Thus, the ``use_classifier`` keyword argument is no longer valid.
* CatPhan 604 support was added thanks to contributions and datasets from `Alan Chamberlain <https://github.com/alanphys>`_.
  More datasets are needed to ensure robust analysis, so please contribute your dataset if it fails analysis.
* The CTP528 slice (High resolution line pairs) behavior was changed to extract the max value from 3 adjacent slices.
  This was done because sometimes the line pair slice selected was slightly offset from the optimum slice. Using the
  mean would lower MTF values. While using the max slightly increases the determined MTF from previous versions,
  the reproducibility was increased across datasets.

Winston-Lutz
^^^^^^^^^^^^

* Certain properties have been deprecated such as gantry/coll/couch vector to iso.
  These are dropped in favor of a cumulative vector.
* A BB shift vector and shift instructions have been added for iterative WL testing.
  I.e. you can get a BB shift to move the BB to the determined iso easily.

  .. code-block:: python

    import pylinac

    wl = pylinac.WinstonLutz.from_demo_images()
    print(wl.bb_shift_instructions())
    # output: RIGHT 0.29mm; DOWN 0.04mm; OUT 0.41mm
    # shift BB and run it again...

* Images taken at nonzero couch angles are now correctly accounted for in the BB shift.
* Images now do not take into account shifts along the axis of the beam (`#116 <https://github.com/jrkerns/pylinac/issues/116>`_).
* The name of the file will now not automatically be interpreted if it can. This could cause issues for valid DICOM files that had sufficient metadata.
  If the image was taken at Gantry of 45 and the file name contained "gantry001" due to, e.g., TrueBeam's default naming convention it would override the DICOM data.
  (`#124 <https://github.com/jrkerns/pylinac/issues/124>`_)

Picket Fence
^^^^^^^^^^^^

* Files can now allow for interpretation by the file name, similar to the WL module. This is helpful for Elekta linacs that may be doing this test (`#126 <https://github.com/jrkerns/pylinac/issues/126>`_).

Core Modules
^^^^^^^^^^^^

* ``is_dicom`` and ``is_dicom_image`` were moved from the ``utilites`` module to the ``io`` module.
* ``field_edges()`` had the parameter ``interpolation`` added so that field edges could be computed more accurately (`#123 <https://github.com/jrkerns/pylinac/issues/123>`_)
* A new class was created called ``LinacDicomImage``. This is a subclass of ``DicomImage`` and currently adds smart gantry/coll/couch angle interpretation but may be extended further in the future.


V 2.0.0
-------

General
^^^^^^^

* Version 2.0 is here! It may or may not be a real major version update worthy of '2.0', but '1.10' just didn't sound as good =)
* A GUI has been added! Most major modules have been added to the GUI. The GUI is a very simple
  interface that will load files and publish a PDF/process files. To start the gui run the `gui()` function like
  so:

  .. code-block:: python

    import pylinac
    pylinac.gui()

  You may also start the GUI from the command line:

  .. code-block:: bash

    pylinac gui

  The GUI is a result of a few causes. Many physicists don't know how to code; this should remove that barrier
  and allow Pylinac to get even more exposure. I have always felt the web was the future, and it likely is, but
  pylinac should be able to run on it's own, and because a rudimentary GUI is relatively easy, I've finally made it.
  The GUI is also free to use and has no hosting costs (unlike assuranceQA.com). Also, due to other ventures, a new job, and a
  newborn, I couldn't devote further time to the assuranceQA site--A native GUI is much easier
  albeit much more primitive.
* Some module PDF methods now don't require filenames. If one is not passed it will default to the name of the file analyzed.
  E.g. "abc123.dcm" would become "abc123.pdf". Modules where multiple images may be passed (e.g. a CBCT directory) still requires a filename.
* PDF methods now have a boolean parameter to open the file after publishing: ``open_file``.
* A number of dependencies have been bumped. Some were for specific reasons and others were just out of good practice.

Watcher
^^^^^^^

* Closes `#84 <https://github.com/jrkerns/pylinac/issues/84>`_ Which would overwrite the resulting zip and PDF of
  initially unzipped CBCTs performed on the same day. I.e. multiple CBCTs would result in only 1 zip/PDF. The image
  timestamp has been edited so that it will include the hour-minute-second of the CBCT to avoid conflict.
* Closes `#86 <https://github.com/jrkerns/pylinac/issues/86>`_ - Which had a discrepancy between the YAML config setting of the file source directories
  and what the watcher was looking for.

CatPhan
^^^^^^^

* Closes `#85 <https://github.com/jrkerns/pylinac/issues/85>`_ Which displayed the nominal CBCT slice width on PDF reports,
  not the detected width for the CatPhan504 & CatPhan600.
* Closes `#89 <https://github.com/jrkerns/pylinac/issues/89>`_ which had variables swapped in the CatPhan503 PDF.
* The ``contrast_threshold`` parameter has been renamed to ``cnr_threshold``. The meaning and values are the same, but has been
  renamed to be consistent with other changes to the ``roi`` module.
* Due to various problems with the SVM classifier, the default setting of the classifier has been set to ``False``.

Planar Phantoms
^^^^^^^^^^^^^^^

* The Las Vegas phantom has been added to the planar imaging module. It's use case is very similar to the existing planar
  phantoms:

  .. code-block:: python

    from pylinac import LasVegas

    lv = LasVegas('myfile.dcm')
    lv.analyze()
    lv.publish_pdf()
    ...

* The :meth:`pylinac.planar_imaging.LeedsTOR.analyze` method has an additional parameter: `angle_offset`. From analyzing multiple Leeds images, it has become
  apparent that the low contrast ROIs are not always perfectly set relative to the phantom. This parameter will allow the user
  to fine-tune the analysis to perfectly overlay the low contrast ROIs by adding an additional angle offset to the analysis.

Winston-Lutz
^^^^^^^^^^^^

* Closes enhancement `#63 <https://github.com/jrkerns/pylinac/issues/63>`_ Files can now have the axis settings interpreted via the file name.
  E.g: "myWL_gantry90_coll0_couch340.dcm". See :ref:`using_file_names_wl` for further info.
* The `x/y/z_offset` properties of the WLImages which were deprecated many versions ago have finally been removed.
* The `collimator/gantry_sag` and associated `plot_gantry_sag` methods have been deprecated. A similar method has been implemented that utilizes the RMS deviation.
  To achieve the "gantry sag" using RMS errors use the method `axis_rms_deviation` with parameter `value='range'`.

TG-51
^^^^^

* The Electron class has been adjusted to reflect the `Muir & Rodgers 2014`_ kecal data which allows the user to calculate kQ from just R50 data.
* The `kq` function now accepts an `r_50` parameter to calculate kQ based on the above data.

.. _Muir & Rodgers 2014: http://onlinelibrary.wiley.com/doi/10.1118/1.4893915/abstract

Core Modules
^^^^^^^^^^^^

* The `Image` class has been fully depricated and is no longer available. Use the functions available in the :module:`pylinac.core.image` module instead.
  See the version 1.4.0 release notes for further details.
* The `remove_edges` method has been deprecated and is now an alias for `crop`. The `crop` method should be used instead. Parameters are exactly the same.

V 1.9.0
-------

General Changes
^^^^^^^^^^^^^^^

* This release introduces PDF reports for most major modules. All classes with this functionality
  have been given a ``publish_pdf`` method. This method takes an output filename and other optional
  data like the author, machine/unit, and any custom notes. See e.g. :meth:`pylinac.starshot.Starshot.publish_pdf`
  or :meth:`pylinac.picketfence.PicketFence.publish_pdf`.
* The watch/process functions have been tweaked to best work on one unit per run. Multiple units/machines should
  have their own config files. A new article :ref:`task_scheduler` describes how to use the process function with Windows Task
  Scheduler to regularly pull and analyze files.

CatPhan
^^^^^^^

* The CatPhan classes, when passed a directory during instantiation, will search through the DICOM files
  for Series UIDs and analyze the files of the most numerous UID. E.g. if a folder has 80 DICOM images including
  one set of 60 CBCT images and a total of 20 VMAT and picket fence images, it will find the CBCT files via UID and analyze
  those, leaving the other images/files alone. This is useful for when all QA images are simply dumped into one folder.
* Raw, uncompressed CatPhan DICOM files can optionally be compressed to a ZIP file after analysis using the new ``zip_after``
  argument in the ``analyze`` method.

Watcher/Processer
^^^^^^^^^^^^^^^^^

* The ``watcher``/``process`` functions have been reworked to produce PDF files rather than PNG/txt files.
* If upgrading the watch/process function from a previous pylinac version be sure to copy/amend the new default YAML config file
  as new keywords have been added and using old YAML files will error out.
* Several new configuration keywords have been changed/added. In the general section, ``use-classifier``
  has been deprecated in favor of individual module keywords of the same name. This allows a user to use a
  classifier for, say, picket fence images but not for winston lutz images. A ``unit`` keyword has been added
  that specifies which unit the files should be considered to be from. This unit name is passed to the PDF
  reports that are generated. If you have multiple units, make individual YAML configuration files, one for each
  unit.
* CatPhan, VMAT, and Winston-Lutz can now take raw, unzipped images as well as the usual ZIP archive. ZIP archives
  are detected only by keywords as usual. For uncompressed CatPhan images, the analyzer will look for any CatPhan DICOM
  file groups via UID (see above CatPhan section), analyze them, and then ZIP the images until no further sets can be found.
  For VMAT and Winston-Lutz if the ``use-classifier`` setting is true their respective sections in the YAML configuration
  then an image classifier is used to group images of the given type and then analyze them.

v 1.8.0
-------

General Changes
^^^^^^^^^^^^^^^

* This release focuses solely on the CBCT/CatPhan module.
* Pylinac now has a logo! Check out the readme on github or landing page on ReadTheDocs.

Watcher/Processer
^^^^^^^^^^^^^^^^^

* The cbct analysis section has been renamed to ``catphan``. Thus, the YAML config file needs to look like the
  following::

    # other sections
    ...

    catphan:  # not cbct:
        ...

    ...


CBCT/CatPhan
^^^^^^^^^^^^

* The Python file/module has been renamed to ``ct`` from ``cbct``. E.g.::

    from pylinac.ct import ...

  Most users import directly from pylinac, so this should affect very few people. This was done to generalize
  the module to make way for other CT/CBCT phantoms that pylinac may support in the future.
* The CBCT module can now support analysis of the CatPhan 600.
* Automatic detection of the phantom is no longer be performed. Previously, it depended on the
  manufacturer to determine the phantom (Varian->504, Elekta->503), but that did not consider users scanning the
  CatPhan in their CT scanners, which would give inconsistent results.
* Due to the above, separate classes have been made for the CatPhan models. I.e. flow looks like this now::

    # old way
    from pylinac import CBCT
    ...

    # new way
    from pylinac import CatPhan504, CatPhan600
    cat504 = CatPhan504('my/folder')
    cat600 = CatPhan600.from_zip('my/zip.zip')

* A classifier has been generated for each CatPhan. Thus, if loading a 503, a 503 classifier will be used, rather
  than a general classifier for all phantoms.
* The ``use_classifier`` parameter has been moved from the ``analyze()`` method to the class instantiation
  methods like so::

    from pylinac import CatPhan504
    cat504 = CatPhan504('my/folder', use_classifier=True)
    cat504.analyze()  # no classifier argument

* MTF is now more consistently calculated. Previously, it would simply look at the first 6 line pair regions.
  In cases of low mA or very noisy images, finding the last few regions would error out or give inconsistent results.
  Contrarily, high dose/image quality scans would only give MTF down to ~50% since the resolution was so good.
  Now, MTF is searched for region-by-region until it cannot find the correct amount of peaks and valleys, meaning it
  is now lost in the noise. This means high-quality scans will find and calculate MTF over more regions and fewer for
  low-quality scans. In general, this makes the MTF plot much more consistent and usually always gives the RMTF down to
  0-20%.
* Individual modules are now only composed of 1 slice rather than averaging the nearby slices. Previously, for consistency,
  a given module (e.g. CTP404) would find the correct slice and then average the pixel values of the slices on either side
  of it to reduce noise and give more consistent results. The drawback of this method is that results that depend on the
  noise of the image are not accurate, and signal/noise calculations were always higher than reality if only looking at
  one slice.


v 1.7.2
-------

* Fixed `(#78) <https://github.com/jrkerns/pylinac/issues/78>`_ - Certain CBCT datasets have irregular background
  values. Additionally, the dead space in the square CT dataset outside the field of view can also be very different
  from the air background. This fix analyzes the dataset for the air background value and uses that as a baseline value
  to use as a CatPhan detection threshold.

V 1.7.0
-------

General Changes
^^^^^^^^^^^^^^^

* The underlying structure of the watcher script has been changed to use a different framework. This change allows
  for analysis of existing files within the directory of interest.
* A new module has been introduced: ``tg51``, handling several common equations and data processing for things
  relating to TG-51 absolute dose calibration such as Kq, PDDx, Dref, pion, ptp, etc. It also comes with classes for
  doing a full TG-51 calculation for photons and electrons with cylindrical chambers.

Log Analyzer
^^^^^^^^^^^^

* The log analyzer has changed from having a main class of ``MachineLog``, to the two distinct log types:
  ``Dynalog`` and ``TrajectoryLog``. These classes are used the same way as machinelog, but obviously is meant for
  one specific type of log. This allows for cleaner source code as the ``MachineLog`` class had large swaths of
  if/else clauses for the two log types. But don't worry! If you're unsure of the log type or need to handle both
  types then a helper function has been made: ``load_log``. This function will load a log just like the ``MachineLog``
  did and as the new classes. The difference is it will do automatic log type detection, returning either a Dynalog
  instance or TrajectoryLog instance. The ``MachineLogs`` class remains unchanged.
* More specific errors have been introduced; specifically ``NogALogError``, ``NotADynalogError``, and ``DynalogMatchError``
  which are self-explanatory and more specific than ``IOError``.
* Fixed `(#74) <https://github.com/jrkerns/pylinac/issues/74>`_ which was causing Dynalogs with patient names containing
  a "V" to be classified as Trajectory logs.
* Fixed `(#75) <https://github.com/jrkerns/pylinac/issues/75>`_ which was skewing gamma pass percent values.

Planar Imaging
^^^^^^^^^^^^^^

* The ``PipsProQC3`` class/phantom has been refactored to correctly reflect its manufacturer to Standard Imaging,
  thus the class has been renamed to ``StandardImagingQC3``.

Directory Watching
^^^^^^^^^^^^^^^^^^

* The ``watch`` command line argument now has a sister function, available in a regular Python program:
  :func:`~pylinac.watcher.watch`.
  With this command you can run the directory watcher programmatically, perfect for continuous log monitoring.
* A new command line argument is available: ``process``. This command is also available in Python as
  :func:`~pylinac.watcher.process`
  which can be called on a directory either through the command line or programmatically and will analyze a
  folder once and then exit, perfect for analyzing a new monthly dataset.
* The structure of querying for files has been changed significantly. Instead of triggering on file changes (e.g. adding a
  new file to the directory), the watcher now constantly queries for new files at a specified interval. This means that
  when started, the watcher will analyze existing files in the folder, not just new ones.
* Information given in the email has been modified for logs, which may potentially contain PHI. Instead of the
  entire log file name given, only the timestamp is given. Additionally, the logs are no longer attached to the email.


V 1.6.0
-------

General Changes
^^^^^^^^^^^^^^^

* Changed the default colormap of dicom/grayscale images to be "normal" gray vs the former inverted gray.
  Brought up in `(#70) <https://github.com/jrkerns/pylinac/issues/70>`_ .
* Added a colormap setting that can be changed. See :ref:`changing_colormaps`
* Added a utility function :func:`~pylinac.core.utilities.clear_data_files` to clear demo files and classifier files.
  This may become useful for classifier updates. I.e. the classifier for a given algorithm can be cleared and updated as need be, without the
  need for a new package release. More information on this will follow as the use of classifiers becomes normal.
* Added a dependency to the pylinac requirements: `scikit-learn <http://scikit-learn.org/stable/>`_. This library will allow for machine learning
  advancements to be used with pylinac. I am aware of the increasing number of dependencies; pylinac has reached
  a plateau I believe in terms of advancement and I hope that this is the last major dependency to be added.

Winston-Lutz
^^^^^^^^^^^^

* `(#69) <https://github.com/jrkerns/pylinac/issues/69>`_ Added EPID position tracking. Now the EPID location will show up in images and will
  give an output value when printing the summary. Relevant methods like :meth:`~pylinac.winston_lutz.WinstonLutz.cax2epid_distance` and
  :meth:`~pylinac.winston_lutz.WinstonLutz.epid_sag`, and :meth:`~pylinac.winston_lutz.WinstonLutz.plot_epid_sag` have been added.
  The summary plot has also been changed to include two sag plots: one for the gantry and one for the EPID.
* Certain properties of WL images have been deprecated. ``x_offset`` has been replaced by :func:`~pylinac.winston_lutz.WLImage.bb_x_offset` and respectively
  for the other axes. Usage of the old properties will raise a deprecation warning and will be removed in v1.7.

  .. note::

    The deprecation warnings may not show up, depending on your python version and/or warning settings. See
    the `python docs <https://docs.python.org/3.5/library/warnings.html#warning-categories>`_ for more info.

CBCT
^^^^

* Added a Support Vector Machine classifier option for finding the HU slice. The classifier is faster (~30%) than
  the brute force method. This option is available as a parameter in the :meth:`~pylinac.cbct.CBCT.analyze` method as ``use_classifier``.
  In the event the classifier does not find any relevant HU slices, it will gracefully fall back to the brute force
  method with a runtime warning. Because of the fallback feature, the classifier is now used first by default.
  Using the classifier requires a one-time download to the demo folder, which happens automatically; just make sure
  you're connected to the internet.

Picket Fence
^^^^^^^^^^^^

* An ``orientation`` keyword argument was added to the :meth:`~pylinac.picketfence.PicketFence.analyze` method. This defaults to ``None``,
  which does an automatic determination (current behavior). In the event that the determined orientation was wrong, this argument can be utilized.

Watcher Service
^^^^^^^^^^^^^^^

* A new option has been added to the ``general`` section: ``use-classifier``. This option tells pylinac whether
  to use an SVM image classifier to determine the type of image passed. This allows the user not to worry about the
  file names; the images can be moved to the monitored folder without regard to naming. The use of the classifier
  does not exclude file naming conventions. If the classifier does not give a good prediction, the algorithm will
  gracefully fall back to the file name convention.

  The following image types currently support automatic detection:

  - Picket Fence
  - Starshot
  - Leeds TOR
  - PipsPro QC-3

V 1.5.6
-------

* Adds the ``dtype`` keyword to ``DicomImage``'s init method.
* `(#66) <https://github.com/jrkerns/pylinac/issues/66>`_ - Fixed an issue with Winston-Lutz
  isocenters not calculating correctly.
* `(#68) <https://github.com/jrkerns/pylinac/issues/68>`_ - Fixed the order of the Winston-Lutz images when plotted.
* Many thanks to Michel for noting the WL errors and `submitting the first external pull request <https://github.com/jrkerns/pylinac/pull/67>`_ !
* Fixed several small bugs and runtime errors.

V 1.5.5
-------

* `(#65) <https://github.com/jrkerns/pylinac/issues/65>`_ - Fixed the FlatSym demo file usage.

V 1.5.4
-------

* `(#64) <https://github.com/jrkerns/pylinac/issues/64>`_ - Fixed the Picket Fence offset from CAX value, which previously were all the same value.

V 1.5.1-3
---------

General Changes
^^^^^^^^^^^^^^^

* Fixed conda entry points so that the user can use pylinac console scripts.
* Moved demo images outside the package to save space. Files are downloaded when relevant methods are invoked.

V 1.5.0
-------

General Changes
^^^^^^^^^^^^^^^

* The pylinac directory watcher service got a nice overhaul. Now, rather than running the watcher script file directly, you
  can use it via the console like so:

  .. code-block:: bash

        $ pylinac watch "path/to/dir"

  This is accomplished through the use of console scripts in the Python setup file.
  Once you upgrade to v1.5, this console command immediately becomes available. See the updated docs on `Directory Watching <http://pylinac.readthedocs.org/en/latest/watcher.html>`_.
  Previously, customizing behavior required changing the watcher script directly. Now, a YAML file can be generated that contains all the
  analysis configurations. Create and customize your own to change tolerances and even to trigger emails on analyses.
* You can now anonymize logs via console scripts:

  .. code-block:: bash

       $ pylinac anonymize "path/to/log/dir"

  This script is a simple wrapper for the log analyzer's `anonymize <http://pylinac.readthedocs.org/en/stable/log_analyzer.html#pylinac.log_analyzer.anonymize>`_ function.

* Pylinac is now on `anaconda.org <https://anaconda.org/jrkerns/pylinac>`_ -- i.e. you can install via ``conda`` and forget about dependency & installation issues.
  This is the recommended way to install pylinac now. To install, add the proper channel to the conda configuration settings.

  .. code-block:: bash

        $ conda config --add channels jrkerns

  Then, installation and upgrading is as simple as:

  .. code-block:: bash

        $ conda install pylinac

  The advantage of saving the channel is that upgrading or installing in other environments is always as easy as ``conda install pylinac``.
* Pylinac's core modules (``image``, ``io``, etc) are now available via the root package level.

  .. code-block:: python

        # old way
        from pylinac.core import image
        # new way
        from pylinac import image

Starshot
^^^^^^^^

* Relative analysis is no longer allowed. I.e. you can no longer pass images that do not have a DPI or SID. If the image does not
  have these values inherently (e.g. jpg), you must pass it explicitly to the Starshot constructor. No changes are required for EPID images
  since those tags are in the image file.
* Added a ``.from_zip()`` class method. This can contain a single image (to save space) or a set of images that will be combined.

Log Analyzer
^^^^^^^^^^^^

* The `anonymize <http://pylinac.readthedocs.org/en/stable/log_analyzer.html#pylinac.log_analyzer.anonymize>`_ function received
  an optimization that boosted anonymization speed by ~3x for Trajectory logs and ~2x for Dynalogs. This function is *very* fast.
* Trajectory log subbeam fluences are now available. This works the same way as for the entire log:

  .. code-block:: python

    log = MachineLog.from_demo_dynalog()
    # calculate & view total actual fluence
    log.fluence.actual.calc_map()
    log.fluence.actual.plot_map()
    # calculate & view the fluence from the first subbeam
    log.subbeams[0].fluence.actual.calc_map()
    log.subbeams[0].fluence.actual.plot_map()

* The gamma calculation has been refactored to use the `image.gamma() <http://pylinac.readthedocs.org/en/stable/core_modules.html#pylinac.core.image.BaseImage.gamma>`_ method.
  Because of this, all ``threshold`` parameters have been changed to fractions:

  .. code-block:: python

    log = MachineLog.from_demo_trajectorylog()
    # old way
    log.fluence.gamma.calc_map(threshold=10)  # <- this indicates 10% threshold
    # new way
    log.fluence.gamma.calc_map(threshold=0.1)  # <- this also indicates 10% threshold

  The gamma threshold parameter requires the value to be between 0 and 1, so any explicit thresholds will raise an error that should be addressed.
* The ``.pixel_map`` attribute of the actual, expected, and gamma fluence structures have been renamed to ``array`` since they are numpy arrays. This
  attribute is not normally directly accessed so few users should be affected.

Bug Fixes
^^^^^^^^^

* Fixed a bug that would not cause certain imaging machine logs (CBCT setup, kV setups) to be of the "Imaging" treatment type.


V 1.4.1
-------

* `(#56) <https://github.com/jrkerns/pylinac/issues/56>`_ - Fixes a starshot issue where if the SID wasn't 100 it was corrected for twice.
* `(#57) <https://github.com/jrkerns/pylinac/issues/57>`_ - CR images sometimes have an RTImageSID tag, but isn't numeric; this caused SID calculation errors.


V 1.4.0
-------

General Changes
^^^^^^^^^^^^^^^

* Nearly all instance-based loading methods (e.g. ``Starshot().load('myfile')``) have been deprecated.
  Essentially, you can no longer do empty constructor calls (``PicketFence()``).
  The only way to load data is through the existing class-based methods (e.g. ``Starshot('myfile')``, ``Starshot.from_url('http...')``, etc).
  The class-based methods have existed for several versions, and they are now the preferred and only way as there is
  no use case for an empty instance.
* Since v1.2 most URLs were downloaded and then the local (but temporary) files were loaded. This practice has now been
  standardized for all modules. I.e. any ``from_url()``-style call downloads a temporary file and loads that. Because the
  downloads are to a temporary directory, then are removed upon exit.
* Loading images using the ``Image`` class has been deprecated (but still works) in favor of the new functions in the same module with the same name.
  Where previously one would do::

        from pylinac.core.image import Image

        img = Image.load('my/file.dcm')

  One should now do::

       from pylinac.core.image import load

       img = load('my/file.dcm')

  Functionality is exactly the same, but supports a better abstraction (there is no reason for a class for just behaviors).
  The same change applies for the other loading methods of the Image class: ``load_url`` and ``load_multiples``. The ``Image``
  class is still available but will be removed in v1.5.

Picket Fence
^^^^^^^^^^^^

* ``PicketFence`` can now load a machine log along with the image to use the expected fluence to determine error. This
  means if an MLC bank is systematically shifted it is now detectable, unlike when the pickets are fitted to the MLC peaks.
  Usage is one extra parameter::

      pf = PicketFence('my/pf.dcm', log='my/pf_log.bin')

Winston-Lutz
^^^^^^^^^^^^

* A ``from_url()`` method has been added.
* Upon loading, all files are searched within the directory, not just the root level.
  This allows for nested files to be included.

CBCT
^^^^

* The ``from_zip_file()`` class constructor method has been renamed to ``from_zip()`` to be consistent with the rest
  of pylinac's similar constructors.

Log Analyzer
^^^^^^^^^^^^

* A new ``treatment_type`` has been added for CBCT and kV logs: ``Imaging``.
* A new function has been added to the module: ``anonymize()``. This function is similar to the ``.anonymize()`` method,
  but doesn't require you to load the logs manually. The function is also threaded so it's very fast for mass anonymization::

     from pylinac.log_analyzer import anonymize

     anonymize('my/log/folder')
     anonymize('mylog.bin')

Starshot
^^^^^^^^

* The starshot minimization algorithm has been changed from `differential evolution <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution>`_ to the
  more predictable `minimize <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize>`_.
  Previously, results would *often* be predictable, but would occasionally give really good or really bad results even though no input
  was changed. This was due to the algorithm; now that a stable algorithm is being used, results are reproducible.

VMAT
^^^^

* The VMAT loading scheme got a few changes. The `Naming Convention <http://pylinac.readthedocs.org/en/latest/vmat_docs.html#naming-convention>`_
  is still the same, but images are always loaded upon instantiation (see General Changes). Also, if the naming convention isn't used,
  image delivery types can be passed in during construction; e.g.::

      VMAT(images=(img1, img2), delivery_types=['open', 'dmlc']

* Loading from a URL has been renamed from ``from_urls()`` to ``from_url()`` and assumes it points to a ZIP archive with the images inside.

Bug Fixes
^^^^^^^^^

* `(#47) <https://github.com/jrkerns/pylinac/issues/47>`_ - Fixes the trajectory log number of beam holds calculation. Thanks, Anthony.
* `(#50) <https://github.com/jrkerns/pylinac/issues/50>`_ - Fixes RMS calculations for "imaging" trajectory logs. Previously,
  the RMS calculation would return ``nan``, but now returns 0.
* `(#51) <https://github.com/jrkerns/pylinac/issues/51>`_ - Results of the starshot wobble were sometimes extremely high or low.
  This has been fixed by using a more stable minimization function.
* `(#52) <https://github.com/jrkerns/pylinac/issues/52>`_ - The starshot wobble diameter was incorrect. A recent change
  of the point-to-line algorithm from 2D to 3D caused this issue and has been fixed.
* `(#53) <https://github.com/jrkerns/pylinac/issues/53>`_ - The Winston-Lutz BB-finding algorithm would sometimes pick up noise, mis-locating the BB.
  A size criteria has been added to avoid detecting specks of noise.
* `(#54) <https://github.com/jrkerns/pylinac/issues/54>`_ - Imaging Trajectory logs, besides having no RMS calculation, was producing warnings when calculating
  the fluence. Since there is no fluence for kV imaging logs, the fluence now simply returns an 0'd fluence array.
* `(#55) <https://github.com/jrkerns/pylinac/issues/55>`_ - Dead pixels outside the field were throwing off the thresholding algorithm and not detecting
  the field and/or BB.

V 1.3.1
-------

* `(#46) <https://github.com/jrkerns/pylinac/issues/46>`_ - Fixes CBCT analysis where there is a ring artifact outside the phantom.
  Incidentally, analysis is sped up by ~10%.

V 1.3.0
-------

General Changes
^^^^^^^^^^^^^^^

* A new dependency has been added: `scikit-image <http://scikit-image.org/>`_. Given that pylinac is largely an image
  processing library, this is actually overdue. Several extremely helpful functions exist that are made use
  of in both the new modules and will slowly be incorporated into the old modules as needed.
  The package is easily installed via pip (``pip install scikit-image``)
  or via conda (``conda install scikit-image``) if using the Anaconda distribution. Finally, if simply upgrading
  pylinac scikit-image will automatically install via pip. For the sake of installation speed I'd recommend conda.
* ROI sampling for CBCT and Leeds classes have been sped up ~10x, making analysis moderately to much faster.
* All user-interface dialog functions/methods have been deprecated. E.g. ``PicketFence.from_UI()`` is
  no longer a valid method. To retain similar functionality use Tk to open your own dialog box and
  then pass in the file name. Specifically, this applies to the VMAT, Starshot, PicketFence, MachineLog(s),
  FlatSym, and CBCT classes. The original goal of pylinac was to be used for a standalone desktop application.
  The assuranceqa.com web interface is the successor to that idea and does not need those UI methods.

Planar Imaging
^^^^^^^^^^^^^^

* A new planar imaging class has been added:
  `PipsProQC3 <http://pylinac.readthedocs.org/en/latest/planar_imaging.html#pipspro-phantom>`_.
  This class analyzes the PipsPro QC-3 MV imaging phantom. The class locates and analyzes low and high contrast ROIs.
* The Leeds phantom utilizes the scikit-image library to do a canny edge search to find the phantom.
  This will bring more stability for this class.

V 1.2.2
-------

* `(#45) <https://github.com/jrkerns/pylinac/issues/45>`_ Fixes various crashes of Leeds analysis.

V 1.2.1
-------

* `(#44) <https://github.com/jrkerns/pylinac/issues/44>`_ Fixed a stale wheel build causing ``pip install`` to install v1.1.

V 1.2.0
-------

General Changes
^^^^^^^^^^^^^^^

* CatPhan 503 (Elekta) analysis is now supported.
* A new planar imaging module has been added for 2D phantom analysis; currently the Leeds TOR phantom is available.
* The ``requests`` package is no longer needed for downloading URLs; the urllib stdlib module is now used instead.
* Requirements were fixed in the docs and setup.py; a numpy function was being used that was introduced in
  v1.9 even though v1.8 was stated as the minimum; the new requirement is v1.9.
* Demonstration methods for the main classes have been fully converted to static methods. This means, for example,
  the following are equivalent: ``CBCT().run_demo()`` and ``CBCT.run_demo()``.

Core Modules
^^^^^^^^^^^^
* A tutorial on the use of the core modules is now available.
* A new ``mask`` core module was created for binary array operations.
* `(#42) <https://github.com/jrkerns/pylinac/issues/42>`_ The Image classes now have a :class:`~pylinac.core.image.ImageMixin.gamma` method available.
* The Image classes' ``median_filter()`` method has been renamed to :meth:`~pylinac.core.image.ImageMixin.filter`, which allows for different types
  of filters to be passed in.
* The Image class can now load directly from a URL: :meth:`~pylinac.core.image.Image.load_url`.

CBCT
^^^^

* CatPhan 503 (Elekta) is now supported. Usage is exactly the same except for the low-contrast module, which
  is not present in the 503.
* The low contrast measurements now use two background bubbles on either side of each contrast ROI. The default contrast
  threshold has been bumped to 15, which is still arbitrary but fits most eyeball values.

Starshot
^^^^^^^^

* `(#43) <https://github.com/jrkerns/pylinac/issues/43>`_ Keyword arguments can be passed to the init and class methods regarding the image info. For example,
  if a .tif file is loaded but the DPI is not in the image header it can be passed in like so:

  .. code-block:: python

     star = Starshot('mystar.tif', dpi=100, sid=1000)

Planar Imaging
^^^^^^^^^^^^^^

* 2D analysis of the Leeds TOR phantom is available. Tests low and high contrast.
  A new :ref:`planar_imaging` doc page has been created.

Winston-Lutz
^^^^^^^^^^^^

* A :meth:`~pylinac.winston_lutz.WinstonLutz.save_summary` method has been added for saving the plot to file.

V 1.1.1
-------

* Winston-Lutz demo images were not included in the pypi package.

V 1.1.0
-------

General Changes
^^^^^^^^^^^^^^^

* This release debuts the new Winston-Lutz module, which easily loads any number of EPID images,
  finds the field CAX and the BB, and can plot various metrics.

Log Analyzer
^^^^^^^^^^^^

* Logs can now be anonymized using the ``.anonymize()`` method for both MachineLog and MachineLogs.
* The ``.to_csv()`` methods for MachineLog and MachineLogs returns a list of the newly created files.
* MachineLogs can now load from a zip archive using ``.from_zip()``.

V 1.0.3
-------

* Fixes #39. MachineLog fluence was inverted in the left-right direction.
* Fixes #40. MachineLog fluence calculations from dynalogs were dependent on the load order (A-file vs. B-file).

V 1.0.2
-------

* Fixes #38. MachineLog fluence calculations would crash if there was no beam-on snapshots (e.g. kV images).

V 1.0.1
-------

* Fixes #37. Reading in a trajectory log txt file with a blank line caused a crash.

V 1.0.0
-------

General Changes
^^^^^^^^^^^^^^^

* This release debuts the new interactive plotting for certain figures.
  Quickly, matplotlib line/bar plots (althouth not yet images/arrays) can be plotted and saved in HTML using the MPLD3 library.
  This is less of interest to users doing interactive work, but this adds the ability to embed HTML plots in web pages.
* Several numpy array indexing calls were converted to ints from floats to avoid the new 1.9 numpy type-casting warnings.
  This also speeds up indexing calls slightly.

Picket Fence
^^^^^^^^^^^^

* The analyzed image now has the option of showing a leaf error subplot beside the image. The image is aligned
  to the image such that the leaves align with the image.

Starshot
^^^^^^^^

* Plotting the analyzed starshot image now shows both the zoomed-out image and a second, zoomed-in view of the wobble.
* Each subplot can be plotted and saved individually.

VMAT
^^^^

* Plotting the analyzed image now shows the open and dmlc images and the segment outlines as well as a profile comparison
  between the two images. Each subplot can also be plotted and saved individually.
* ``MLCS`` is no longer a test option; ``DRMLC`` should be used instead.


V 0.9.1
-------

* Fixed a bug with the log analyzer treatment type property.


V 0.9.0
-------

General Changes
^^^^^^^^^^^^^^^

* This release has a few new features for the CBCT class, but is mostly an internal improvement.
  If you only use the main classes (CBCT, PicketFence, Starshot, etc), there should be no changes needed.

CBCT
^^^^

* The CBCT analysis now examines low contrast ROIs and slice thickness.
* CBCT components have been renamed. E.g. the HU linearity attr has been renamed ``hu`` from ``HU``.

Starshot
^^^^^^^^

* Fixes #32 which was causing FWHM peaks on starshots to sometimes be erroneous for uint8/uint16 images.

PicketFence
^^^^^^^^^^^

* Adds #31, a method for loading multiple images into PicketFence.

Log Analyzer
^^^^^^^^^^^^

* Fixes a bug which sometimes caused the parsing of the associated .txt log file for trajectory logs
  to crash.


V 0.8.2
-------

* Fixed a bug with the picket fence overlay for left-right picket patterns.
* Plots for starshot, vmat, and picketfence now have a larger DPI, which should mean some more
  detail for saved images.


V 0.8.1
-------

* Fixed an import bug


V 0.8.0
-------

General Changes
^^^^^^^^^^^^^^^

* An upgrade for the robustness of the package. A LOT of test images were added for the Starshot, CBCT, PicketFence, and VMAT modules and
  numerous bugs were caught and fixed in the process.
* The debut of the "directory watcher". Run this script to tell pylinac to watch a directory; if a file with certain keywords is placed in the directory,
  pylinac will analyze the image and output the analyzed image and text file of results in the same directory.
* A generic troubleshooting section has been added to the documentation, and several modules have specific troubleshooting sections to help identify common errors
  and how to fix them.

VMAT
^^^^

* Added a ``from_zip()`` and ``load_zip()`` method to load a set of images that are in a zip file.
* Added an ``x_offset`` parameter to ``analyze()`` to make shifting segments easier.

PicketFence
^^^^^^^^^^^

* Fixed #30, which wasn't catching errors on one side of the pickets, due to a signed error that should've been absolute.
* Two new parameters have been added to ``analyze()``: ``num_pickets`` and ``sag_adjustment``, which are somewhat self-explanatory.
  Consult the docs for more info.

Starshot
^^^^^^^^

* Fixed #29, which was causing analysis to fail for images with a pin prick.

CBCT
^^^^

* Fixed #28, which was applying the phantom roll adjustment the wrong direction.


V 0.7.1
-------

General Changes
^^^^^^^^^^^^^^^

* Added ``.from_url()`` class method and ``.load_url()`` methods to most modules.

PicketFence
^^^^^^^^^^^

* Fixed #23, which was not properly detecting pickets for picket patterns that covered less than half the image.
* Fixed #24, which was failing analysis from small but very large noise. A small median filter is now applied to images upon loading.


V 0.7.0
-------

General Changes
^^^^^^^^^^^^^^^

* The scipy dependency has been bumped to v0.15 to accommodate the new differential evolution function using in the Starshot module.

CBCT
^^^^

* Whereas v0.6 attempted to fix an issue where if the phantom was not centered in the scan it would error out by adding
  a z-offset, v0.7 is a move away from this idea. If the offset given was not correct then analysis would error disgracefully.
  It is the point of automation to automatically detect things like where the phantom is in the dataset. Thus, v0.7 is a move
  towards this goal. Briefly, upon loading all the images are scanned and the HU linearity slice is searched for. Of the detected
  slices, the median value is taken. Other slices are known relative to this position.
* As per above, the z-offset idea is no longer used or allowed.
* Plots are now all shown in grayscale.
* If the phantom was not completely scanned (at least the 4 modules of analysis) analysis will now error out more gracefully.


V 0.6.0
-------

General Changes
^^^^^^^^^^^^^^^

* Pylinac now has a wheel variation. Installation should thus be quicker for users with Python 3.4.
* Most main module classes now have a save method to save the image that is plotted by the plot method.

Class-based Constructors
########################

* This release presents a normalized and new way of loading and initializing classes for the PicketFence, Starshot, VMAT and CBCT classes.
  Those classes all now accept the image path (folder path for CBCT) in the initialization method. Loading other types of data
  should be delegated to class-based constructors (e.g. to load a zip file into the CBCT class, one would use
  ``cbct = CBCT.from_zip_file('zfiles.zip')``). This allows the user to both initialize and load the images/data
  in one step. Also prevents user from using methods before initialization (i.e. safer). See ReadTheDocs page for more info.

Dependencies
############

* Because the VMAT module was reworked and is now based on Varian specs, the pandas package will no longer be required. FutureWarnings have been removed.

CBCT
^^^^

* Bug #18 is fixed. This bug did not account for slice thickness when determining the slice positions of the
  relevant slices.
* Bug #19 is fixed. This bug allowed the loading of images that did not belong to the same study. An error is now raised
  if such behavior is observed.
* Demo files are now read from the zipfile, rather than being extracted and then potentially cleaning up afterward. Behavior
  is now quicker and cleaner.
* Individual plots of certain module/slices can now be done. Additionally, the MTF can be plotted.
* The user can now adjust the relative position of the slice locations in the event the phantom is not set up to calibration
  conditions.

Log Analyzer
^^^^^^^^^^^^

* Keys in the ``txt`` attr dict weren't stripped and could have trailing spaces. Keys are now stripped.

VMAT
^^^^

* Ability to offset the segments has been added.
    Complete overhaul to conform to new Varian RapidArc QA specs. This includes the following:
* Rather than individual samples, 4 or 7 segments are created, 5x100mm each.
* Deviation is now calculated for each segment, based on the average segment value.
* The ``DRMLC`` test has changed name to ``MLCS``. E.g. passing a test should be:
  ``myvmat.analyze('mlcs')``, not ``myvmat.analyze('drmlc')``; the latter will still work but raises a future warning.

Starshot
^^^^^^^^

* Fixed a bug where an image that did not have pixels/mm information would error out.
* Added a tolerance parameter to the analyze method.


V 0.5.1
-------

Log Analyzer
^^^^^^^^^^^^

* Axis limits are now tightened to the data when plotting log_analyzer.Axis data.
* Gamma map plot luminescence is now normalized to 1 and a colorbar was added.
* Bug #14 fixed, where Tlogs v3 were not loading couch information properly.
* Trajectory log .txt files now also load along with the .bin file if one is around.

Starshot
^^^^^^^^

* Multiple images can now be superimposed to form one image for analysis.

VMAT
^^^^

* ``load_demo_image()`` parameter changed from ``test_type`` to ``type``

V 0.5.0
-------

* A new flatness & symmetry module allows for film and EPID image analysis.
* The ``log_analyzer`` module now supports writing trajectory logs to CSV.
* A FutureWarning that pandas will be a dependency in later versions if it's not installed.

V 0.4.1
-------

* Batch processing of logs added via a new class.
* ~4x speedup of fluence calculations.

V 0.4.0
-------

* A Varian MLC picket fence analysis module was added;
  this will analyze EPID PF images of any size and either orientation.


V 0.3.0
-------

* Log Analyzer module added; this module reads Dynalogs and Trajectory logs from Varian linear accelerators.

Starshot
^^^^^^^^
* The profile circle now aligns with the lines found.
* Recursive option added to analyze for recursive searching of a reasonable wobble.

* Image now has a cleaner interface and properties

V 0.2.1
-------

* Demo files were not included when installed from pip

V 0.2.0
-------

* Python 2.7 support dropped.
  Python 3 has a number of features that Python 2 does not,
  and because this project is just getting started, I didn't want to support Python 2,
  and then eventually drop it as Python 3 becomes more and more mainstream.
* Internal overhaul.
  Modules are now in the root folder.
  A core module with specialized submodules was created with a number of various tools.
* Demo files were assimilated into one directory with respective subdirectories.
* VMAT module can now handle HDMLC images.
* CBCT module was restructured and is much more reliable now.
* method names normalized, specifically the `return_results` method, which had different names
  in different modules.
* Lots of tests added; coverage increased dramatically.

V 0.1.3
-------

Overall

A module for analyzing CBCT DICOM acquisitions of a CatPhan 504 (Varian) has been added.
The starshot demo files have been compressed to zip files to save space.
A value decorator was added for certain functions to enforce, e.g., ranges of values that are acceptable.
The "Files" directory was moved outside the source directory.
-Starshot now reports the diameter instead of radius

V 0.1.2
-------

A PyPI setup.py bug was not properly installing pylinac nor including demo files.
Both of these have been fixed.


V 0.1.1
-------

Several small bugs were fixed and small optimizations made.
A few methods were refactored for similarity between modules.


V 0.1.0
-------

This is the initial release of Pylinac. It includes two modules for doing TG-142-related tasks:
Starshot & VMAT QA

Versioning mostly follows standard semantic revisioning. However, each new module will result in a bump in minor release, while bug fixes
will bump patch number.