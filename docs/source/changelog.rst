
=========
Changelog
=========

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