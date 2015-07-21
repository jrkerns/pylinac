
=========================
VMAT module documentation
=========================

Overview
--------

.. automodule:: pylinac.vmat
    :no-members:

Running the Demo
----------------

To run the VMAT demo, create a script or start an interpreter session and input::

    from pylinac.vmat import VMAT
    VMAT().run_demo_drgs()
    # alternatively, you can run the MLC Speed demo by:
    VMAT().run_demo_mlcs()

Results will be printed to the console and a figure showing both the Open field and MLC field image will pop up::

    Dose Rate & MLC Speed
    Test Results (Tol. +/-1.5%): PASS
    Max Deviation: 0.437%
    Absolute Mean Deviation: 0.382%

.. image:: images/vmat_analyzed.png

Image Acquisition
-----------------

If you want to perform these specific QA tests, you'll need DICOM plan files that control the linac precisely to deliver the test fields.
These can be downloaded from my.varian.com. Once logged in, search for RapidArc and you should see two items called "RapidArc QA Test Procedures
and Files for TrueBeam"; there will be a corresponding one for C-Series. Use the RT Plan files and follow the instructions, not including the assessment procedure,
which is the point of this module. Save & move the images to a place you can use pylinac.

Typical Use
-----------

The VMAT QA analysis follows what is specified in the Varian RapidArc QA tests and assumes your tests will run the exact same way.
Import the class::

    from pylinac.vmat import VMAT

The minimum needed to get going is to:

* **Load images** -- Loading the EPID DICOM images into your VMAT class object can be done by passing the file path(s) or by using a UI to
  find and get the file(s). While each image can be loaded directly, using an easy naming convention can simplify the process. Just make
  sure the string 'open' (case-insensitive) is in the Open field filename::

      open_img = "C:/QA Folder/VMAT/open_field.dcm"  # note the 'open'
      dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"  # no 'open'
      myvmat = VMAT((open_img, dmlc_img))  # the order doesn't matter

  Even easier, you can load the files using a UI dialog box::

    myvmat = VMAT.from_images_UI()  # a dialog box will pop up for you to choose both images

  You can also load the images one at a time. This is helpful when the file names do not follow the naming convention::

      open_img = "path/to/img1.dcm"
      dmlc_img = "path/to/img2.dcm"
      myvmat = VMAT()
      myvmat.load_image(open_img, im_type='open')
      myvmat.load_image(dmlc_img, im_type='mlc')

  Or load the files individually from a UI dialog box::

      myvmat = VMAT()
      myvmat.load_image_UI(im_type='open')  # tkinter UI box will pop up
      myvmat.load_image_UI(im_type='mlc')

  .. note::
    In previous versions of pylinac, loading images was instance-method based and only allowed one image at a time,
    meaning loading looked like the 3rd example above, no matter the name. This behavior has been simplified in favor
    of initialization normalization and adding class-method constructors (``VMAT.from_X``). The reason for this is that
    certain actions should only be allowed until after the images are loaded. Furthermore, loading the images should always be
    the first action of the analysis sequence. By using class constructors, certain pitfalls and errors can be avoided.
    Don't worry though, the old behavior still works.

* **Analyze the images** -- This is where pylinac does its work. Once the images are loaded, tell VMAT to analyze the images. See the
  Algorithm section for details on how this is done. The test to run (whether DRGS or MLCS) needs to be specified. Tolerance
  can also be passed, but has a default value of 1.5%::

      # analyze
      myvmat.analyze(test='mlcs', tolerance=1.5)

* **View the results** -- The VMAT module can print out the summary of results to the console as well as draw a matplotlib image to show where the
  segments were placed and their values::

      # print results to the console
      print(myvmat.return_results())
      # view analyzed images
      myvmat.plot_analyzed_image()

Tips & Tricks
-------------

The VMAT test can be customized using the ``tolerance`` parameter to whatever tolerance is needed, whether simply
using a default (1.5), or a clinical value.

Furthermore, the location of the segments can be adjusted. Older QA tests had the segments slightly offset from center.
The new revision of tests is now centered. To account for this, a :class:`~pylinac.vmat.Settings` class has been made
and can be used to customize the position of segments. Both an ``x_offset`` and ``y_offset`` attribute can be adjusted::

    myvmat = VMAT.from_demo_images('drgs')
    myvmat.settings.x_offset = 20
    myvmat.analyze('drgs')

The test results are also not constricted to just printing out. Important results can be accessed directly.
Continuing from the example above::

    myvmat.avg_abs_r_deviation
    0.357
    myvmat.max_r_deviation  # regardless of sign
    -1.423
    myvmat.passed
    True

Algorithm
---------

The VMAT analysis algorithm is based on the Varian RapidArc QA Test Procedures for C-Series and Truebeam. Two tests
(besides Picket Fence, which has its own module) are specified. Each test takes 10x0.5cm samples, each corresponding to
a distinct section of radiation. A corrected reading of each segment is made, defined as:
:math:`R_{corr}(x) = \frac{R_{DRGS}(x)}{R_{open}(x)} * 100`. The reading deviation of each segment is calculated as:
:math:`R_{deviation}(x) = \frac{R_{corr}(x)}{R_{corr}} * 100 - 100`, where :math:`R_{corr}` is the average of all segments.

The algorithm works like such:

**Allowances**

* The images can be acquired at any SID.
* The images can be acquired with any EPID (aS500, aS1000, aS1200).

**Restrictions**

    .. warning:: Analysis can fail or give unreliable results if any Restriction is violated.

* The tests must be delivered using the DICOM RT plan files provided by Varian which follow the test layout of Ling et al.
* The images must be acquired with the EPID.

**Pre-Analysis**

* **Determine image scaling** -- Segment determination is based on offsets from the center pixel of the image. However,
  some physicists use 150 cm SID and others use 100 cm, and others can use a clinical setting that may be different
  than either of those. To account for this, the SID is determined and then scaling factors are determined to be
  able to perform properly-sized segment analysis.

**Analysis**

.. note::
    Calculations tend to be lazy, computed only on demand. This represents a nominal analysis
    where all calculations are performed.

* **Calculate sample boundaries** -- Segment positions are always the same and are determined by offsets from the center pixel given by Varian.
  These values are then scaled with the image scaling factor determined above.
* **Calculate the corrected reading** -- For each segment, the mean pixel value is determined for both the open and DMLC image.
  These values are used to determine the corrected reading: :math:`R_{corr}`.
* **Calculate sample and segment ratios** -- Once the images are normalized to themselves, the sample values of the DMLC
  field are divided by their corresponding open field values.
* **Calculate segment deviations** -- Segment deviation is then calculated once all the corrected readings are determined.
  The average absolute deviation can also be calculated.

**Post-Analysis**

* **Test if segments pass tolerance** -- Each segment is checked to see if it was within the specified tolerance. If any samples
  fail, the whole test is considered failing.


API Documentation
-----------------

.. autoclass:: pylinac.vmat.VMAT
    :no-show-inheritance:

.. autoclass:: pylinac.vmat.Settings
    :no-show-inheritance:

.. autoclass:: pylinac.vmat.SegmentHandler
    :no-show-inheritance:

.. autoclass:: pylinac.vmat.Segment
