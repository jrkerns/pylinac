
=========================
VMAT module documentation
=========================

Overview
--------

.. automodule:: pylinac.vmat
    :no-members:


.. note::
    There are two classes in the VMAT module: ``DRGS`` and ``DRMLC``. Each have the exact same methods.
    Anytime one class is used here as an example, the other class can be used the same way.

Running the Demos
-----------------

For this example we will use the DRGS class:

.. plot::

    from pylinac import DRGS
    DRGS.run_demo()

Results will be printed to the console and a figure showing both the Open field and MLC field image will pop up::

    Dose Rate & Gantry Speed
    Test Results (Tol. +/-1.5%): PASS
    Max Deviation: 1.01%
    Absolute Mean Deviation: 0.459%

Image Acquisition
-----------------

If you want to perform these specific QA tests, you'll need DICOM plan files that control the linac precisely to deliver the test fields.
These can be downloaded from my.varian.com. Once logged in, search for RapidArc and you should see two items called "RapidArc QA Test Procedures
and Files for TrueBeam"; there will be a corresponding one for C-Series. Use the RT Plan files and follow the instructions, not including the assessment procedure,
which is the point of this module. Save & move the VMAT images to a place you can use pylinac.

Typical Use
-----------

The VMAT QA analysis follows what is specified in the Varian RapidArc QA tests and assumes your tests will run the exact same way.
Import the appropriate class:

.. code-block:: python

    from pylinac import DRGS, DRMLC

The minimum needed to get going is to:

* **Load images** -- Loading the EPID DICOM images into your VMAT class object can be done by passing the file paths,
  passing a ZIP archive, or passing a URL:

  .. code-block:: python

      open_img = "C:/QA Folder/VMAT/open_field.dcm"
      dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"
      mydrgs = DRGS(image_paths=(open_img, dmlc_img))  # use the DRMLC class the exact same way

      # from zip
      mydrmlc = DRMLC.from_zip(r'C:/path/to/zip.zip')

      # from a URL
      mydrgs = DRGS.from_url('http://myserver.org/vmat.zip')


  Finally, if you don't have any images, you can use the demo ones provided:

  .. code-block:: python

     mydrgs = DRGS.from_demo_images()
     mydrmlc = DRMLC.from_demo_images()

* **Analyze the images** -- Once the images are loaded, tell the class to analyze the images.
  See the Algorithm section for details on how this is done. Tolerance can also be passed and has a default value of 1.5%:

  .. code-block:: python

      mydrgs.analyze(tolerance=1.5)

* **View/Save the results** -- The VMAT module can print out the summary of results to the console as well as draw a matplotlib image to show where the
  segments were placed and their values:

  .. code-block:: python

      # print results to the console
      print(mydrgs.results())
      # view analyzed images
      mydrgs.plot_analyzed_image()

  .. plot::
      :include-source: false

      import pylinac

      pylinac.DRGS.run_demo()

  PDF reports can also be generated:

  .. code-block:: python

    myvmat.publish_pdf('drgs.pdf')


Accessing Data
--------------

As with most other modules, the raw data can be examined and use as part of a larger project, e.g. QATrack+.
Perusing the documentation will show what is available, but a quick list is shown here by example:

.. code-block:: python

    mydrgs.avg_abs_r_deviation  # float
    mydrgs.avg_r_deviation  # float; usually artificially low due to positive and negative R values
    mydrgs.max_r_deviation  # float; regardless of sign
    mydrgs.dmlc_image  # DicomImage instance
    mydrgs.open_image  # DicomImage instance
    mydrgs.passed  # bool
    mydrgs.r_devs  # numpy array of all R_deviation values
    mydrgs.segments[1].r_corr  # the 1st segment's R ratio
    mydrgs.segments[0].passed  # whether the 0th segment passed based on the tolerance of ``analyze()``


Algorithm
---------

The VMAT analysis algorithm is based on the Varian RapidArc QA Test Procedures for C-Series and Truebeam. Two tests
(besides Picket Fence, which has its own module) are specified. Each test takes 10x0.5cm samples, each corresponding to
a distinct section of radiation. A corrected reading of each segment is made, defined as:
:math:`M_{corr}(x) = \frac{M_{DRGS}(x)}{M_{open}(x)} * 100`. The reading deviation of each segment is calculated as:
:math:`M_{deviation}(x) = \frac{M_{corr}(x)}{M_{corr}} * 100 - 100`, where :math:`M_{corr}` is the average of all segments.

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

* **Calculate sample boundaries** -- Segment positions are always the same within the image. The x-positions are based on the
  FWHM of the detected field. This allows for old and new style tests that have an x-offset from each other.
  These values are then scaled with the image scaling factor determined above.
* **Calculate the corrected reading** -- For each segment, the mean pixel value is determined for both the open and DMLC image.
  These values are used to determine the corrected reading: :math:`M_{corr}`.
* **Calculate sample and segment ratios** -- The sample values of the DMLC
  field are divided by their corresponding open field values.
* **Calculate segment deviations** -- Segment deviation is then calculated once all the corrected readings are determined.
  The average absolute deviation is also calculated.

**Post-Analysis**

* **Test if segments pass tolerance** -- Each segment is checked to see if it was within the specified tolerance. If any samples
  fail, the whole test is considered failing.


API Documentation
-----------------

.. autoclass:: pylinac.vmat.DRGS

.. autoclass:: pylinac.vmat.DRMLC

.. autoclass:: pylinac.vmat.VMATBase

.. autoclass:: pylinac.vmat.Segment
