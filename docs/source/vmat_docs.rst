
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

    from pylinac import VMAT
    VMAT.run_demo_drgs()
    # alternatively, you can run the MLC Speed demo by:
    VMAT.run_demo_drmlc()

Results will be printed to the console and a figure showing both the Open field and MLC field image will pop up::

    Dose Rate & MLC Speed
    Test Results (Tol. +/-1.5%): PASS
    Max Deviation: 0.437%
    Absolute Mean Deviation: 0.382%

.. plot::

    import pylinac
    pylinac.VMAT.run_demo_drmlc()

Image Acquisition
-----------------

If you want to perform these specific QA tests, you'll need DICOM plan files that control the linac precisely to deliver the test fields.
These can be downloaded from my.varian.com. Once logged in, search for RapidArc and you should see two items called "RapidArc QA Test Procedures
and Files for TrueBeam"; there will be a corresponding one for C-Series. Use the RT Plan files and follow the instructions, not including the assessment procedure,
which is the point of this module. Save & move the images to a place you can use pylinac.

.. _naming_convention:

Naming Convention
-----------------

For VMAT analysis, pylinac needs to know which image is the DMLC image. You can specify this when loading an image and pass the test
type directly. However, if you follow a simple naming convention, you can skip that step.

**If the file has 'open' or 'dmlc' (case-insenitive) in the file name, pylinac will automatically assign the images.
Also, if the filenames have 'drmlc' or 'drgs' in the file names, then the test type is automatically registered as well.**
The following example shows what this would look like::

    # well-named files
    img1 = r'C:/path/drmlc_open.dcm'  # note the drmlc, which is the test type and the open which is the delivery
    img2 = r'C:/path/drmlc_dmlc.dcm'
    # reading these in registers both image type and test type
    vmat = VMAT(images=(img1, img2))
    vmat.analyze() # no need for test type

If your files don't follow this convention it's okay; you just need to specify the image and test types manually. See below
and :meth:`~pylinac.vmat.VMAT.analyze` for more info.

Typical Use
-----------

The VMAT QA analysis follows what is specified in the Varian RapidArc QA tests and assumes your tests will run the exact same way.
Import the class::

    from pylinac import VMAT

The minimum needed to get going is to:

* **Load images** -- Loading the EPID DICOM images into your VMAT class object can be done by passing the file paths,
  passing a ZIP archive, or passing a URL. While the image delivery types can be passed, using an easy :ref:`naming_convention` can simplify the process::

      open_img = "C:/QA Folder/VMAT/open_field.dcm"  # note the 'open'
      dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"  # note the 'dmlc'
      myvmat = VMAT(images=(open_img, dmlc_img))  # the order doesn't matter in this case

  However, if our names weren't so clearly named, we can just pass the delivery types, which assign the images::

      img1 = "C:/QA/vmat1.dcm"  # the open field
      img2 = "C:/QA/vmat2.dcm"  # the dmlc field
      myvmat = VMAT(images=[img1, img2], delivery_types=['open', 'dmlc'])

  A zip file holding both the images can be used, but must follow the :ref:`naming_convention`::

      myvmat = VMAT.from_zip(r'C:/path/to/zip.zip')

  A URL can be passed which must point to a ZIP archive hosted on a server::

     myvmat = VMAT.from_url('http://myserver.org/vmat.zip')

  Finally, if you don't have any images, you can use the demo ones provided::

     drgs_vmat = VMAT.from_demo_images(type='drgs')
     drmlc_vmat = VMAT.from_demo_images(type='drmlc')

* **Analyze the images** -- This is where pylinac does its work. Once the images are loaded, tell VMAT to analyze the images. See the
  Algorithm section for details on how this is done. Tolerance can also be passed, but has a default value of 1.5%::

      myvmat.analyze(test='drmlc', tolerance=1.5)

  .. note::
    The ``test`` parameter does not need to be specified if the :ref:`naming_convention` for tests is followed.

* **View the results** -- The VMAT module can print out the summary of results to the console as well as draw a matplotlib image to show where the
  segments were placed and their values::

      # print results to the console
      print(myvmat.return_results())
      # view analyzed images
      myvmat.plot_analyzed_image()

  .. plot::

    import pylinac
    pylinac.VMAT.run_demo_drmlc()

  The individual plots can also be plotted and saved::

      myvmat.plot_analyzed_subimage('profile')
      myvmat.save_analyzed_subimage('myprofile.png', subimage='profile')

  .. raw:: html
      :file: images/vmat_profiles.html

  PDF reports can also be generated:

  .. code-block:: python

    myvmat.publish_pdf('drmlc.pdf')

Tips & Tricks
-------------

The VMAT test can be customized using the ``tolerance`` parameter to whatever tolerance is needed, whether simply
using a default (1.5), or a clinical value.

Furthermore, the location of the segments can be adjusted. Older QA tests had the segments slightly offset from center.
The new revision of tests is now centered. To account for this, just add an offset to analyze::

    myvmat = VMAT.from_demo_images('drgs')
    myvmat.analyze(x_offset=20)

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

* **Calculate sample boundaries** -- Segment positions are always the same and are determined by offsets from the center pixel given by Varian.
  These values are then scaled with the image scaling factor determined above.
* **Calculate the corrected reading** -- For each segment, the mean pixel value is determined for both the open and DMLC image.
  These values are used to determine the corrected reading: :math:`M_{corr}`.
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

.. autoclass:: pylinac.vmat.Settings

.. autoclass:: pylinac.vmat.SegmentManager

.. autoclass:: pylinac.vmat.Segment
