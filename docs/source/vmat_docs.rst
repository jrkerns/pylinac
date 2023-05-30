
====
VMAT
====

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
      mydrgs = DRGS(
          image_paths=(open_img, dmlc_img)
      )  # use the DRMLC class the exact same way

      # from zip
      mydrmlc = DRMLC.from_zip(r"C:/path/to/zip.zip")

      # from a URL
      mydrgs = DRGS.from_url("http://myserver.org/vmat.zip")


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

    myvmat.publish_pdf("drgs.pdf")

.. _customizing_vmat_analysis:

Customizing the analysis
------------------------

You can alter both the segment size and segment positions as desired.

To change the segment size:

.. code-block:: python

    drgs = DRGS.from_demo_image()
    drgs.analyze(
        ..., segment_size_mm=(10, 150)
    )  # ROI segments will now be 10mm wide by 150mm tall
    # same story for DRMLC

To change the x-positions of the ROI segments or change the number of ROI, use a custom ROI config dictionary and pass it to the ``analyze`` method.

.. code-block:: python

    from pylinac import DRGS, DRMLC

    # note the keys are the names of the ROIs and can be anything you like
    custom_roi_config = {
        "200 MU/min": {"offset_mm": -100},
        "300 MU/min": {"offset_mm": -80},
    }  # add more as needed

    my_drgs = DRGS(...)  # works the same way for DRMLC
    my_drgs.analyze(..., roi_config=custom_roi_config)

.. plot::
    :include-source: false

    from pylinac import DRGS

    # note the keys are the names of the ROIs and can be anything you like
    custom_roi_config = {'200 MU/min': {'offset_mm': -20}, '300 MU/min': {'offset_mm': 20}}

    my_drgs = DRGS.from_demo_images()
    my_drgs.analyze(roi_config=custom_roi_config)
    my_drgs.plot_analyzed_image()

Accessing Data
--------------

.. versionchanged:: 3.0

Using the VMAT module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.vmat.VMATBase.results_data` method
which returns a :class:`~pylinac.vmat.VMATResult` instance.

.. note::
    While the pylinac tooling may change under the hood, this object should remain largely the same and/or expand.
    Thus, using this is more stable than accessing attrs directly.

Continuing from above:

.. code-block:: python

    data = my_drmlc.results_data()
    data.test_type
    data.passed
    # and more

    # return as a dict
    data_dict = my_drmlc.results_data(as_dict=True)
    data_dict["test_type"]
    ...

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

* **Calculate sample boundaries** -- The Segment x-positions are based on offsets from the center of the
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

Benchmarking the Algorithm
--------------------------

With the image generator module we can create test images to test the VMAT algorithm on known results. This is useful to isolate what is or isn't working
if the algorithm doesn't work on a given image and when commissioning pylinac.

.. note::

    The below examples are for the DRMLC test but can equally be applied to the DRGS tests as well.

Perfect Fields
^^^^^^^^^^^^^^

In this example, we generate a perfectly flat set of images.

The script will generate the files, but you can also download them here:
:download:`perfect_open_drmlc.dcm <files/perfect_open_drmlc.dcm>` :download:`perfect_dmlc_drmlc.dcm <files/perfect_dmlc_drmlc.dcm>`.

.. plot::

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, PerfectFieldLayer, AS1200Image

    # open image
    open_path = 'perfect_open_drmlc.dcm'
    as1200 = AS1200Image()
    as1200.add_layer(PerfectFieldLayer(field_size_mm=(150, 110), cax_offset_mm=(0, 5)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.generate_dicom(file_out_name=open_path)

    # DMLC image
    dmlc_path = 'perfect_dmlc_drmlc.dcm'
    as1200 = AS1200Image()
    for offset in (-40, -10, 20, 50):
        as1200.add_layer(PerfectFieldLayer((150, 20), cax_offset_mm=(0, offset)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.generate_dicom(file_out_name=dmlc_path)

    # analyze it
    vmat = pylinac.DRMLC(image_paths=(open_path, dmlc_path))
    vmat.analyze()
    print(vmat.results())
    vmat.plot_analyzed_image()

with output::

    Dose Rate & MLC Speed
    Test Results (Tol. +/-1.5%): PASS
    Max Deviation: 0.0%
    Absolute Mean Deviation: 0.0%

Noisy, Realistic
^^^^^^^^^^^^^^^^

We now add a horn effect and random noise to the data:

.. plot::

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, FilteredFieldLayer, AS1200Image, RandomNoiseLayer

    # open image
    open_path = 'noisy_open_drmlc.dcm'
    as1200 = AS1200Image()
    as1200.add_layer(FilteredFieldLayer(field_size_mm=(150, 110), cax_offset_mm=(0, 5)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.add_layer(RandomNoiseLayer(sigma=0.03))
    as1200.generate_dicom(file_out_name=open_path)

    # DMLC image
    dmlc_path = 'noisy_dmlc_drmlc.dcm'
    as1200 = AS1200Image()
    for offset in (-40, -10, 20, 50):
        as1200.add_layer(FilteredFieldLayer((150, 20), cax_offset_mm=(0, offset)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.add_layer(RandomNoiseLayer(sigma=0.03))
    as1200.generate_dicom(file_out_name=dmlc_path)

    # analyze it
    vmat = pylinac.DRMLC(image_paths=(open_path, dmlc_path))
    vmat.analyze()
    print(vmat.results())
    vmat.plot_analyzed_image()

with output::

    Dose Rate & MLC Speed
    Test Results (Tol. +/-1.5%): PASS
    Max Deviation: 0.0332%
    Absolute Mean Deviation: 0.0257%

Erroneous data
^^^^^^^^^^^^^^

Let's now get devious and randomly adjust the height of each ROI (effectively changing the apparent MLC speed):

.. note::

    Due to the purposely random nature shown below, this exact result is likely not reproducible, nor was it intended to be.
    To get reproducible behavior, use numpy with a seed value.

.. plot::

    import random

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, FilteredFieldLayer, AS1200Image, RandomNoiseLayer

    # open image
    open_path = 'noisy_open_drmlc.dcm'
    as1200 = AS1200Image()
    as1200.add_layer(FilteredFieldLayer(field_size_mm=(150, 110), cax_offset_mm=(0, 5)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.add_layer(RandomNoiseLayer(sigma=0.03))
    as1200.generate_dicom(file_out_name=open_path)

    # DMLC image
    dmlc_path = 'noisy_dmlc_drmlc.dcm'
    as1200 = AS1200Image()
    for offset in (-40, -10, 20, 50):
        as1200.add_layer(FilteredFieldLayer((150, 20), cax_offset_mm=(0, offset), alpha=random.uniform(0.93, 1)))
    as1200.add_layer(GaussianFilterLayer(sigma_mm=2))
    as1200.add_layer(RandomNoiseLayer(sigma=0.03))
    as1200.generate_dicom(file_out_name=dmlc_path)

    # analyze it
    vmat = pylinac.DRMLC(image_paths=(open_path, dmlc_path))
    vmat.analyze()
    print(vmat.results())
    vmat.plot_analyzed_image()

with an output of::

    Dose Rate & MLC Speed
    Test Results (Tol. +/-1.5%): FAIL
    Max Deviation: 2.12%
    Absolute Mean Deviation: 1.13%

API Documentation
-----------------

Main classes
^^^^^^^^^^^^

These are the classes a typical user may interface with.

.. autoclass:: pylinac.vmat.DRGS
    :inherited-members:
    :members:

.. autoclass:: pylinac.vmat.DRMLC
    :inherited-members:
    :members:

.. autoclass:: pylinac.vmat.VMATResult
    :inherited-members:
    :members:

.. autoclass:: pylinac.vmat.SegmentResult
    :inherited-members:
    :members:


Supporting Classes
^^^^^^^^^^^^^^^^^^

You generally won't have to interface with these unless you're doing advanced behavior.

.. autoclass:: pylinac.vmat.VMATBase
    :inherited-members:
    :members:

.. autoclass:: pylinac.vmat.Segment
    :inherited-members:
    :members:
