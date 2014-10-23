
=========================
VMAT module documentation
=========================

Overview
--------

.. automodule:: pylinac.vmatqa.vmat

Running the Demo
----------------

To run the VMAT demo, run vmat.py, or create a script and run::

    from pylinac import VMAT
    VMAT().run_demo_drgs()
    # alternatively, you can run the MLC Speed demo by:
    VMAT().run_demo_mlcs()

Results will be printed to the console and a figure showing both the Open field and MLC field image will pop up::

    Dose Rate & Gantry Speed
    Test Results (Tol. +/-3.0%): PASS
    Max Positive Deviation: 0.911%
    Max Negative Deviation: -1.636%
    Absolute Mean Deviation: 0.632%

.. image:: /images/vmat_analyzed.png
   :height: 340px
   :width: 800px


Image Acquisition
-----------------

If you want to perform these specific QA tests, you'll need DICOM plan files that control the linac precisely to deliver the test fields.
These can be downloaded from my.varian.com, but conveniently, they are included in this package #TODO: add link to files.

Typical Use
-----------

The VMAT QA analysis follows what is specified in Jorgensen et al. and assumes your tests will run the exact same way. Let us assume
you've made a VMAT object as follows::

    from pylinac import VMAT
    myvmat = VMAT()

The minimum needed to get going is to:

* **Load images** -- Loading the EPID DICOM images into your VMAT class object can be done by passing the file path or by using a UI to
  find and get each file. The code might look like either of the following::

    # set the file path
    open_img = "C:/QA Folder/VMAT/openfield.dcm"
    mlc_img = "C:/QA Folder/VMAT/mlcfield.dcm"
    # load the images from the file path
    myvmat.load_image(open_img, im_type='open')
    myvmat.load_image(mlc_img, im_type='mlc')

    # *OR*

    # Identify the images using a UI
    myvmat.load_image_UI(im_type='open')
    myvmat.load_image_UI(im_type='mlc')


* **Analyze the images** -- This is where pylinac does its work. Once the images are loaded, tell VMAT to analyze the images. See the
  Algorithm section for details on how this is done. The test to run (whether DRGS or DRMLC) needs to be specified. Tolerance
  can also be passed, but has a default value of 3%::

    # analyze!
    myvmat.analyze(test='drmlc', tolerance=3) # define the test and tolerance explicitly

* **View the results** -- The VMAT module can print out the summary of results to the console as well as draw a matplotlib image to show where the
  samples were taken and their values::

      # print results to the console
      print(myvmat.get_string_results())
      # view analyzed images
      myvmat.show_img_results()

Algorithm
---------

The VMAT analyze algorithm is based on an interpretation of Jorgensen et al.'s descriptions. Two terms are used
throughout this module and are worth defining. *Segment* is the area covered while the linac is performing 1 piece of
the test (see Jorgensen Table I and II for more info). For example, the first *segment* of the DRGS test is the area of
radiation exposed while the dose rate is 105 MU/min. *Sample* is the individual area exposed by one MLC pair for any
given *segment*. Thus, a *segment* is composed of a number of individual *samples*.

The algorithm works like such:

**Allowances**

* The images can be acquired at any SID.

**Restrictions**

* The tests must be delivered using the DICOM RT plan files provided by Varian which follow the test layout of Jorgensen et al.
* The images must be acquired with the EPID.
* The MLC model must be a regular Millennium 120 (not HD).

**Pre-Analysis**

* *Determine image scaling* -- Segment determination is based on offsets from the center pixel of the image. However,
  some physicists use 150 cm SID and others use 100 cm, and others can use a clinical setting that may be different
  than either of those. To account for this, the SID is determined and then scaling factors are determined to be
  able to produce properly-sized sample analysis.

**Analysis**

* *Calculate sample boundaries & extract* -- Because the Jorgensen tests are always the
  same in terms of where the radiation gets delivered, these values are hardcoded as offsets from the center pixel.
  These values are then scaled with the image scaling factor determined above. The mean of the pixel values determined
  by the boundaries is saved as the sample value.
* *Normalize images* -- Once the sample values (and thus segments) are determined, both images are normalized.
  Depending on the test, images are normalized by the 4th (DRGS) or mean of 2nd and 3rd (MLCS) segment for each image.
* *Calculate sample and segment ratios* -- Once the images are normalized to themselves, the sample values of the MLC
  field are divided by their corresponding open field values.
* *Calculate sample deviations* -- Sample deviations are calculated using Jorgensen's "deviation" equation. This divides
  each sample by the mean of all the sample values for that MLC pair (i.e. the mean of the segments for that MLC pair).

**Post-Analysis**

* *Test if samples pass tolerance* -- Each sample is checked to see if it was within the specified tolerance. If any samples
  fail, the whole test is considered failing.


API Documentation
-----------------

.. autoclass:: pylinac.vmatqa.vmat.VMAT
    :members:
