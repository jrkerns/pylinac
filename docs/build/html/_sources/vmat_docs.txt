
=========================
VMAT module documentation
=========================

Overview
========

.. automodule:: pylinac.vmatqa

Running the Demo
================

To run the VMAT demo, run vmat.py, or create a script and run::

    from pylinac import VMAT
    VMAT().run_demo_drgs()
    # alternatively, you can run the MLC Speed demo by:
    VMAT().run_demo_mlcs()

Results will be printed to the console and a figure showing both the Open field and MLC field image will pop up.
#TODO: show screenshots of demo


Image Acquisition
=================

If you want to perform these specific QA tests, you'll need DICOM plan files that control the linac precisely to deliver the test fields.
These can be downloaded from my.varian.com, but conveniently, they are included in this package #TODO: add link to files.

Typical Use
===========

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


* **Analyze the images** -- This is where pylinac does its work. Once everything is set up, tell VMAT to analyze the images. See the
  Algorithm section for details on how this is done::

    # analyze!
    myvmat.analyze()

* **View the results** -- VMAT can print out the summary of results to the console as well as draw a matplotlib image to show where the
  ROIs were taken and their values::

      # print results to the console
      print(myvmat.get_string_results())
      # view analyzed images
      myvmat.show_img_results()

For future reference, you can actually set the test and tolerance as arguments to the analyze method for even shorter code.

Algorithm
=========

The VMAT analyze algorithm works like such:

**Assumptions**

* The tests were delivered using the DICOM RT plan files provided by Varian which follow the layout of Jorgensen et al.
* The images were acquired with the EPID.

**Pre-Analysis**

* *Determine image scaling* -- ROI determination is based on offsets from the center pixel of the image. However,
  some physicists use 150 cm SID and others use 100 cm, and others can use a clinical setting that may be between the two. To account for
  this, the SID is determined and then ROI scaling factors are determined to be able to produce properly-sized ROIs.
* *Calculate ROI regions* -- Once the scaling is determined, the ROI regions are calculated. Because the Jorgensen tests are always the
  same in terms of where the radiation gets delivered, these values are hardcoded. These values are scaled with the scaling factor
  determined above.
* *ROI values extracted & evaluated* --

API Documentation
=================

.. autoclass:: pylinac.vmatqa.vmat.VMAT
    :members:
