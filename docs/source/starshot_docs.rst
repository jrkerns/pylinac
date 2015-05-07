
.. _starshot_doc:

=============================
Starshot module documentation
=============================

.. _star_overview:

Overview
--------

.. automodule:: pylinac.starshot
    :no-members:

Running the Demo
----------------

To run the Starshot demo, create a script or start an interpreter and input::

    from pylinac.starshot import Starshot
    Starshot().run_demo()

Results will be printed to the console and a matplotlib figure showing the analyzed starshot image will pop up::

    Result: PASS

    The minimum circle that touches all the star lines has a diameter of 0.385 mm.

    The center of the minimum circle is at 1269.3, 1437.8

.. image:: images/starshot_analyzed.png

Image Acquisition
-----------------

To capture starshot images, film is often used, but superimposed EPID images can also work for collimator measurements.
See the literature mentioned in the :ref:`star_overview` for more info.

Typical Use
-----------

The Starshot analysis can be run first by importing the Starshot class::

    from pylinac.starshot import Starshot
    mystar = Starshot()

A typical analysis sequence looks like so:

* **Load image(s)** -- Loading film or superimposed EPID DICOM images into your Starshot class object can be done by
  passing the file path or by using a UI to find and get the file. The code might look like any of the following::

    # set the file path
    star_img = "C:/QA Folder/gantry_starshot.tif"
    # load the images from the file path
    mystar.load_image(star_img)

  or::

    # Load the image using a UI file dialog box
    mystar.load_image_UI()

  Multiple images can be easily superimposed and used; e.g. collimator shots at various angles::

    # Select all the images from a UI box
    mystar.load_multiple_images_UI()

  or::

    # set the file paths from a list
    star_imgs = ['path/star0.tif', 'path/star45.tif', 'path/star90.tif']
    mystar.load_multiple_images()

* **Analyze the images** -- After loading the image, all that needs to be done is analyze the image with a few
  settings passed in::

    # analyze
    mystar.analyze(radius=50) # see API docs for more parameter info

* **View the results** -- Starshot can print out the summary of results to the console as well as draw a matplotlib image to show the
  detected radiation lines and wobble circle (zoom in to see the wobble circle)::

      # print results to the console
      print(mystar.return_results())
      # view analyzed image
      mystar.plot_analyzed_image()

Algorithm
---------

**Allowances**

* The image can be either inversion (radiation is darker or brighter).
* The image can be any size.
* The image can be DICOM (from an EPID) or most image formats (scanned film).
* If multiple images are used, they must all be the same size.

**Restrictions**

    .. warning:: Analysis can fail or give unreliable results if any Restriction is violated.

* The center of the "star" must be in the central 1/3 of the image.
* The radiation spokes must extend to both sides of the center. I.e. the spokes must not end at the center of the circle.

**Pre-Analysis**

* **Check image inversion** -- The image is checked for proper inversion by summing the image along each axis and then
  effectively finding the point of maximum value. If the point is not in the central 1/3 of the image, it is thought to be inverted.
  This check can be skipped but is enabled by default.
* **Set algorithm starting point** -- Unless the user has manually set the pixel location of the start point,
  it is automatically found by summing the image along each axis and finding the
  center of the full-width, 80%-max of each sum.

**Analysis**

* **Extract circle profile** -- A circular profile is extracted from the image centered around the starting point
  and at the radius given.
* **Find spokes** -- The circle profile is analyzed for peaks. Once the peaks are found, the pixels around the peak
  are reanalyzed to find the center of the full-width, half-max. An even number of spokes must be found (1 for each
  side. E.g. 3 collimator angles should produce 6 spokes, one for each side of the CAX).
* **Match peaks** -- Peaks are matched to their counterparts opposite the CAX to compose a line using a simple peak offset.
* **Find wobble** -- Starting at the initial starting point, an evolutionary gradient method is utilized like this:
  for a 3x3 matrix around a starting point, the sum of distances to each line is calculated. If the center matrix value
  is not the lowest value (within a tolerance), the calculation is performed again, starting at the point of lowest value of the previous
  calculation. Thusly, the algorithm moves 1 pixel at a time toward a point of minimum distance to all lines. Once the
  nearest pixel is found, the algorithm switches to sub-pixel precision and repeats.

**Post-Analysis**

* **Check if passed** -- Once the wobble is calculated, it is tested against the tolerance given, and passes if below the
  tolerance.

.. _star_apidoc:

API Documentation
-----------------

.. autoclass:: pylinac.starshot.Starshot

.. autoclass:: pylinac.starshot.StarProfile

.. autoclass:: pylinac.starshot.Wobble


