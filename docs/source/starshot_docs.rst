
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

To capture starshot images, film is often used, but a sequence of EPID images can also work for collimator measurements. Pylinac can automatically superimpose the images.
See the literature mentioned in the :ref:`star_overview` for more info on acquisition.

Typical Use
-----------

The Starshot analysis can be run first by importing the Starshot class::

    from pylinac.starshot import Starshot

A typical analysis sequence looks like so:

* **Load image(s)** -- Loading film or superimposed EPID DICOM images can be done by
  passing the file path or by using a UI to find and get the file. The code might look like any of the following::

    star_img = "C:/QA Folder/gantry_starshot.tif"
    mystar = Starshot(star_img)

  Or, use a dialog box::

    mystar = Starshot.from_image_UI()

  Multiple images can be easily superimposed and used; e.g. collimator shots at various angles::

    mystar = Starshot.from_multiple_images_UI()

  Or, pass them programatically::

    star_imgs = ['path/star0.tif', 'path/star45.tif', 'path/star90.tif']
    mystar = Starshot.from_multiple_images()

  .. note::
    In previous versions of pylinac, loading images was instance-method based. This behavior has been simplified in favor
    of initialization normalization and adding class-method constructors (``Starshot.from_X``). The reason for this is that
    certain actions should only be allowed until after the image is loaded. Furthermore, loading the image should always be
    the first action of the analysis sequence. By using class constructors, certain pitfalls and errors can be avoided.
    Don't worry though, the old behavior still works.

* **Analyze the image** -- After loading the image, all that needs to be done is analyze the image. You may optionally
  pass in some settings::

    # analyze
    mystar.analyze(radius=50, tolerance=0.8) # see API docs for more parameter info

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

* The image must have at least 6 spokes (3 angles).
* The center of the "star" must be in the central 1/3 of the image.
* The radiation spokes must extend to both sides of the center. I.e. the spokes must not end at the center of the circle.

**Pre-Analysis**

* **Check image inversion** -- The image is checked for proper inversion by summing the image along each axis and then
  effectively finding the point of maximum value. If the point is not in the central 1/3 of the image, it is thought to be inverted.
  This check can be skipped but is enabled by default.
* **Set algorithm starting point** -- Unless the user has manually set the pixel location of the start point,
  it is automatically found by summing the image along each axis and finding the
  center of the full-width, 80%-max of each sum. The maximum value point is also located. Of the two points, the
  one closest to the center of the image is chosen as the starting point.

**Analysis**

* **Extract circle profile** -- A circular profile is extracted from the image centered around the starting point
  and at the radius given.
* **Find spokes** -- The circle profile is analyzed for peaks. Optionally, the profile is reanalyzed to find the center
  of the FWHM. An even number of spokes must be found (1 for each side. E.g. 3 collimator angles should produce 6
  spokes, one for each side of the CAX).
* **Match peaks** -- Peaks are matched to their counterparts opposite the CAX to compose a line using a simple peak number offset.
* **Find wobble** -- Starting at the initial starting point, an evolutionary gradient method is utilized to find the
  point of minimum distance to all lines.

**Post-Analysis**

* **Check if passed** -- Once the wobble is calculated, it is tested against the tolerance given, and passes if below the
  tolerance. If the image carried a pixel/mm conversion ratio, the tolerance and result are in mm, otherwise they
  will be in pixels.

Troubleshooting
---------------

First, check the general :ref:`general_troubleshooting` section, especially if an image won't load. Specific to the starshot
analysis, there are a few things you can do.

* **Set recursive to True** - This easy step in :meth:`~pylinac.starshot.Starshot.analyze` allows pylinac to search for a reasonable
  wobble even if the conditions you passed don't for some reason give one.
* **Make sure the center of the star is in the central 1/3 of the image** - Otherwise, pylinac won't find it.
* **Make sure there aren't egregious artifacts** - Pin pricks can cause wild pixel values; crop them out if possible.


.. _star_apidoc:

API Documentation
-----------------

.. autoclass:: pylinac.starshot.Starshot
    :no-show-inheritance:

.. autoclass:: pylinac.starshot.StarProfile

.. autoclass:: pylinac.starshot.Wobble

.. autoclass:: pylinac.starshot.LineManager

.. autoclass:: pylinac.starshot.Tolerance
