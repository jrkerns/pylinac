
.. _starshot_doc:

========
Starshot
========

.. _star_overview:

Overview
--------

.. automodule:: pylinac.starshot
    :no-members:

Running the Demo
----------------

To run the Starshot demo, create a script or start an interpreter and input:

.. plot::

    from pylinac import Starshot

    Starshot.run_demo()

Results will be printed to the console and a matplotlib figure showing the analyzed starshot image will pop up::

    Result: PASS

    The minimum circle that touches all the star lines has a diameter of 0.381 mm.

    The center of the minimum circle is at 1270.0, 1437.2

Image Acquisition
-----------------

To capture starshot images, film is often used, but a sequence of EPID images can also work for collimator measurements. Pylinac can automatically superimpose the images.
See the literature mentioned in the :ref:`star_overview` for more info on acquisition.

Typical Use
-----------

The Starshot analysis can be run first by importing the Starshot class:

.. code-block:: python

    from pylinac import Starshot

A typical analysis sequence looks like so:

* **Load image(s)** -- Loading film or superimposed EPID DICOM images can be done by
  passing the file path or by using a UI to find and get the file. The code might look like any of the following:

  .. code-block:: python

    star_img = "C:/QA Folder/gantry_starshot.tif"
    mystar = Starshot(star_img)

  Multiple images can be easily superimposed and used; e.g. collimator shots at various angles:

  .. code-block:: python

    star_imgs = ["path/star0.tif", "path/star45.tif", "path/star90.tif"]
    mystar = Starshot.from_multiple_images(star_imgs)

* **Analyze the image** -- After loading the image, all that needs to be done is analyze the image. You may optionally
  pass in some settings:

  .. code-block:: python

    mystar.analyze(radius=0.5, tolerance=0.8)  # see API docs for more parameter info

* **View the results** -- Starshot can print out the summary of results to the console as well as draw a matplotlib image to show the
  detected radiation lines and wobble circle:

  .. code-block:: python

      # print results to the console
      print(mystar.results())
      # view analyzed image
      mystar.plot_analyzed_image()

  Additionally, the data can be accessed through a convenient :class:`~pylinac.starshot.StarshotResults` class
  which comes in useful when using pylinac through an API or for passing data to other scripts/routines.

  .. code-block:: python

    # return a dataclass with introspection
    data = mystar.results_data()
    data.tolerance_mm
    data.passed
    ...

    # return as a dict
    data_dict = mystart.results_data(as_dict=True)
    data_dict["passed"]
    ...

  Each subplot can be plotted independently as well:

  .. code-block:: python

      # just the wobble plot
      mystar.plot_analyzed_subimage("wobble")
      # just the zoomed-out plot
      mystar.plot_analyzed_subimage("whole")

  Saving the images is also just as easy:

  .. code-block:: python

      mystar.save_analyzed_image("mystar.png")

  You may also save to PDF:

  .. code-block:: python

      mystar.publish_pdf("mystar.pdf")

Analysis Parameters
-------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      See :meth:`~pylinac.starshot.Starshot.analyze` for details.

   .. tab-item:: RadMachine
      :sync: radmachine

      * **Image dots per inch**: The DPI of the image. Only needed for TIFF images and only if the resolution tag is
        not in the TIFF tags.
      * **Source to film distance**: The distance from the radiation source to the film **in mm**. Only needed for TIFF images.
      * **Normalized Radius**: The radius of the circle used to detect the spokes. 0 is the center of the spoke intersection
        and 1 is the extent of the "shortest" detected spoke. E.g. the spokes may be partially cut off from the image edge.
      * **Tolerance**: The isocenter diameter tolerance in mm.
      * **Normalized minimum peak height**: The percentage minimum height a peak must be to be considered a valid peak. A lower value catches
        radiation peaks that vary in magnitude (e.g. different MU delivered or gantry shot), but could also pick up noise.
        If necessary, lower value for gantry shots and increase for noisy images.

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
* Markings like pin pricks or marker dots can cause wild pixel values and analysis errors and should be cropped out if possible.

**Pre-Analysis**

* **Ground** -- The image is "grounded" by setting the lowest pixel value to 0. This is done to avoid issues with
  high-background pixel values.
* **Check image inversion** -- The image is checked for proper inversion using histogram analysis.
* **Set algorithm starting point** -- If the user provided a ``start_point`` to the ``analyze`` method,
  this value is used. Otherwise, the start point is automatically found by examining the central 1/3 of the image and taking the max of that image window along each axis and finding the
  center of the full-width, 80%-max of each profile.

**Analysis**

* **Extract circle profile** -- A circular profile is extracted from the image centered around the starting point
  and at the radius given.
* **Find spokes** -- The circle profile is analyzed for peaks (either simple peaks or FWHM peaks depending on the ``fwhm`` parameter).
  An even number of spokes (1 for each side; e.g. 3 collimator angles should produce 6
  spokes, one for each side of the CAX) and at least 6 spokes (3 angles) must be found or an error is raised. If ``recursive`` is set to true, this process
  is repeated for the radii range 0.1-0.95 and from 0.05-0.95 ``min_peak_height``. If no combination of these
  parameters produces a valid set of spokes, an error is raised.
* **Match peaks** -- Peaks are matched to their counterparts opposite the CAX to compose a line using a simple peak number offset.
  After matching, the distances from the lines to the start point is measured. If the distance is greater than 10mm of
  any line, an error is raised and follows the same logic regarding ``recursive`` as the previous step. This is to
  avoid issues specifically with gantry starshots. This can also be corrected ahead of time by using a lower ``min_peak_height``.

  .. figure:: images/bad_line_match.png

      The above image has an even number of spokes, but has not caught the low-dose end of the gantry spokes on the right.


* **Find wobble** -- Starting at the initial starting point, a Nelder-Mead gradient method is utilized to find the
  point of minimum distance to all lines. The wobble is assumed to be found if the wobble diameter is <2mm and is <10 mm from the
  starting point. If not, the same logic as the previous steps is followed regarding the ``recursive`` parameter.

**Post-Analysis**

* **Check if passed** -- Once the wobble is calculated, it is tested against the tolerance given, and passes if below the
  tolerance. If the image carried a pixel/mm conversion ratio, the tolerance and result are in mm, otherwise they
  will be in pixels.

.. _interpreting-starshot-results:

Interpreting Results
--------------------

This section explains what is returned in the ``results_data`` object.
This is also the same information that is given in the RadMachine results
section.

* ``tolerance_mm`` -- The tolerance used for the analysis in mm.
* ``circle_radius_mm`` -- The radius of the minimum circle that touches all the star lines in mm.
* ``circle_diameter_mm`` -- The diameter of the minimum circle that touches all the star lines in mm.
* ``circle_center_x_y`` -- The center position of the minimum circle in pixels.
* ``passed`` -- Whether the analysis passed or failed.

Troubleshooting
---------------

First, check the general :ref:`general_troubleshooting` section, especially if an image won't load. Specific to the starshot
analysis, there are a few things you can do.

* **Set recursive to True** - This easy step in :meth:`~pylinac.starshot.Starshot.analyze` allows pylinac to search for a reasonable
  wobble even if the conditions you passed don't for some reason give one.
* **Make sure the center of the star is in the central 1/3 of the image** - Otherwise, pylinac won't find it.
* **Make sure there aren't egregious artifacts** - Pin pricks can cause wild pixel values; crop them out if possible.
* **Set ``invert`` to True** - While right most of the time, it's possible the inversion checker got it wrong. This would
  look like peak locations in the "valley" regions of the image. If so, pass ``invert=True`` to the ``analyze`` method.

Benchmarking the Algorithm
--------------------------


With the image generator module we can create test images to test the starshot algorithm on known results. This is useful to isolate what is or isn't working
if the algorithm doesn't work on a given image and when commissioning pylinac.

Perfect shot
^^^^^^^^^^^^

.. note::

    Due to the rounding of pixel positions of the star lines an absolutely perfect (0.0000mm wobble) is not achievable. The uncertainty of the algorithm is ~0.05mm.

Let's create a perfect irradiation of a starshot pattern:

.. plot::

    from scipy import ndimage

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, FilteredFieldLayer, AS1200Image, RandomNoiseLayer


    star_path = 'perfect_starshot.dcm'
    as1200 = AS1200Image()
    for _ in range(6):
        as1200.add_layer(FilteredFieldLayer((270, 5), alpha=0.5))
        as1200.image = ndimage.rotate(as1200.image, 30, reshape=False, mode='nearest')
    as1200.add_layer(GaussianFilterLayer(sigma_mm=3))
    as1200.generate_dicom(file_out_name=star_path)

    # analyze it
    star = pylinac.Starshot(star_path)
    star.analyze()
    print(star.results())
    star.plot_analyzed_image()

with an output of::

    Result: PASS

    The minimum circle that touches all the star lines has a diameter of 0.045 mm.

    The center of the minimum circle is at 639.5, 639.5

Note that there is still an identified wobble of ~0.045mm due to pixel position rounding of the generated image star lines.
The center of the star is dead on at 639.5 (AS1200 image of shape 1278 and going to the middle of the pixel).

We can also evaluate the effect of changing the radius:

.. plot::

    from scipy import ndimage

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, FilteredFieldLayer, AS1200Image, RandomNoiseLayer


    star_path = 'perfect_starshot.dcm'
    as1200 = AS1200Image()
    for _ in range(6):
        as1200.add_layer(FilteredFieldLayer((270, 5), alpha=0.5))
        as1200.image = ndimage.rotate(as1200.image, 30, reshape=False, mode='nearest')
    as1200.add_layer(GaussianFilterLayer(sigma_mm=3))
    as1200.generate_dicom(file_out_name=star_path)

    # analyze it
    star = pylinac.Starshot(star_path)
    star.analyze(radius=0.6)  # radius changed
    print(star.results())
    star.plot_analyzed_image()

which results in::

    Result: PASS

    The minimum circle that touches all the star lines has a diameter of 0.036 mm.

    The center of the minimum circle is at 639.5, 639.5

The center hasn't moved but we do have a diameter of ~0.03mm now. Again, this is a limitation of both the algorithm and image generation.

Offset
^^^^^^

We can also generate an offset starshot:

.. note::

    This image is completely generated and depending on the angle and number of spokes, this result may change due to the fragility of rotating the image.

.. plot::

    from scipy import ndimage

    import pylinac
    from pylinac.core.image_generator import GaussianFilterLayer, FilteredFieldLayer, AS1200Image, RandomNoiseLayer


    star_path = 'offset_starshot.dcm'
    as1200 = AS1200Image()
    for _ in range(6):
        as1200.add_layer(FilteredFieldLayer((270, 5), alpha=0.5, cax_offset_mm=(1, 1)))
        as1200.image = ndimage.rotate(as1200.image, 60, reshape=False, mode='nearest')
    as1200.add_layer(GaussianFilterLayer(sigma_mm=3))
    as1200.generate_dicom(file_out_name=star_path)

    # analyze it
    star = pylinac.Starshot(star_path)
    star.analyze()
    print(star.results())
    star.plot_analyzed_image()

with an output of::

    Result: FAIL

    The minimum circle that touches all the star lines has a diameter of 1.035 mm.

    The center of the minimum circle is at 637.8, 633.3

Note that we still have the 0.035mm error from the algorithm uncertainty but that we have caught the 1mm offset appropriately.

.. _star_apidoc:

API Documentation
-----------------


.. autoclass:: pylinac.starshot.Starshot
    :members:

.. autopydantic_model:: pylinac.starshot.StarshotResults

.. autoclass:: pylinac.starshot.StarProfile
    :members:

.. autoclass:: pylinac.starshot.Wobble
    :members:

.. autoclass:: pylinac.starshot.LineManager
    :members:
