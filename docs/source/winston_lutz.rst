
=================================
Winston-Lutz module documentation
=================================

Overview
--------

.. automodule:: pylinac.winston_lutz
    :no-members:

Running the Demo
----------------

To run the Winston-Lutz demo, create a script or start an interpreter session and input:

.. code-block:: python

    from pylinac import WinstonLutz
    WinstonLutz.run_demo()

Results will be printed to the console and a figure showing the zoomed-in images will be generated::

    Winston-Lutz Analysis
    =================================
    Number of images: 17
    Maximum 2D CAX->BB distance: 1.23mm
    Median 2D CAX->BB distance: 0.69mm
    Shift to iso: facing gantry, move BB: RIGHT 0.36mm; OUT 0.36mm; DOWN 0.20mm
    Gantry 3D isocenter diameter: 1.05mm (9/17 images considered)
    Maximum Gantry RMS deviation (mm): 1.03mm
    Maximum EPID RMS deviation (mm): 1.31mm
    Gantry+Collimator 3D isocenter diameter: 1.11mm (13/17 images considered)
    Collimator 2D isocenter diameter: 1.09mm (7/17 images considered)
    Maximum Collimator RMS deviation (mm): 0.79
    Couch 2D isocenter diameter: 2.32mm (7/17 images considered)
    Maximum Couch RMS deviation (mm): 1.23

.. plot::
    :include-source: false

    from pylinac import WinstonLutz

    WinstonLutz.run_demo()

Image Acquisition
-----------------

The Winston-Lutz module will only load EPID images. The images can be from any EPID however, and any SID. To ensure
the most accurate results, a few simple tips should be followed. Note that these are not unique to pylinac; most
Winston-Lutz analyses require these steps:

* The BB should be fully within the field of view.
* The MLC field should be symmetric.
* The BB should be <2cm from the isocenter.

.. _coordinate_space:

Coordinate Space
----------------

.. note::

   In pylinac 2.3, the coordinates changed to be compliant with IEC 61217. Compared to previous versions,
   the Y and Z axis have been swapped. The new Z axis has also flipped which way is positive.


When interpreting results from a Winston-Lutz test, it's important to know the coordinates, origin, etc. Pylinac uses
IEC 61217 coordinate space. Colloquial descriptions are as if standing at the foot of the couch looking at the gantry.

.. image:: images/IEC61217.svg

* **X-axis** - Lateral, or left-right, with right being positive.
* **Y-axis** - Superior-Inferior, or in-out, with sup/in being positive.
* **Z-axis** - Anterior-Posterior, or up-down, with up/anterior being positive.


Typical Use
-----------

Analyzing a Winston-Lutz test is as simple as loading in your images. So, let's import the class:

.. code-block:: python

    from pylinac import WinstonLutz

From here, you can load a directory:

.. code-block:: python

    my_directory = 'path/to/wl_images'
    wl = WinstonLutz(my_directory)

You can also load a ZIP archive with the images in it:

.. code-block:: python

    wl = WinstonLutz.from_zip('path/to/wl.zip')

And that's it! Once loaded you can view images, print the results, or publish a PDF report:

.. code-block:: python

    # plot all the images
    wl.plot_images()
    # plot an individual image
    wl.images[3].plot()
    # save a figure of the image plots
    wl.save_plots('wltest.png')
    # print to PDF
    wl.publish_pdf('mywl.pdf')

If you want to shift the BB based on the results and perform the test again there is a method for that:

.. code-block:: python

    print(wl.bb_shift_instructions())
    # LEFT: 0.1mm, DOWN: 0.22mm, ...

You can also pass in your couch coordinates and the new values will be generated:

.. code-block:: python

    print(wl.bb_shift_instructions(couch_vrt=0.41, couch_lng=96.23, couch_lat=0.12))
    New couch coordinates (mm): VRT: 0.32; LNG: 96.11; LAT: 0.11

.. _using_file_names_wl:


Accessing data
--------------

All the data is easily reachable. For attributes that are inclusive of all the images (e.g. isocenter size), they can be
reached like so:

.. code-block:: python

    wl = WinstonLutz(...)
    wl.gantry_iso_size
    wl.collimator_iso_size
    wl.couch_iso_size
    wl.gantry_coll_iso_size
    wl.bb_shift_vector
    wl.cax2bb_distance('max')  # this is the maximum scalar between CAX and BB across all images

Accessing individual images
---------------------------

Each image can be plotted and otherwise accessed easily:

.. code-block:: python

    wl = WinstonLutz(...)
    # access first image
    wl.images[0]  # these are subclasses of the pylinac.core.image.DicomImage class, with a few special props
    # plot 3rd image
    wl.images[0].plot()  # the plot method is special to the WL module and shows the BB, EPID, and Field CAX.
    # get 2D x/y vector of an image
    wl.images[4].cax2bb_vector  # this is a Vector with a .x and .y attribute. Note that x and y are in respect to the image, not the fixed room coordinates.

Using File Names
----------------

If your linac EPID images do not include axis information (such as Elekta) you can specify it in the file name.
Any and all of the three axes can be defined. If one is not defined and is not in the DICOM tags, it will default to 0.
The syntax to define the axes: "<*>gantry0<*>coll0<*>couch0<*>". There can be any text before, after, or in between each axis definition.
However, the axes numerical value **must** immediately follow the axis name. Axis names are also fixed. The following examples
are valid:

* MyWL-gantry0-coll90-couch315.dcm
* gantry90_stuff_coll45-couch0.dcm
* abc-couch45-gantry315-coll0.dcm
* 01-gantry0-abcd-coll30couch10abc.dcm
* abc-gantry30.dcm
* coll45abc.dcm

The following are invalid:

* mywl-gantry=0-coll=90-couch=315.dcm
* gan45_collimator30-table270.dcm

Using the filenames within the code is done by passing the ``use_filenames=True`` flag to the init method:

.. code-block:: python

    my_directory = 'path/to/wl_images'
    wl = WinstonLutz(my_directory, use_filenames=True)

.. note:: If using filenames any relevant axes must be defined, otherwise they will default to zero. For example,
          if the acquisition was at gantry=45, coll=15, couch=0 then the filename must include both the gantry and collimator
          in the name (<...gantry45...coll15....dcm>). For this example, the couch need not be defined since it is 0.

Image types & output definitions
--------------------------------

The following terms are used in pylinac's WL module and are worth defining.

**Image axis definitions/Image types**
Images are classified into 1 of 6 image types, depending on the position of the axes. The image type is then
used for determining whether to use the image for the given calculation. Image types allow the module to isolate the
analysis to a given axis if needed. E.g. for gantry iso size, as opposed to overall iso size, only the gantry should be moving
so that no other variables influence it's calculation.

* **Reference**: This is when all axes are at value 0 (gantry=coll=couch=0).
* **Gantry**: This is when all axes but gantry are at value 0, e.g. gantry=45, coll=0, couch=0.
* **Collimator**: This is when all axes but collimator are at value 0.
* **Couch**: This is when all axes but the couch are at value 0.
* **GB Combo**: This is when either the gantry or collimator are non-zero but the couch is 0.
* **GBP Combo**: This is where the couch is kicked and the gantry and/or collimator are rotated.

**Analysis definitions**
Given the above terms, the following calculations are performed.

* **Maximum 2D CAX->BB distance (scalar, mm)**: Analyzes all images individually for the maximum 2D distance from rad field center to the BB.
* **Median 2D CAX->BB distance (scalar, mm)**: Same as above but the median.
* **Shift of BB to isocenter (vector, mm)**: The instructions of how to move the BB/couch in order to place the BB at the determined isocenter.
* **Gantry 3D isocenter diameter (scalar, mm)**: Analyzes only the gantry axis images (see above image types). Applies backprojection of the
  CAX in 3D and then minimizes a sphere that touches all the 3D backprojection lines.
* **Gantry+Collimator 3D isocenter diameter (scalar, mm)**: Same as above but also considers Collimator and GB Combo images.
* **[Couch, Collimator] 2D isocenter diameter (scalar, mm)**: Analyzes only the collimator or couch images to determine
  the size of the isocenter according to the axis in question. The maximum distance between any of the points is the isocenter size.
  The couch and collimator are treated separately for obvious reasons. If no
  images are given that rotate about the axis in question (e.g. cardinal gantry angles only) the isocenter size will default to 0.
* **[Maximum, All][Gantry, Collimator, Couch, GB Combo, GBP Combo, EPID] RMS deviation (array of scalars, mm)**: Analyzes the images for the axis in question to determine the overall RMS
  inclusive of all 3 coordinate axes (vert, long, lat). I.e. this is the overall displacement as a function of the axis in question.
  For EPID, the displacement is calculated as the distance from image center to BB for all images with couch=0. If no
  images are given that rotate about the axis in question (e.g. cardinal gantry angles only) the isocenter size will default to 0.


Algorithm
---------

The Winston-Lutz algorithm is based on the works of `Winkler et al`_, `Du et al`_, and `Low et al`_.
Winkler found that the collimator and couch iso could be found using a minimum optimization
of the field CAX points. They also found that the gantry isocenter could by found by "backprojecting"
the field CAX as a line in 3D coordinate space, with the BB being the reference point. This method is used to find the
gantry isocenter size.

Low determined the geometric transformations to apply to 2D planar images to calculate the shift to apply to the BB.
This method is used to determine the shift instructions. Specifically, equations 6 and 9.

.. note::

    If doing research, it is very important to note that Low implicitly used the "Varian" coordinate system.
    This is an old coordinate system and any new Varian linac actually uses IEC 61217. However, because the
    gantry and couch definitions are different, the matrix definitions are technically incorrect when using
    IEC 61217. Pylinac assumes the images are in IEC 61217 scale and will internally convert it to varian scale
    to be able to use Low's equations. The gantry is also reversed in Varian scale, but this can be compensated
    by simply reversing the sign of the Z axis.

The algorithm works like such:

**Allowances**

* The images can be acquired with any EPID (aS500, aS1000, aS1200) at any SID.
* The BB does not need to be near the real isocenter to determine isocenter sizes,
  but does affect the 2D image analysis.

**Restrictions**

    .. warning:: Analysis can fail or give unreliable results if any Restriction is violated.

* The BB must be fully within the field of view.
* The BB must be within 2.0cm of the real isocenter.
* The images must be acquired with the EPID.
* The linac scale should be IEC 61217.

**Analysis**

* **Find the field CAX** -- The spread in pixel values (max - min) is divided by 2, and any pixels above
  the threshold is associated with the open field. The pixels are converted to black & white and
  the center of mass of the pixels is assumed to be the field CAX.

* **Find the BB** -- The image is converted to binary based on pixel values *both* above the 50% threshold as above,
  and below the upper threshold. The upper threshold is an iterative value, starting at the image maximum value,
  that is lowered slightly when the BB is not found. If the binary image has a reasonably circular ROI, the BB is
  considered found and the pixel-weighted center of mass of the BB is considered the BB location.

.. note:: Strictly speaking, the following aren't required analyses, but are explained for fullness and clarity.

* **Backproject the CAX for gantry images** -- Based on the vector of the BB to the field CAX and the gantry angle,
  a 3D line projection of the CAX is constructed. The BB is considered at the origin. Only images where the
  couch was at 0 are used for CAX projection lines.

* **Determine gantry isocenter size** - Using the backprojection lines, an optimization function is run
  to minimize the maximum distance to any line. The optimized distance is the isocenter radius.

* **Determine collimator isocenter size** - The maximum distance between any two field CAX locations
  is calculated for all collimator images.

* **Determine couch isocenter size** - Instead of using the BB as the non-moving reference point, which is now
  moving with the couch, the Reference image (gantry = collimator = couch = 0) CAX location is the reference. The
  maximum distance between any two BB points is calculated and taken as the isocenter size.

.. note::
    Collimator iso size is always in the plane normal to the gantry, while couch iso size is always in
    the x-z plane.

.. _Winkler et al: http://iopscience.iop.org/article/10.1088/0031-9155/48/9/303/meta;jsessionid=269700F201744D2EAB897C14D1F4E7B3.c2.iopscience.cld.iop.org
.. _Du et al: http://scitation.aip.org/content/aapm/journal/medphys/37/5/10.1118/1.3397452
.. _Low et al: https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.597475

API Documentation
-----------------

.. autoclass:: pylinac.winston_lutz.WinstonLutz

.. autoclass:: pylinac.winston_lutz.ImageManager

.. autoclass:: pylinac.winston_lutz.WLImage
