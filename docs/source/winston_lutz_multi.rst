
=========================
Winston-Lutz Multi-Target
=========================

Overview
--------

The Multi-Target Winston-Lutz (MTWL) is an advanced test meant to measure multiple locations away from isocenter,
typically to represent multi-lesion SRS cases. The MTWL module can analyze images with any number of BBs in any
arrangement. It is generalizable such that new phantom analyses can be created quickly.

.. versionadded:: 3.7

Differences from single-target WL
---------------------------------

.. warning::

    The MTWL algorithm is new and provisional. There are a number of limitations with the algorithm. Hopefully,
    these are removed in future updates. The algorithm is still considered valuable even with these limitations
    which is why it is released.

.. important::

    In a nutshell, the MTWL analyzes BB positions only, whereas vanilla WL provides more machine-related data as well as
    BB position data.

Unlike the single-target WL algorithm (aka "vanilla" WL), there are more limitations to acquisition and
outputs. This should improve over time, but for now you can think of the MTWL as a subset of the
vanilla WL algorithm.

* Utility methods such as loading images are the same.
* Outputs related to the BBs are different.
* BB size is not a parameter but is part of the BB arrangement.
* Single images cannot be analyzed.
* Axis deviations (Gantry wobble, etc) are not available (yet).
* A 3D plot of the BB measured and expected positions can be created.

See the following sections for more info.

* :ref:`mtwl_image_acquisition`
* :ref:`mtwl_phantoms`

Running the Demo
----------------

To run the multi-target Winston-Lutz demo, create a script or start an interpreter session and input:

.. code-block:: python

    from pylinac import WinstonLutzMultiTarget

    WinstonLutzMultiTarget.run_demo()

Results will be printed to the console and a figure showing the zoomed-in images will be generated::

    Winston-Lutz Multi-Target Analysis
    ==================================
    Number of images: 4
    Maximum distance between nominal & measured BB locations: 1.16mm
    Median distance between nominal & measured BB locations: 0.63mm
    Mean distance between nominal & measured BB locations: 0.65mm
    BB deviations (mm)
    ==================
    BB #0: Left 0.00 (0.12); In 30.00 (0.01); Up 0.00 (0.16)
    BB #1: Left 30.00 (0.72); In 15.00 (-0.32); Up 0.00 (0.24)
    BB #2: Left 0.00 (-0.41); In 0.00 (0.22); Up 0.00 (-0.01)
    BB #3: Left 0.00 (0.53); In -30.00 (-0.45); Up 0.00 (-0.28)
    BB #4: Left -30.00 (-0.37); In -50.00 (-0.36); Up 0.00 (-1.04)
    BB #5: Left 0.00 (-0.48); In -70.00 (-0.16); Up 0.00 (0.02)

.. plot::
    :include-source: false

    from pylinac import WinstonLutzMultiTarget

    WinstonLutzMultiTarget.run_demo()

.. _mtwl_image_acquisition:

Image Acquisition
-----------------

The Winston-Lutz module will only load EPID images. The images can be from any EPID however, and any SID. To ensure
the most accurate results the following should be noted.

* All BBs must be visible in all images.
* The couch cannot be rotated.
* The BBs should not occlude each other.

Coordinate Space
----------------

The MTWL algorithm uses the same coordinate system as the vanilla WL. :ref:`coordinate_space`.

Passing a coordinate system
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike the vanilla WL, there is no shift instructions to move the phantom to the optimum location, thus
no coordinate system is passed or used (yet).

.. note::

    This is a target for the MTWL algorithm, so expect this to change in the future.

.. _mtwl_phantoms:

Supported Phantoms
------------------

Currently, only the `MultiMet-WL <https://www.sunnuclear.com/products/multimet-wl-cube>`__ cube from SNC is supported. However, the algorithm is generalized and
can be easily adapted to analyze other phantoms. See :ref:`custom-bb-arrangements`.

Typical Use
-----------

Analyzing a multi-target Winston-Lutz test is simple. First, let's import the class:

.. code-block:: python

    from pylinac import WinstonLutzMultiTarget

From here, you can load a directory:

.. code-block:: python

    my_directory = 'path/to/wl_images'
    wl = WinstonLutzMultiTarget(my_directory)

You can also load a ZIP archive with the images in it:

.. code-block:: python

    wl = WinstonLutzMultiTarget.from_zip('path/to/wl.zip')

Now, analyze it. Unlike the vanilla WL algorithm, we have to pass the BB arrangement to know where the BBs should be in space.
:ref:`Preset phantoms <mtwl_phantoms>` exist, or a custom arrangement can be passed.

.. code-block:: python

    wl.analyze(bb_arrangement=BBArrangement.SNC_MULTIMET)

And that's it! You can now view images, print the results, or publish a PDF report:

.. code-block:: python

    # plot all the images
    wl.plot_images()
    # plot a 3D graph of the BB locations
    wl.plot_locations()
    # plot an individual image
    wl.images[3].plot()
    # save a figure of the image plots
    wl.save_plots('wltest.png')
    # print to PDF
    wl.publish_pdf('mywl.pdf')

Passing in Axis values
----------------------

Passing in the axis values work :ref:`the same <passing-in-axis-values>` as for vanilla WL.

Changing BB detection size
--------------------------

To change the size of BB pylinac is expecting you must change it in the BB arrangement. This
allows phantoms with multiple BB sizes to still be analyzed.

.. _custom-bb-arrangements:

Custom BB Arrangements
----------------------

The MTWL algorithm uses a priori BB arrangements. I.e. you need to know where the BBs **should** exist in
space relative to isocenter. The MTWL algorithm is flexible to accommodate any reasonable arrangement of BBs.

To create a custom arrangement, say for an in-house phantom or commercial phantom not yet supported, define the
BB offsets and size like so. Use negative values to move the other direction:

.. code-block:: python

    my_special_phantom_bbs = [
        {'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 5},  # 5mm BB at iso
        {'offset_left_mm': 30, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 4},  # 4mm BB 30mm to left of iso
        {'offset_left_mm': 0, 'offset_up_mm': -20, 'offset_in_mm': 10, 'bb_size_mm': 5},  # BB DOWN 20mm and in 10mm
        ...  # keep going as needed
        )
    ]

Pass it to the algorithm like so:

.. code-block:: python

    wl = WinstonLutzMultiTarget(...)
    wl.analyze(bb_arrangement=my_special_phantom_bbs)
    wl.plot_locations()
    ...

Algorithm
---------

The MTWL algorithm is loosely based on the vanilla WL algorithm. Each image
is searched for the BBs. The BB locations in the EPID plane are then used
to create backprojections of rays emanating from the linac source. Each
BB will thus have a number of projections equal to the number of images.
The intersection of these projections give the perceived location of the BB.

The algorithm works like such:

**Allowances**

* The images can be acquired with any EPID (aS500, aS1000, aS1200) at any SID.
* The image can have any number of BBs.
* The BBs can be at any 3D location.

**Restrictions**

    .. warning:: Analysis can fail or give unreliable results if any Restriction is violated.

* All BBs must be visible in all images.
* Each BB must be within 5mm of the expected position in x and y in the EPID plane. I.e. it must be <=7mm in scalar distance.
* BBs must not occlude or touch each other in the 2D image.

**Analysis**

* **Find the field CAX** -- The spread in pixel values (max - min) is divided by 2, and any pixels above
  the threshold is associated with the open field. The pixels are converted to black & white and
  the center of mass of the pixels is assumed to be the field CAX.

* **Find every BB** -- The image is converted to binary based on pixel values *both* above the 50% threshold as above,
  and below the upper threshold. The upper threshold is an iterative value, starting at the image maximum value,
  that is lowered slightly when the BB is not found. If the binary image has a reasonably circular ROI,
  is approximately the right size, and is within 2cm of the expected BB position, the BB is
  considered found and the pixel-weighted center of mass of the BB is considered the BB location. This is
  repeated for each BB given in the BB arrangement. Any failure to find a BB will cause an analysis error.

* **Create BB ray traces** -- After the BBs are found on each image, "rays" are created from the linac source through
  the location of the BB. This creates a 3D ray, one for each BB on each image.

* **Determine BB measured position** -- For each BB, the location that minimizes the distance to all rays of that BB
  is considered to be the BB position.

* **Evaluate against nominal position** -- Once the measured BB position is known, both the scalar distance and vector
  from the nominal position to the measured position in 3D space is determined.

Benchmarking the Algorithm
--------------------------

With the image generator module we can create test images to test the WL algorithm on known results. This is useful to isolate what is or isn't working
if the algorithm doesn't work on a given image and when commissioning pylinac. It is common, especially with the WL module,
to question the accuracy of the algorithm. Since no linac is perfect and the results are sub-millimeter, discerning what
is true error vs algorithmic error can be difficult. The image generator module is a perfect solution since it can remove or reproduce the former error.

.. note::

    With the introduction of the MTWL algorithm, so to a multi-target synthetic image generator has been created.

2-BB Perfect Delivery
^^^^^^^^^^^^^^^^^^^^^

Create a perfect set of fields with 1 BB at iso and another 20mm inward.

.. plot::

    import pylinac
    from pylinac.core.image_generator import simulators, layers, generate_winstonlutz_multi_bb

    wl_dir = 'wl_dir'
    generate_winstonlutz_multi_bb(
        simulator=simulators.AS1200Image(),
        field_layer=layers.PerfectFieldLayer,
        dir_out='wl_dir',
        offsets=((0, 0, 0), (0, 0, 20)),
        bb_size_mm=4, field_size_mm=(100, 100))

    arrangement = ({'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 5},
                   {'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 20, 'bb_size_mm': 5})
    wl = pylinac.WinstonLutzMultiTarget(wl_dir)
    wl.analyze(arrangement)
    print(wl.results())
    wl.plot_locations()

which has an output of::

    Winston-Lutz Multi-Target Analysis
    ==================================
    Number of images: 4
    Maximum distance between nominal & measured BB locations: 0.00mm
    Median distance between nominal & measured BB locations: 0.00mm
    Mean distance between nominal & measured BB locations: 0.00mm
    BB deviations (mm)
    ==================
    BB #0: Left 0.00 (0.00); In 0.00 (0.00); Up 0.00 (0.00)
    BB #1: Left 0.00 (-0.00); In 20.00 (-0.00); Up 0.00 (0.00)

As shown, we have perfect results.

Offset BBs
^^^^^^^^^^

Let's now offset both BBs by 1mm to the left:

.. plot::

    import pylinac
    from pylinac.core.image_generator import simulators, layers, generate_winstonlutz_multi_bb

    wl_dir = 'wl_dir'
    generate_winstonlutz_multi_bb(
            simulator=simulators.AS1200Image(),
            field_layer=layers.PerfectFieldLayer,
            dir_out=wl_dir,
            offsets=((1, 0, 0), (1, 0, 20)),  # here's the offset
            bb_size_mm=4, field_size_mm=(100, 100))

    arrangement = ({'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 5},
                   {'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 20, 'bb_size_mm': 5})
    wl = pylinac.WinstonLutzMultiTarget(wl_dir)
    wl.analyze(arrangement)
    print(wl.results())
    wl.plot_locations()

with an output of::

    Winston-Lutz Multi-Target Analysis
    ==================================
    Number of images: 4
    Maximum distance between nominal & measured BB locations: 1.01mm
    Median distance between nominal & measured BB locations: 1.00mm
    Mean distance between nominal & measured BB locations: 1.00mm
    BB deviations (mm)
    ==================
    BB #0: Left 0.00 (-1.01); In 0.00 (-0.00); Up 0.00 (-0.00)
    BB #1: Left 0.00 (-1.00); In 20.00 (0.01); Up 0.00 (-0.00)

Both BBs report a shift of 1mm to the left (the value in parentheses).

SNC MultiMet
^^^^^^^^^^^^

We can even simulate the SNC MultiMet phantom. The ``generate_winstonlutz_multi_bb``
function can take in the same arrangement format as the ``analyze`` method.

Additionally, we add jitter to the synthetic images to give some realism.

.. warning::

    BB size is still set in the image generator via the ``bb_size_mm`` parameter; the arrangement value is not respected.

.. plot::

    import pylinac
    from pylinac.winston_lutz import BBArrangement
    from pylinac.core.image_generator import simulators, layers, generate_winstonlutz_multi_bb

    wl_dir = 'wl_dir'
    generate_winstonlutz_multi_bb(
            simulator=simulators.AS1200Image(),
            field_layer=layers.PerfectFieldLayer,
            dir_out=wl_dir,
            offsets=BBArrangement.SNC_MULTIMET,
            bb_size_mm=4, field_size_mm=(120, 100), jitter_mm=1)

    wl = pylinac.WinstonLutzMultiTarget(wl_dir)
    wl.analyze(BBArrangement.SNC_MULTIMET)
    print(wl.results())
    wl.plot_locations()

with output of::

    Winston-Lutz Multi-Target Analysis
    ==================================
    Number of images: 4
    Maximum distance between nominal & measured BB locations: 0.76mm
    Median distance between nominal & measured BB locations: 0.44mm
    Mean distance between nominal & measured BB locations: 0.49mm
    BB deviations (mm)
    ==================
    BB #0: Left 0.00 (0.09); In 30.00 (0.18); Up 0.00 (-0.73)
    BB #1: Left 30.00 (0.24); In 15.00 (-0.33); Up 0.00 (-0.12)
    BB #2: Left 0.00 (-0.34); In 0.00 (0.27); Up 0.00 (0.10)
    BB #3: Left 0.00 (0.19); In -30.00 (0.05); Up 0.00 (-0.38)
    BB #4: Left -30.00 (-0.18); In -50.00 (-0.02); Up 0.00 (-0.20)
    BB #5: Left 0.00 (0.40); In -70.00 (-0.18); Up 0.00 (0.37)

Here we've purposely added some random noise to the BB locations between 0 and 1mm.

API Documentation
-----------------

.. autoclass:: pylinac.winston_lutz.WinstonLutzMultiTarget
    :members:

.. autoclass:: pylinac.winston_lutz.WinstonLutzResult
    :members:
    :inherited-members:

.. autoclass:: pylinac.winston_lutz.WinstonLutz2DMultiTarget
    :members:
