.. _cheese:

=================
"Cheese" Phantoms
=================

.. versionadded:: 3.9

.. warning::

    These algorithms have only a limited amount of testing data and results should be scrutinized.
    Further, the algorithm is more likely to change in the future when a more robust test suite is built up.
    If you'd like to submit data, enter it `here <https://forms.gle/RBR5ubFvjogE9iC67>`_. Thanks!

The Cheese module provides routines for automatically analyzing
DICOM images of phantoms commonly called "cheese" phantoms, defined
by round phantoms with holes where the user can insert plugs, usually of known density.
The primary use case is performing HU calibration or verification although some plugs
allow for additional functionality such as spatial resolution.
It can load a folder or zip file of images, correcting for
translational and rotational offsets.

Phantoms Supported
------------------

The following phantoms are supported:

* Tomotherapy Cheese :class:`~pylinac.cheese.TomoCheese`
* CIRS Electron Density  :class:`~pylinac.cheese.CIRS062M`

Image Acquisition
-----------------

To be able to properly analyze the phantom, images should be
acquired in such a way that the couch and phantom are separated.
If the phantom it sitting in such a way that it is directly touching the
couch analysis will likely fail. If the phantom does not come with a
support device that keeps it off the table, it should be set on a low-density material such as foam blocks.

Running the Demo
----------------

To run one of the Cheese phantom demos, create a script or start an interpreter and input:

.. plot::

    from pylinac import TomoCheese

    TomoCheese.run_demo()

Results will be also be printed to the console::

    - TomoTherapy Cheese Phantom Analysis -
    - HU Module -
    ROI 1 median: 17.0, stdev: 37.9
    ROI 2 median: 20.0, stdev: 44.2
    ROI 3 median: 23.0, stdev: 36.9
    ROI 4 median: 1.0, stdev: 45.7
    ROI 5 median: 17.0, stdev: 37.6
    ROI 6 median: -669.0, stdev: 39.6
    ROI 7 median: 14.5, stdev: 45.8
    ROI 8 median: 26.0, stdev: 38.6
    ROI 9 median: 653.0, stdev: 47.4
    ROI 10 median: 25.0, stdev: 36.7
    ROI 11 median: 24.0, stdev: 35.3
    ROI 12 median: 102.0, stdev: 46.2
    ROI 13 median: 8.0, stdev: 38.1
    ROI 14 median: -930.0, stdev: 43.8
    ROI 15 median: 23.0, stdev: 36.3
    ROI 16 median: 15.0, stdev: 37.1
    ROI 17 median: -516.0, stdev: 45.1
    ROI 18 median: 448.0, stdev: 38.1
    ROI 19 median: 269.0, stdev: 45.3
    ROI 20 median: 15.0, stdev: 37.9


Typical Use
-----------

The cheese phantom analyses follows a similar pattern of load/analyze/output as the rest of the library.
Unlike the CatPhan analysis, tolerances are not applied and comparison to known values is not the goal.
There are two reasons for this: 1) The plugs are interchangable and thus the reference values are not
necessarily constant. 2) Evaluation against a reference is not the end goal as described in the :ref:`philosophy <philosophy>`.
Thus, measured values are provided; what you do with them is your business.

To use the Tomo Cheese analysis, import the class:

.. code-block:: python

    from pylinac import TomoCheese

And then load, analyze, and view the results:

* **Load images** -- Loading can be done with a directory or zip file:

  .. code-block:: python

    cheese_folder = r"C:/TomoTherapy/QA/September"
    cheese = TomoCheese(cheese_folder)

  or load from zip:

  .. code-block:: python

    cheese_zip = r"C:/TomoTherapy/QA/September.zip"
    cheese = TomoCheese.from_zip(cheese_zip)

* **Analyze** -- Analyze the dataset:

  .. code-block:: python

    cheese.analyze()

* **View the results** -- Reviewing the results can be done in text or dictionary format
  as well as images:

  .. code-block:: python

    # print text to the console
    print(cheese.results())
    # return a dictionary or dataclass
    results = cheese.results_data()
    # view analyzed image summary
    cheese.plot_analyzed_image()
    # save the images
    cheese.save_analyzed_image()
    # finally, save a PDF
    cheese.publish_pdf()

.. _plotting_tomo_density:

Plotting density
----------------

An HU-to-density curve can be plotted if an ROI configuration is passed to the ``analyze`` parameter like so:

.. code-block:: python

    import pylinac

    density_info = {
        "1": {"density": 1.0},
        "3": {"density": 3.05},
    }  # add more as needed. all keys must have a dict with 'density' defined
    tomo = pylinac.TomoCheese(...)
    tomo.analyze(roi_config=density_info)
    tomo.plot_density_curve()  # in this case, ROI 1 and 3 will be plotted vs the stated density

This will plot a simple HU vs density graph.

.. plot::
    :include-source: false

    import pylinac

    density_info = {'1': {'density': 1.0}, '3': {'density': 3.05}}
    tomo = pylinac.TomoCheese.from_demo_images()
    tomo.analyze(roi_config=density_info)
    tomo.plot_density_curve()

.. note::

    The keys of the configuration must be strings matching the ROI number on the phantom. I.e. ``1`` matches to "ROI 1", etc.

.. note::

    Not all ROI densities have to be defined. Any ROI between 1 and 20 can be set.

Adjusting ROI locations
-----------------------

To adjust ROI locations, see the sister section for CT analysis: :ref:`adjusting-roi-locations`.

.. _extending_cheese_phantom:

Extending for other phantoms
----------------------------

While new commercial cheese-like phantoms will continue to be added to this module,
creating new classes is relatively easy. The following steps show how this can be accomplished.

#. Create a new class "module" that inherits from ``CheeseModule``. This
   class contains information about the ROIs, such as the distance and angle away from the center.
   You can use the ``TomoCheeseModule`` as a guide in the source code. An example:

   .. code-block:: python

        from pylinac.cheese import CheeseModule


        class SwissCheeseModule(CheeseModule):
            common_name = "Swiss cheese phantom"
            roi_settings = {  # configuration of each ROI.
                "1": {  # each ROI should have a string key and the following keys
                    "angle": 90,
                    "distance": 45,
                    "radius": 6,
                },
                "2": {
                    "angle": 45,
                    "distance": 80,
                    "radius": 6,
                },
                "3": {...},
            }

   .. note::

        Not all ROIs have to be defined. E.g. if you are only interested in 5 ROIs out of 20 then simply configure those 5.

#. Create a new class that inherits from ``CheesePhantomBase``. This will define the phantom itself:

   .. code-block:: python

        from pylinac.cheese import CheesePhantomBase


        class SwissCheesePhantom(CheesePhantomBase):
            model = "Swiss Cheese Phantom"
            # generally this is just the radius of a normal ROI
            air_bubble_radius_mm = 14
            # This is the radius in mm to a "ring" of ROIs that is used for localization and roll determination.
            # Generally speaking, set the value to the ring that contains the highest ROI HUs.
            localization_radius = 110
            # minimum number of images that should be in the dataset
            min_num_images = 10
            # the radius of the phantom itself
            catphan_radius_mm = 150
            # set this to the module we just created above
            module_class = SwissModule
            # Optional: for the best type inference when using an IDE, set this as well to the new module. Note it's only a type annotation!!
            module: SwissModule

#. Use the class as normal. The base classes contain all the infrastructure code for analysis and plotting.

   .. code-block:: python

        swiss = SwissCheesePhantom("my/swiss/data")
        swiss.analyze()
        swiss.plot_analyzed_image()

Analysis Parameters
-------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      See :meth:`pylinac.cheese.TomoCheese.analyze` for details.

   .. tab-item:: RadMachine
      :sync: radmachine

      * **ROI <N> Density**: The density of the plug in the ROI. This is used to plot the HU-to-density curve.
        The number of inputs depend on the phantom. All densities are optional. If at least one
        density is provided, a density-to-HU curve will be plotted. Not all densities need to be provided. I.e.
        any subset of the ROI densities can be provided.
      * **X adjustment**: A fine-tuning adjustment to the detected x-coordinate of the phantom center. This will move the
        detected phantom position by this amount in the x-direction in mm. Positive values move the phantom to the right.
      * **Y adjustment**: A fine-tuning adjustment to the detected y-coordinate of the phantom center. This will move the
        detected phantom position by this amount in the y-direction in mm. Positive values move the phantom down.
      * **Angle adjustment**: A fine-tuning adjustment to the detected angle of the phantom. This will rotate the phantom by this amount in degrees.
        Positive values rotate the phantom clockwise.
      * **ROI size factor**: A fine-tuning adjustment to the ROI sizes of the phantom. This will scale the ROIs by this amount.
        Positive values increase the ROI sizes. In contrast to the scaling adjustment, this
        adjustment effectively makes the ROIs bigger or smaller, but does not adjust their position.
      * **Scaling factor**: A fine-tuning adjustment to the detected magnification of the phantom. This will zoom the ROIs and phantom outline (if applicable) by this amount.
        In contrast to the roi size adjustment, the scaling adjustment effectively moves the phantom and ROIs
        closer or further from the phantom center. I.e. this zooms the outline and ROI positions, but not ROI size.
      * **Origin slice**: The slice number that corresponds to the slice to analyze.
        This is a fallback mechanism in case the automatic detection fails.

Algorithm
---------

The TomoCheese algorithm leverages a lot of infrastructure from the CatPhan algorithm.
It is not based on a manual.

Allowances
^^^^^^^^^^

* The images can be any size.
* The phantom can have significant translation in all 3 directions.
* The phantom can have significant roll and moderate yaw and pitch.

Restrictions
^^^^^^^^^^^^

.. warning:: Analysis can fail or give unreliable results if any Restriction is violated.

* The phantom cannot touch any edge of the FOV.
* There must be at least one ROI in the "outer" ROI ring that is higher than water/background phantom. This
  has to do with the automatic roll compensation.

  .. note::

        This is not strictly required but will assist in accurate sub-degree roll compensation.

Pre-Analysis
^^^^^^^^^^^^

The pre-analysis is almost exactly the same as the :ref:`CatPhan pre-analysis <cbct_pre-analysis>`.

Analysis
^^^^^^^^

* **Determine image properties** -- Automatic roll compensation is attempted by creating a circular
  profile at the radius of the "outer" ROIs. This profile is then searched for peaks, which correspond
  to high-density plugs that have been inserted. If a peak is not found, no correction is applied and the
  phantom is assumed to be at 0.
  This would occur if all plugs have been filled with water/background plugs. If a peak is found,
  i.e. a plug has been inserted with HU detectably above water, the center of the peak is determined
  which would correspond to the center of the ROI. The distance to the nearest nominal ROI is calculated.
  If the value is <5 degrees, the roll compensation is applied. If the value is >5 degrees, the
  compensation is not applied and the phantom is assumed to be at 0.
* **Measure HU values of each plug** -- Based on the nominal spacing and roll compensation (if applied),
  each plug area is sampled for the median and standard deviation values.

Interpreting Results
^^^^^^^^^^^^^^^^^^^^

The outcome from analyzing the phantom available in RadMachine or from
``results_data`` is:

* ``origin_slice``: The slice index that was used for the ROI analysis.
* ``num_images``: The number of images that were in the passed dataset.
* ``phantom_roll``: The roll of the phantom in degrees.
* ``rois``: A dictionary of ROIs. The key is the ROI number and the value
  of each key contains:

  * ``center_x``: The x-coordinate of the center of the ROI in pixels.
  * ``center_y``: The y-coordinate of the center of the ROI in pixels.
  * ``diameter``: The diameter of the ROI in pixels.
  * ``median``: The median HU value of the ROI.
  * ``std``: The standard deviation of the HU values of the ROI.

* ``roi_<n>``: This is the same thing as an individual result from ``rois``, but the name itself has
  the ROI number appended. I.e. ``rois['11'] == roi_11``. It is redundant information and is the older implementation
  of providing ROI data. Some "cheese" analyses may not have this set of keys. It
  is deprecated due to the variable number of ROIs that can be analyzed.

API Documentation
-----------------


.. autoclass:: pylinac.cheese.TomoCheese
    :inherited-members:
    :members:

.. autoclass:: pylinac.cheese.CIRS062M
    :inherited-members:
    :members:

.. autoclass:: pylinac.cheese.TomoCheeseModule
    :inherited-members:
    :members:

.. autoclass:: pylinac.cheese.CIRSHUModule
    :inherited-members:
    :members:

.. autopydantic_model:: pylinac.cheese.CheeseResult

.. autopydantic_model:: pylinac.cheese.TomoCheeseResult

.. autoclass:: pylinac.cheese.CheesePhantomBase
    :inherited-members:
    :members:
