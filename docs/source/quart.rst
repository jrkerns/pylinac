.. _quart:

=====
Quart
=====

.. versionadded:: 3.2

Overview
--------

The Quart module provides routines for automatically
analyzing DICOM images of the Quart DVT phantom typically used with the Halcyon linac system.
It can load a folder or zip file of images, correcting for
translational and rotational offsets.

.. versionadded:: 3.2

.. warning::

    These algorithms have only a limited amount of testing data and results should be scrutinized.
    Further, the algorithm is more likely to change in the future when a more robust test suite is built up.
    If you'd like to submit data, enter it `here <https://forms.gle/RBR5ubFvjogE9iC67>`_.

Typical Use
-----------

The Quart phantom analysis follows a similar pattern of load/analyze/output as the rest of the library.
Unlike the CatPhan analysis, customization is not a goal, as the phantoms and analyses
are much more well-defined. I.e. there's less of a use case for custom phantoms in this
scenario.

To use the Quart analysis, import the class:

.. code-block:: python

    from pylinac import QuartDVT
    from pylinac.quart import QuartDVT  # equivalent import

And then load, analyze, and view the results:

* **Load images** -- Loading can be done with a directory or zip file:

  .. code-block:: python

    quart_folder = r"C:/CT/Quart/Sept 2021"
    quart = QuartDVT(quart_folder)


  or load from zip:

  .. code-block:: python

    quart_folder = (
        r"C:/CT/Quart/Sept 2021.zip"  # this contains all the DICOM files of the scan
    )
    quart = QuartDVT.from_zip(quart_folder)

* **Analyze** -- Analyze the dataset:

  .. code-block:: python

    quart.analyze()

* **View the results** -- Reviewing the results can be done in text or dict format
  as well as images:

  .. code-block:: python

    # print text to the console
    print(quart.results())
    # view analyzed image summary
    quart.plot_analyzed_image()
    # view images independently
    quart.plot_images()
    # save the images
    quart.save_images()
    # finally, save a PDF
    quart.publish_pdf("myquart.pdf")

Hypersight
----------

.. versionadded:: 3.17

The Hypersight variant of the Quart phantom includes a water ROI in the HU module.
A sister class can be used to also analyze this phantom: :class:`~pylinac.quart.HypersightQuartDVT`
and will include an additional ROI analysis of the water bubble.

The class can be used interchangeably with the normal class and throughout this documentation.

Advanced Use
------------

Adjusting ROI locations
^^^^^^^^^^^^^^^^^^^^^^^

To adjust ROI locations, see the sister section for CT analysis: :ref:`adjusting-roi-locations`.

Using ``results_data``
^^^^^^^^^^^^^^^^^^^^^^

Using the Quart module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.quart.QuartDVT.results_data` method
which returns a :class:`~pylinac.quart.QuartDVTResult` instance.

Continuing from above:

.. code-block:: python

    data = quart.results_data()
    data.hu_module.roi_radius_mm
    # and more

    # return as a dict
    data_dict = quart.results_data(as_dict=True)
    data_dict["hu_module"]["roi_radius_mm"]
    ...

Analysis Parameters
-------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      See :meth:`pylinac.quart.QuartDVT.analyze` for details.

   .. tab-item:: RadMachine
      :sync: radmachine

      * **HU Tolerance**: The tolerance in HU for the phantom materials.
      * **Scaling tolerance**: The tolerance in mm for the scaling of the phantom.
      * **Slice thickness tolerance**: The tolerance in mm for the slice thickness.
      * **CNR Ratio**: The required minimum ratio for the contrast-to-noise to be considered passing.

Algorithm
---------

The Quart algorithm is nearly the same as the :ref:`CBCT Algorithm <cbct-algorithm>`. The
image loading and localization use the same type of logic.

.. _quart-high-res:

High-Resolution
^^^^^^^^^^^^^^^

For high-resolution resolvability, the Quart manual does describe an equation for calculating
the MTF using the line-spread function (LSF) of the phantom edge. For simplicity,
we use the Varian Halcyon IPA document, which outlines a similar logic
with specific measurements of the -700 -> -200 HU distance using a vertical
and horizontal profile.

Within pylinac, to reduce the number of input parameters and also match commissioning
values, these are the values used. The result is the distance in mm from
these two HU values.

.. note::

  The images in pylinac are "grounded", meaning -1000 -> 0. So the actual algorithm
  search values are +300 HU (-700 + 1000) and +800 HU (-200 + 1000).


CNR/SNR
^^^^^^^

While normally the :ref:`contrast <contrast>` algorithm is chosen by the user,
for the Quart phantom it is hardcoded based on the equations in the manual.
Specifically, contrast to noise is defined as:

.. math:: \frac{|Polystyrene - Acrylic|}{Acrylic}

where the values are the median pixel value of the given ROI. Poly
was given as a possible recommendations in the Quart user manual.
Acrylic is the base material of the phantom, i.e. background.

.. note:: The numerator is an absolute value.

The signal to noise is defined as:

.. math:: \frac{Polystyrene + 1000}{\sigma_{Polystyrene}}

where :math:`\sigma` is the standard deviation of the Polystyrene ROI pixel values.
The poly ROI was chosen by us to match the selection for the CNR equation.

Interpreting Results
--------------------

The outcome from analyzing the phantom available in RadMachine or from
``results_data`` is:

* ``phantom_model``: The model of the phantom, e.g. "Quart DVT".
* ``phantom_roll_deg``: The roll of the phantom in degrees.
* ``origin_slice``: The slice number of the origin image.
* ``num_images``: The number of images given in the dataset.
* ``hu_module``: A dictionary of the HU module results with the following
  items:

  * ``offset``: The offset of the module slice in mm from the origin slice.
  * ``measured_slice_thickness_mm``: The measured slice thickness in mm.
  * ``signal_to_noise``: The signal to noise ratio.
  * ``contrast_to_noise``: The contrast to noise ratio.
  * ``roi_settings``: A dictionary of the ROI settings. The keys are the HU material such as ``Acrylic`` with the following
    items:

    * ``value``: The mean value of the ROI in HU.
    * ``angle``: The angle of the ROI in degrees.
    * ``distance``: The distance of the ROI from the center of the phantom in mm.
    * ``radius``: The radius of the ROI in mm.
    * ``distance_pixels``: The distance of the ROI from the center of the phantom in pixels.
    * ``radius_pixels``: The radius of the ROI in pixels.
    * ``angle_corrected``: The angle of the ROI corrected for phantom roll in degrees.

  * ``rois``: A dictionary of ROI results where the key is the name of the material. Each material
    has the following items:

    * ``name``: The name of the material.
    * ``value``: The mean value of the ROI in HU.
    * ``stdev``: The standard deviation of the ROI in HU.
    * ``difference``: The difference in HU from the Acrylic ROI.
    * ``nominal_value``: The nominal value of the material in HU.
    * ``passed``: A boolean indicating if the material HU was within tolerance.

* ``uniformity_module``: A dictionary of the uniformity module results with the following
  items:

  * ``offset``: The offset of the module slice in mm from the origin slice.
  * ``passed``: A boolean indicating if the module passed.
  * ``rois``: A dictionary of ROI results where the key is the name of the material. Each material
    has the following items:

    * ``name``: The name of the material.
    * ``value``: The mean value of the ROI in HU.
    * ``stdev``: The standard deviation of the ROI in HU.
    * ``difference``: The difference in HU from the Acrylic ROI.
    * ``nominal_value``: The nominal value of the material in HU.
    * ``passed``: A boolean indicating if the material HU was within tolerance.

  * ``roi_settings``: A dictionary of the ROI settings. The keys are the HU material such as ``Acrylic`` with the following
    items:

    * ``value``: The mean value of the ROI in HU.
    * ``angle``: The angle of the ROI in degrees.
    * ``distance``: The distance of the ROI from the center of the phantom in mm.
    * ``radius``: The radius of the ROI in mm.
    * ``distance_pixels``: The distance of the ROI from the center of the phantom in pixels.
    * ``radius_pixels``: The radius of the ROI in pixels.
    * ``angle_corrected``: The angle of the ROI corrected for phantom roll in degrees.

* ``geometric_module``: A dictionary containing the following items:

  * ``offset``: The offset of the module slice in mm from the origin slice.
  * ``distances``: A dictionary of the phantom size itself in horizontal and vertical dimensions in mm.
    The keys are ``horizontal mm`` and ``vertical mm``.
  * ``high_contrast_distances``: A dictionary of the high contrast distances in mm. See: :ref:`quart-high-res`.
    The key is the region of the line and the value is the distance in mm.
  * ``mean_high_contrast_distance``: The mean of the high contrast distances in mm.


API Documentation
------------------


.. autoclass:: pylinac.quart.QuartDVT
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.HypersightQuartDVT
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartHUModule
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartUniformityModule
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartGeometryModule
    :inherited-members:
    :members:

.. autopydantic_model:: pylinac.quart.QuartDVTResult

.. autopydantic_model:: pylinac.quart.QuartHUModuleOutput

.. autopydantic_model:: pylinac.quart.QuartUniformityModuleOutput

.. autopydantic_model:: pylinac.quart.QuartGeometryModuleOutput
