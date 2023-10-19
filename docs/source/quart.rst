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

Advanced Use
------------

Using ``results_data``
^^^^^^^^^^^^^^^^^^^^^^

Using the Quart module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.quart.QuartDVT.results_data` method
which returns a :class:`~pylinac.acr.QuartDVTResult` instance.

Continuing from above:

.. code-block:: python

    data = quart.results_data()
    data.hu_module.roi_radius_mm
    # and more

    # return as a dict
    data_dict = quart.results_data(as_dict=True)
    data_dict["hu_module"]["roi_radius_mm"]
    ...

Algorithm
---------

The Quart algorithm is nearly the same as the :ref:`CBCT Algorithm <cbct-algorithm>`. The
image loading and localization use the same type of logic.

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

API Documentation
------------------


.. autoclass:: pylinac.quart.QuartDVT
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

.. autoclass:: pylinac.quart.QuartDVTResult
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartHUModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartUniformityModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.quart.QuartGeometryModuleOutput
    :inherited-members:
    :members:
