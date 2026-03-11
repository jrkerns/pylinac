.. _helios:

===============================
GE Helios CT Daily QA Phantom
===============================

Overview
--------

.. versionadded:: 3.X

The Helios module provides routines for automatically
analyzing DICOM images of the **GE Helios CT Daily QA Phantom**
(GE part 2144715). It can load a folder or zip file of images, correcting for
translational offsets.

The phantom is cylindrical with two active sections at fixed scan locations:

* **Section 1** -- Contains a Plexiglass resolution block embedded in water.
  Tests performed: Contrast Scale, High Contrast
  Spatial Resolution (bar patterns).
* **Section 3** -- Uniform water bath. Tests performed: Noise, Uniformity.

Phantom reference information is drawn from the
GE Technical Reference Manual (6060606-1EN Rev.3).

Typical Use
-----------

The GE Helios analysis follows the same load / analyze / output pattern
as the rest of the library.

To use the Helios analysis, import the class:

.. code-block:: python

    from pylinac import GEHeliosCTDaily

And then load, analyze, and view the results:

* **Load images** -- Loading can be done with a directory or zip file:

  .. code-block:: python

    helios_folder = r"C:/CT/Helios/"
    helios = GEHeliosCTDaily(helios_folder)

  or load from zip:

  .. code-block:: python

    helios_zip = r"C:/CT/Helios/helios.zip"
    helios = GEHeliosCTDaily.from_zip(helios_zip)

* **Analyze** -- Analyze the dataset:

  .. code-block:: python

    helios.analyze()

* **View the results** -- Reviewing the results can be done in text or dict format
  as well as images:

  .. code-block:: python

    # print text to the console
    print(helios.results())
    # view analyzed image summary
    helios.plot_analyzed_image()
    # view images independently
    helios.plot_images()
    # save the images
    helios.save_analyzed_image("helios_summary.png")
    # or
    helios.save_images(directory="output/")
    # finally, save a PDF
    helios.publish_pdf("helios_report.pdf")

Customizing Modules
-------------------

To customize aspects of the analysis modules, subclass the relevant module and set the attribute in the
analysis class. E.g. to customize the Contrast Scale module:

.. code-block:: python

    from pylinac.helios import GEHeliosCTDaily, HeliosContrastScaleModule


    class ContrastScaleModified(HeliosContrastScaleModule):
        """Custom location for the contrast scale ROIs"""

        roi_settings = {
            "Plexiglass": {"width": 20, "height": 20, "distance": 35, "angle": -135},
            "Water": {"width": 20, "height": 20, "distance": 75, "angle": 90},
        }


    # now pass to the analysis class
    class MyHelios(GEHeliosCTDaily):
        contrast_scale_module = ContrastScaleModified


    # use as normal
    helios = MyHelios(...)
    helios.analyze(...)

There are 4 modules in ACR MRI Large analysis that can be overridden.
The attribute name should stay the same but the name of the subclassed
module can be anything as long as it subclasses the original module:

.. code-block:: python

    class GEHeliosCTDaily:
        # overload these as you wish. The attribute name cannot change.
        contrast_scale_module = HeliosContrastScaleModule
        high_contrast_module = HeliosHighContrastModule
        noise_uniformity_module = HeliosNoiseUniformityModule

Customizing module offsets
--------------------------

Customizing the module offsets in the Helios module is easier than for the CT module.
To do so, simply override any relevant constant like so:

.. code-block:: python

    import pylinac.helios

    pylinac.helios.SECTION_3_OFFSET_MM = 55

    helios = pylinac.GEHeliosCTDaily(...)  # will use offset above

Advanced Use
------------

Using ``results_data``
^^^^^^^^^^^^^^^^^^^^^^

Using the Helios module in your own scripts? The analysis results they can be accessed
by using the :meth:`~pylinac.helios.GEHeliosCTDaily.results_data` method
which returns a :class:`~pylinac.helios.GEHeliosResult` instance.

Continuing from above:

.. code-block:: python

    data = helios.results_data()
    data.contrast_scale.contrast_difference
    # and more

    # return as a dict
    data_dict = helios.results_data(as_dict=True)
    data_dict["contrast_scale"]["contrast_difference"]

Adjusting ROI Locations
^^^^^^^^^^^^^^^^^^^^^^^

To adjust ROI locations, see the sister section for CT analysis: :ref:`adjusting-roi-locations`.

Analysis Parameters
-------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      See :meth:`pylinac.helios.GEHeliosCTDaily.analyze` for details.

   .. tab-item:: RadMachine
      :sync: radmachine

      To be completed.


Interpreting Results
--------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      The outcome from analyzing the phantom and calling ``.results_data()``
      is a :class:`~pylinac.helios.GEHeliosResult` instance. See the API
      documentation below for details and also :ref:`exporting-results`.

   .. tab-item:: RadMachine
      :sync: radmachine

      The outcome from analyzing the phantom in RadMachine will return an “Entire Result” as follows:


      * ``phantom_model``: The model of the phantom used.
      * ``phantom_roll_deg``: The roll of the phantom in degrees.
      * ``origin_slice``: The slice number of the "origin" slice; for helios this is Section 1.
      * ``num_images``: The number of images in the passed dataset.
      * ``contrast_scale``: The results from the Contrast Scale module, with the following items:

        * ``offset``: Module offset in mm from origin.
        * ``plexiglass_hu``: Mean HU of the Plexiglass ROI.
        * ``plexiglass_stdev``: Standard deviation of the Plexiglass ROI.
        * ``water_hu``: Mean HU of the Water ROI.
        * ``water_stdev``: Standard deviation of the Water ROI.
        * ``contrast_difference``: Plexiglass minus Water.

      * ``high_contrast``: The results from the High Contrast Spatial Resolution module, with the following items:

        * ``offset``: Module offset in mm from origin.
        * ``stdev_per_bar``: Standard deviation for each bar pattern.
        * ``mtf_lp_mm``: Relative MTF at each spatial frequency.

      * ``noise_uniformity``: The results from the Noise & Uniformity module, with the following items:

        * ``offset``: Module offset in mm from origin.
        * ``center_mean_hu``: Center ROI mean HU.
        * ``center_stdev``: Center noise standard deviation.
        * ``edge_rois``: Mean HU values of each edge ROI, keyed by position name.
        * ``uniformity_difference``: Center minus edge mean.

Algorithm
---------

The Helios analysis is based on the Technical Reference Manual.

Analysis
^^^^^^^^

* **Contrast Scale** - The Contrast Scale is measured using two rectangular
  ROIs (10 mm) placed on the origin slice, one over the Plexiglass resolution block and
  one over the surrounding water.
* **High Contrast Spatial Resolution** - Four rectangular ROIs are placed over the
  bar-pattern groups within the Plexiglass block. The standard deviation within each
  ROI quantifies how well the bars are resolved. A relative MTF curve is also
  calculated where each spatial frequency is normalized to the coarsest (1.6 mm) bar
  pattern.
* **Noise & Uniformity** - On the uniform water slice a 25 mm center box ROI
  measures noise (standard deviation). A 15 mm center box ROI and two 15 mm
  edge ROIs measure uniformity.

API Documentation
-----------------
.. autoclass:: pylinac.helios.GEHeliosCTDaily
    :members:
    :inherited-members:

.. autopydantic_model:: pylinac.helios.GEHeliosResult

.. autopydantic_model:: pylinac.helios.HeliosContrastScaleModuleOutput

.. autopydantic_model:: pylinac.helios.HeliosHighContrastModuleOutput

.. autopydantic_model:: pylinac.helios.HeliosNoiseUniformityModuleOutput
