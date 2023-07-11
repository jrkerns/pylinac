.. _acr:

============
ACR Phantoms
============

Overview
--------

.. versionadded:: 3.2

.. warning::

    These algorithms have only a limited amount of testing data and results should be scrutinized.
    Further, the algorithm is more likely to change in the future when a more robust test suite is built up.
    If you'd like to submit data, enter it `here <https://forms.gle/RBR5ubFvjogE9iC67>`_.

The ACR module provides routines for automatically
analyzing DICOM images of the ACR CT 464 phantom and Large MR phantom.
It can load a folder or zip file of images, correcting for
translational and rotational offsets.

Phantom reference information is drawn from the
`ACR CT solution article <https://accreditationsupport.acr.org/support/solutions/articles/11000053945-overview-of-the-ct-phantom>`_
and the analysis is drawn from the
`ACR CT testing article <https://accreditationsupport.acr.org/support/solutions/articles/11000056197-acr-ct-phantom-scanning-instructions>`_.
MR analysis is drawn from the `ACR Guidance document <https://www.acraccreditation.org/-/media/ACRAccreditation/Documents/MRI/LargePhantomGuidance.pdf?la=en>`_.

.. warning::

    Due to the rectangular ROIs on the MRI phantom analysis,
    rotational errors should be <= 1 degree. Translational errors are still
    accounted for however for any reasonable amount.

Typical Use
-----------

The ACR CT and MR analyses follows a similar pattern of load/analyze/output as the rest of the library.
Unlike the CatPhan analysis, customization is not a goal, as the phantoms and analyses
are much more well-defined. I.e. there's less of a use case for custom phantoms in this
scenario. CT is mostly used here but is interchangeable with the MRI class.

To use the ACR analysis, import the class:

.. code-block:: python

    from pylinac import ACRCT, ACRMRILarge

And then load, analyze, and view the results:

* **Load images** -- Loading can be done with a directory or zip file:

  .. code-block:: python

    acr_ct_folder = r"C:/CT/ACR/Sept 2021"
    ct = ACRCT(acr_ct_folder)
    acr_mri_folder = r"C:/MRI/ACR/Sept 2021"
    mri = ACRMRILarge(acr_mri_folder)

  or load from zip:

  .. code-block:: python

    acr_ct_zip = r"C:/CT/ACR/Sept 2021.zip"
    ct = ACRCT.from_zip(acr_ct_zip)

* **Analyze** -- Analyze the dataset:

  .. code-block:: python

    ct.analyze()

* **View the results** -- Reviewing the results can be done in text or dict format
  as well as images:

  .. code-block:: python

    # print text to the console
    print(ct.results())
    # view analyzed image summary
    ct.plot_analyzed_image()
    # view images independently
    ct.plot_images()
    # save the images
    ct.save_analyzed_image()
    # or
    ct.save_images()
    # finally, save a PDF
    ct.publish_pdf()

.. _choosing-mr-echo-number:

Choosing an MR Echo
-------------------

With MRI, a dual echo scan can be obtained. These can result in a combined DICOM dataset but are distinct
acquisitions. To select between multiple echos, use the ``echo_number`` parameter:

.. code-block:: python

  from pylinac import ACRMRILarge

  mri = ACRMRILarge(...)  # load zip or dir with dual echo image set
  mri.analyze(echo_number=2)
  mri.results()

If no echo number is passed, the first and lowest echo number is selected and analyzed.

.. _customizing-acr-modules:

Customizing MR/CT Modules
-------------------------

To customize aspects of the MR analysis modules, subclass the relevant module and set the attribute in the
analysis class. E.g. to customize the "Slice1" MR module:

.. code-block:: python

  from pylinac.acr import ACRMRILarge, MRSlice1Module


  class Slice1Modified(MRSlice1Module):
      """Custom location for the slice thickness ROIs"""

      thickness_roi_settings = {
          "Top": {"width": 100, "height": 4, "distance": -3},
          "Bottom": {"width": 100, "height": 4, "distance": 2.5},
      }


  # now pass to the MR analysis class
  class MyMRI(ACRMRILarge):
      slice1 = Slice1Modified


  # use as normal
  mri = MyMRI(...)
  mri.analyze(...)

There are 4 modules in ACR MRI Large analysis that can be overridden. The attribute name should stay the same
but the name of the subclassed module can be anything as long as it subclasses the original module:

.. code-block:: python

    class ACRMRILarge:
        # overload these as you wish. The attribute name cannot change.
        slice1 = MRSlice1Module
        geometric_distortion = GeometricDistortionModule
        uniformity_module = MRUniformityModule
        slice11 = MRSlice11PositionModule


    class ACRCT:
        ct_calibration_module = CTModule
        low_contrast_module = LowContrastModule
        spatial_resolution_module = SpatialResolutionModule
        uniformity_module = UniformityModule

Customizing module offsets
--------------------------

Customizing the module offsts in the ACR module is easier than for the CT module.
To do so, simply override any relevant constant like so:

.. code-block:: python

  import pylinac

  pylinac.acr.MR_SLICE11_MODULE_OFFSET_MM = 95

  mri = pylinac.ACRMRILarge(...)  # will use offset above

The options for module offsets are as follows along with their default value:

.. code-block:: python

  # CT
  CT_UNIFORMITY_MODULE_OFFSET_MM = 70
  CT_SPATIAL_RESOLUTION_MODULE_OFFSET_MM = 100
  CT_LOW_CONTRAST_MODULE_OFFSET_MM = 30

  # MR
  MR_SLICE11_MODULE_OFFSET_MM = 100
  MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM = 40
  MR_UNIFORMITY_MODULE_OFFSET_MM = 60



Advanced Use
------------

Using ``results_data``
^^^^^^^^^^^^^^^^^^^^^^

Using the ACR module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.acr.ACRCT.results_data` method
which returns a :class:`~pylinac.acr.ACRCTResult` instance. For MRI this is :meth:`~pylinac.acr.ACRMRILarge.results_data` method
and :class:`~pylinac.acr.ACRMRILargeResult` respectively.

Continuing from above:

.. code-block:: python

    data = ct.results_data()
    data.ct_module.roi_radius_mm
    # and more

    # return as a dict
    data_dict = ct.results_data(as_dict=True)
    data_dict["ct_module"]["roi_radius_mm"]
    ...

API Documentation
-----------------


.. autoclass:: pylinac.acr.ACRCT
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.ACRCTResult
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.CTModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.UniformityModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.SpatialResolutionModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.LowContrastModuleOutput
    :inherited-members:
    :members:


.. autoclass:: pylinac.acr.ACRMRILarge
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.ACRMRIResult
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.MRSlice11ModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.MRSlice1ModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.MRUniformityModuleOutput
    :inherited-members:
    :members:

.. autoclass:: pylinac.acr.MRGeometricDistortionModuleOutput
    :inherited-members:
    :members:
