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
`ACR CT solution article <https://accreditationsupport.acr.org/support/solutions/articles/11000053945-overview-of-the-ct-phantom>`__
and the analysis is drawn from the
`ACR CT testing article <https://accreditationsupport.acr.org/support/solutions/articles/11000056197-acr-ct-phantom-scanning-instructions>`__.
MR analysis is drawn from the `ACR Guidance document <https://www.acraccreditation.org/-/media/ACRAccreditation/Documents/MRI/LargePhantomGuidance.pdf?la=en>`__.

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

MRI Algorithm
-------------

The ACR MR analysis is based on the `official guidance document <https://www.acraccreditation.org/-/media/ACRAccreditation/Documents/MRI/LargePhantomGuidance.pdf?la=en>`__.
Because the guidance document is extremely specific (nice job ACR!) only a few highlights are given here. The guidance is followed as reasonably close as possible.

Allowances
^^^^^^^^^^

* Multiple MR sequences can be present in the dataset.
* The phantom can have significant cartesian shifts.

Restrictions
^^^^^^^^^^^^

* There should be 11 slices per scan (although multiple echo scans are allowed) per the guidance document (section 0.3).
* The phantom should have very little pitch, yaw, or roll (<1 degree).

Analysis
^^^^^^^^

Section 0.4 specifies the 7 tests to perform. Pylinac can perform 6 of these 7. It cannot yet perform the
low-contrast visibility test.

* **Geometric Accuracy** - The geometric accuracy is measured using profiles of slice 5. The only difference
  is that pylinac will use the 60th percentile pixel value of the image as a high-pass filter and then take the FWHM
  of several profiles of this new image. The width between the two pixels defining the FWHM is the diameter.
* **High Contrast** - High contrast is hard to measure for the ACR MRI phantom simply because it does not use line pairs,
  but rather offset dots as well as the qualitative description in the guidance document about how to score these.
  Pylinac measures the high-contrast by sampling a circular ROI on the left ROI (phantom right) set. This is the
  baseline which all other measurements will be normalized to. The actual dot-ROIs are sampled by taking a circular
  ROI of the row-based set and the column-based set. Each row-based ROI is evaluated against the other row-based ROIs.
  The same is done for column-based ROIs. The ROIs use the maximum and minimum pixel values inside the sample ROI.
  No dot-counting is performed.

  .. tip::

    It is suggested to visually perform the contrast measurement manually and compare to pylinac values to establish
    a cross-comparison ratio. After a ratio has been established, the pylinac MTF can be used as the baseline value
    moving forward.

* **Slice thickness** - Slice thickness is measured using the FWHM of two rectangular ROIs. This is very similar
  to the guidance document explanation.

  Slice thickness is defined the same as in the guidance document:

  .. math:: Thickness = 0.2 * \frac{Top * Bottom}{Top + Bottom}

* **Slice Position** - Slice position accuracy is measured very similarly to the manual method described in the document:
  "The display level setting ... should be set to a level roughly half that of the signal in
  the bright, all-water portions of the phantom."
  For each vertical bar, the pixel nearest to the mid-value between min and max of the rectangular ROI is used as the bar position:

  .. math:: position_{bar} = \frac{ROI_{max} - ROI_{min}}{2} + ROI_{min}

  The difference in positions between the bars is the value reported.

* **Uniformity** - Uniformity is measured using a circular ROI at the center of the phantom and ROIs to the top, bottom,
  left, and right of the phantom, very similar to the guidance document.

  The percent integral uniformity (PIU) is defined as:

  .. math:: PIU = 100 * (1 - \frac{high-low}{high+low})

  Instead of using the WL/WW to find the low and high 1cm\ :sup:`2` ROI, pylinac uses the 1st and 99th percentile of pixel values
  inside the central ROI.

  The ghosting ratio is defined the same as the ACR guidance document:

  .. math:: ghosting_{ratio} = |\frac{(top + bottom) - (left + right)}{2*ROI_{large}}|

  where all values are the median pixel values of their respective ROI. The percent-signal ghosting (PSG) is:

  .. math:: PSG = ghosting_{ratio} * 100


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
