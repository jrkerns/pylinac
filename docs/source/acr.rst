
============
ACR Phantoms
============

Overview
--------

The ACR module provides routines for automatically
analyzing DICOM images of the ACR CT 464 phantom.
It can load a folder or zip file of images, correcting for
translational and rotational offsets.

Phantom reference information is drawn from the
`ACR solution article <https://accreditationsupport.acr.org/support/solutions/articles/11000053945-overview-of-the-ct-phantom>`_
and the analysis is drawn from the
`ACR testing article <https://accreditationsupport.acr.org/support/solutions/articles/11000056197-acr-ct-phantom-scanning-instructions>`_.

Typical Use
-----------

ACR CT analysis follows a similar pattern of load/analyze/output as the rest of the library.
Unlike the CatPhan analysis, customization is not a goal, as the phantom and analysis
are much more well-defined. I.e. there's less of a use case for custom phantoms in this
scenario.

To use the ACR analysis, import the class:

.. code-block:: python

    from pylinac import CT464

And then load, analyze, and view the results:

* **Load images** -- Loading can be done with a directory or zip file:

  .. code-block:: python

    acr_ct_folder = r"C:/CT/ACR/Sept 2021"
    ct = CT464(acr_ct_folder)

  or load from zip:

  .. code-block:: python

    acr_ct_zip = r"C:/CT/ACR/Sept 2021.zip"
    ct = CT464.from_zip(acr_ct_zip)

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

Advanced Use
------------

Using ``results_data``
^^^^^^^^^^^^^^^^^^^^^^

Using the acr module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.acr.CT464.results_data` method
which returns a :class:`~pylinac.acr.ACRCTResult` instance.

.. note::
    While the pylinac tooling may change under the hood, this object should remain largely the same and/or expand.
    Thus, using this is more stable than accessing attrs directly.

Continuing from above:

.. code-block:: python

    data = mycbct.results_data()
    data.catphan_model
    data.ctp404.measured_slice_thickness_mm
    # and more

    # return as a dict
    data_dict = mycbct.results_data(as_dict=True)
    data_dict['ctp404']['measured_slice_thickness_mm']
    ...

 API Documentation
------------------


.. autoclass:: pylinac.acr.CT464
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
