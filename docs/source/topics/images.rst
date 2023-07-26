
.. _image_loading:

Images
------

Pylinac deals nearly exclusively with DICOM image data. Film has been actively avoided where possible because of 1)
the increased use and technological advances of EPIDs. EPID data also contains useful tags that give contextual information
about the acquisition (unless you use Elekta). And 2) film images tend to be much messier in general; they often have
markings on them such as a pin prick, marker writing to identify the image, or flash on the edges of the image where
the scanner and film edge did not line up. That being said, pylinac can generally handle DICOM images and general
images (PNG, JPG, etc) relatively well.

The ``image`` module within pylinac is quite powerful and flexible to do arbitrary operations
for use in custom algorithms. For example, images can be loaded easily, filters applied, cropped,
rotated, and more with straightforward methods.

How data is loaded
^^^^^^^^^^^^^^^^^^

Pylinac uses the excellent ``pydicom`` library to load DICOM images. The pydicom dataset is stored in pylinac images
under the ``metadata`` attribute.

For non-DICOM images (JPG, PNG, TIFF, etc), the ``Pillow`` library is used.

To load an image, the easiest way is like so:

.. code-block:: python

    from pylinac import image

    my_dcm = image.load("path/to/my/image.dcm")
    my_dcm.metadata.GantryAngle  # the GantryAngle tag of the DICOM file
    # these won't have the metadata property as they aren't DICOM
    my_tiff = image.load("path/to/my/image.tiff")
    my_jpg = image.load("path/to/my/image.jpg")

See the :func:`~pylinac.core.image.load` function for details. This will return an image-like
object ready for plotting or manipulation. Note that :ref:`XIM <xim-images>` images are handled separately.

We can also test whether a file is image-like without causing an error if it's not:

.. code-block:: python

    from pylinac import image

    is_image = image.is_image("path/to/questionable.file")  # True or False

Image Manipulation
^^^^^^^^^^^^^^^^^^

To manipulate an image, such as cropping, simply run the method. Some examples:

.. code-block:: python

    from pylinac import image

    my_dcm = image.load(...)
    my_dcm.filter(size=0.01, kind="median")
    my_dcm.fliplr()  # flip the image left-right
    my_dcm.ground()  # set minimum value to 0; useful for images with short dynamic range
    my_dcm.crop(pixels=30, edges=("top", "left"))
    my_dcm.normalize()  # normalize values to 1.0
    my_dcm.rot90(n=1)  # rotate the image by 90 degrees
    my_dcm.bit_invert()  # flip the image so that dark is light and light is dark. Useful for EPID images.
    my_dcm.plot()  # plot the image for visualization

These and similar methods are available to all types of images. However, some image types
have additional properties and methods. For a DICOM that is from a linac EPID, we have
a few extras. We need to load it specifically:

.. code-block:: python

    from pylinac import image

    my_linac_dcm = image.LinacDicomImage("path/to/image.dcm")
    my_linac_dcm.cax()  # a Point instance. E.g. (x=550, y=550)
    my_linac_dcm.dpmm()  # the dots/mm at isocenter. Will account for the SID.


TIFF to DICOM
^^^^^^^^^^^^^

Pylinac will often internally convert TIFF images to pseudo-DICOM files so that
the same methods are available as a DICOM. To do so:

.. code-block:: python

    from pylinac import image

    image.tiff_to_dicom(
        tiff_file="path/to/image.tiff",
        dicom_file="my_new_dicom.dcm",
        sid=1000,
        gantry=0,
        coll=0,
        couch=0,
        dpi=280,
    )

We will now have a file in our working directory named ``my_new_dicom.dcm`` that is, for all intents and purposes,
a DICOM file. It can be loaded with ``image.load()`` or ``pydicom`` like any normal DICOM.

Gamma
^^^^^

We can compute the gamma between two arrays or images using :func:`~pylinac.core.image.gamma_2d`:

.. code-block:: python

    import matplotlib.pyplot as plt
    from pylinac import image

    ref = image.load("reference_dicom.dcm")
    eval = image.load("eval_dicom.dcm")

    gamma = image.gamma_2d(
        reference=ref,
        evaluation=eval,
        dose_to_agreement=2,
        distance_to_agreement=3,
        global_dose=True,
        ...,
    )

    # gamma is a numpy array the same size as the reference/eval image
    plt.imshow(gamma)

.. _pixel_inversion:

Pixel Data & Inversion
^^^^^^^^^^^^^^^^^^^^^^

This is the most common issue when dealing with image analysis. The inversion, meaning the pixel value to radiation fluence relationship,
of pylinac images used to be a simple imcompliment, meaning inverting the data while respecting the bit ranges, since
most images' raw pixel data was inverted. However, to handle newer EPID images that included more and better pixel relationships,
this has changed in v3.0.

.. note:: The axiom for pylinac (for v3.0+) is that higher pixel values == more radiation == lighter/whiter display

Image pixel values will proceed through the following conditions. The
first condition that matches will be executed:

* If the ``raw_pixels`` parameter is set to ``True``, no tags will be searched and
  the values from the DICOM file will be used directly. E.g.

  .. code-block:: python

    from pylinac.core import image

    dcm = image.load("my_dcm_file.dcm", raw_pixels=True)
    # OR
    dcm = image.DicomImage("my_dcm_file.dcm", raw_pixels=True)

  .. versionadded:: 3.13

* If the image has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_,
  `Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052>`_ and the `Pixel Intensity Relationship Sign <https://dicom.innolitics.com/ciods/rt-image/rt-image/00281041>`_
  attributes, all of them are applied with a simple linear correction: :math:`P_{corrected} = Sign * Slope * P_{raw} + Intercept`
  Images from newer linac platforms appear more likely to have this attribute.

* If the image only has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_ and
  `Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052>`_ but not the relationship tag then it is applied as:
  :math:`P_{corrected} = Slope * P_{raw} + Intercept`. This is the most common scenario encountered to date.

  .. note:: It is possible that the slope has a negative value which is implicitly applying a relationship and would be equivalent to the first case, however, older images often have a simple positive slope relationship.

* If the image does not have these two tags, then an imcompliment is applied: :math:`new array = -old array + max(old array) + min(old array)`.
  Very old images will likely reach this condition.

.. note::

    If your image appears to be incorrectly inverted, missing tags are likely why.
    Pylinac has parameters to force the inversion of the image if the end result is wrong.
    Furthermore, some modules perform another inversion check at runtime.
    This is mostly historical but was done because some images were always expected to have a certain relationship and
    the tag logic above was not applied consistently (both new and old images were imcomplimented, causing differences).
    For those modules, tags were not used but a simple histogram analysis which expects the irradiated part of the image to be either centrally located
    or most of the image to NOT be irradiated. This is how pylinac historically worked around this issue and got reliable results across image eras.
    However with this new logic, there may be analysis differences for those images. It is more correct to follow the tags but
    for backwards compatibility the module-specific inversion checks remain.
