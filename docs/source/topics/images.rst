
.. _image_loading:

Images
------

Pylinac deals nearly exclusively with DICOM image data. Film has been actively avoided where possible because of 1)
the increased use and technological advances of EPIDs. EPID data also contains useful tags that give contextual information
about the acquisition (unless you use Elekta). And 2) film images tend to be much messier in general; they often have
markings on them such as a pin prick, marker writing to identify the image, or flash on the edges of the image where
the scanner and film edge did not line up.

How data is loaded
^^^^^^^^^^^^^^^^^^

Pylinac uses the excellent pydicom library to load DICOM images. The pydicom dataset is actually stored in pylinac images
under the ``metadata`` attribute, so if want to access them, they're there.

Pixel Data & Inversion
^^^^^^^^^^^^^^^^^^^^^^

This is the most common issue when dealing with image analysis. The inversion, meaning the pixel value to radiation fluence relationship,
of pylinac images used to be a simple imcompliment, meaning inverting the data while respecting the bit ranges, since
most images' raw pixel data was inverted. However, to handle newer EPID images that included more and better pixel relationships,
this has changed in v3.0.

.. note:: The axiom for pylinac (for v3.0+) is that higher pixel values == more radiation == lighter/whiter display

Assigned pixel values now have the following logic:

If the image has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_,
`Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052>`_ and the `Pixel Intensity Relationship Sign <https://dicom.innolitics.com/ciods/rt-image/rt-image/00281041>`_
attributes, all of them are applied with a simple linear correction: :math:`P_{corrected} = Sign * Slope * P_P{raw} + Intercept`
Images from newer linac platforms appear more likely to have this attribute.

If the image only has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_ and
`Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052>`_ but not the relationship tag then it is applied as:
:math:`P_{corrected} = Slope * P_{raw} + Intercept`. This is the most common scenario encountered to date.

.. note:: It is possible that the slope has a negative value which is implicitly applying a relationship and would be equivalent to the first case, however, older images often have a simple positive slope relationship.

If the image does not have these two tags, then an imcompliment is applied: :math:`new array = -old array + max(old array) + min(old array)`.
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
