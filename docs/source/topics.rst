
======
Topics
======

.. _image_loading::

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
under the `metadata` attribute, so if want to access them, they're there.

Pixel Data & Inversion
^^^^^^^^^^^^^^^^^^^^^^

This is the most common issue when dealing with image analysis. The inversion, meaning the pixel value to radiation fluence relationship,
of pylinac images used to be a simple imcompliment, meaning inverting the data while respecting the bit ranges, since
most images' raw pixel data was inverted. However, to handle newer EPID images that included more and better pixel relationships,
this has changed in v3.0.

.. note:: The axiom for pylinac (for v3.0+) is that higher pixel values == more radiation == lighter/whiter display

Assigned pixel values now have the following logic:

If the image has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_,
`Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052`_ and the `Pixel Intensity Relationship Sign <https://dicom.innolitics.com/ciods/rt-image/rt-image/00281041>`_
attributes, all of them are applied with a simple linear correction: :math:`P_corrected = Sign * Slope * P_raw + Intercept`
Images from newer linac platforms appear more likely to have this attribute.

If the image only has the `Rescale Slope <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281053>`_ and
`Rescale Intercept <https://dicom.innolitics.com/ciods/ct-image/ct-image/00281052`_ but not the relationship tag then it is applied as:
:math:`P_corrected = Slope * P_raw + Intercept`. This is the most common scenario encountered to date.

.. note:: It is possible that the slope has a negative value which is implicitly applying a relationship and would be equivalent to the first case, however, older images often have a simple positive slope relationship.

If the image does not have these two tags, then an imcompliment is applied: :math:`new_array = -old_array + max(old_array) + min(old_array)`.
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


.. _contrast:

Contrast
--------

Contrast is used in the catphan and planar imaging modules. There are two contrasts that are evaluated: high contrast
and low contrast. High contrast is also called spatial resolution, and refers to the ability of the device to resolve
high contrast objects that are abutting. This is usually measured with line pairs or a high-contrast point. Low contrast
referes to the ability of the device to measure differences between two similarly-attenuating materials. The materials
and regions need not be abutting as for high contrast.

Depending on who you ask/read, there are multiple definitions of contrast. For high contrast, this is less contentious than
low contrast. We describe here the equations used or offered in pylinac to calculate contrast.

High contrast
^^^^^^^^^^^^^

High contrast calculations are performed by analyzing multiple ROIs and calculating the maximum and minimum pixel value from each ROI.
An ROI is used for each high contrast region (e.g. each line pair region). The contrast is first calculated, then normalized.
The high contrast calculation uses the Michelson contrast, aka visibility. See here for more comparisons: https://en.wikipedia.org/wiki/Display_contrast

.. math:: \frac{ \frac{I_{max} - I_{min}}{I_{max} + I_{min}}}{\max{\left( \frac{I_{max} - I_{min}}{I_{max} + I_{min}}\right)}}

where :math:`I = {1, ..., n}` line pair ROIs.

Low contrast
^^^^^^^^^^^^

Low contrast calculations are also performed by analyzing multiple ROIs, but each ROI has only one value: the mean pixel value.
These pixel values are the compared to a reference ROI. However, that comparison is different depending on who you ask.
Previously, pylinac gave only the Michelson contrast as the low contrast option. However, there are now multiple options available.

.. note:: The combination of low contrast and ROI size is handled in the next section. Do not confuse low contrast with visibility/perception.


For all below :math:`I` is the given ROI and :math:`Ref` is the reference ROI.

Michelson (default; good choice)

.. math:: \frac{I_{mean} - R_{mean}}{I_{mean} + R_{mean}}

Weber

.. math:: \frac{I_{mean} - R_{mean}}{I_{mean}}

Ratio

.. math:: \frac{I_{mean}}{R_{mean}}

.. _visibility:

Visibility
^^^^^^^^^^

Visibility is the ability for humans to detect signal against noise. Visibility is a component of low contrast detectability.
Typically, low contrast is evaluated irrespective of the size of the object. However, as a phantom like Las Vegas shows,
a large-sized object with small contrast might be seen, but a small-sized object of the same contrast might not. This
is referred to as visibility. Visibility in pylinac follows the `Rose <https://www.osapublishing.org/josa/abstract.cfm?uri=josa-38-2-196>`_ model, which is the following formula:

.. math:: Visibility(I) = Contrast(I) * \sqrt{Area(I) * DQE(I)} = Contrast(I) * \frac{\sqrt{\pi * I_{radius}^2}}{I_{std}}

where contrast is an option from the low contrast methods.

.. note::
     What is meant by "noise" is unclear in the literature. Technically, it was meant to be the detective quantum efficiency.
     For simplicity and easy of understanding, the standard deviation does just fine.

.. note::
    Pylinac ROIs are smaller than that actual size of the contrast ROI on the phantom. Uncertainty in the phantom detection
    algorithm means that the ROIs must be smaller to allow tolerance in the algorithm. Thus, visibility is a very specific
    number that depends on the size of the **sampling** ROI.

Contrast-to-noise ratio
^^^^^^^^^^^^^^^^^^^^^^^

The contrast to noise ratio (CNR) is defined as follows:

.. math:: CNR(I) = \frac{Contrast(I)}{noise(I)} = \frac{Contrast(I)}{stdev(I)}

where contrast is an option from the low contrast methods.


Modulation Transfer Function (MTF)
----------------------------------

The MTF is used in CBCT and planar imaging metrics to describe high-contrast characteristics of the imaging system.
An excellent introduction is here: https://www.edmundoptics.com/knowledge-center/application-notes/optics/introduction-to-modulation-transfer-function/
In pylinac, MTF is calculated using equation 3 of the above reference:

.. math:: contrast = \frac{I_max - I_min}{I_max + I_min}

Then, all the contrasts are normalized to the largest one, resulting in a normalized MTF or rMTF (relative).
Pylinac only reports rMTF values. This is the first of two inputs. The other is the line pair spacing. The spacing
is usually provided by the phantom manufacturer. The rMTF is the plotted against the line pair/mm values. Also from
this data the MTF at a certain percentage (e.g. 50%) can be determined in units of lp/mm.

However, it's important to know what :math:`I_max` and :math:`I_min` means here. For a line pair set, each bar and space-between
is one contrast value. Thus, one contrast value is calculated for each bar/space combo. For phantoms with areas of the
same spacing (e.g. the Leeds), all bars and spaces are the same and thus we can use an area-based ROI for the input to
the contrast equation.
