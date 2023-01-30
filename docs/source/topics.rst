
.. _topics:

======
Topics
======

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

.. _xim-images:

XIM images
----------

Images ending in ``.xim`` are generally produced by a Varian TrueBeam or newer linac. They are images with additional
tags. Unfortunately, they are written in binary into a custom format so using a typical image library will not work.

The binary file specification appears to be unofficial, but it does work. You can find the spec `here <https://bitbucket.org/dmoderesearchtools/ximreader/raw/4900d324d5f28f8b6b57752cfbf4282b778a4508/XimReader/xim_readme.pdf>`__
which comes from this repo: https://bitbucket.org/dmoderesearchtools/ximreader/src/master/

.. warning::

    Rant ahead.

The XIM images used a custom compression format. Why they chose to use a custom format is beyond me. Moreso, the
format they chose was that of a PNG algorithm. So, XIM images are just PNG images but with a custom lookup table
and property tags. A TIFF format would've worked just as well. It's possible this is security by obscurity or NIH syndrome.

Loading an XIM image
^^^^^^^^^^^^^^^^^^^^

To load an XIM images use the :class:`~pylinac.core.image.XIM` class:

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    # plot the image
    xim_img.plot()

    # see the XIM properties
    print(xim_img.properties)

Reconstructing the image pixels is relatively slow (~1s for AS1200 image) thanks to the custom compression format,
so if you are only searching through the properties you can skip reconstructing the pixels. Skipping the
pixels and only reading the properties is relatively fast (order of milliseconds):

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_files = [r"C:\TDS\H12345\QA\image.xim", ...]
    files_to_analyze = []
    for file in my_xim_files:
        # will load relatively fast
        xim_img = XIM(file, read_pixels=False)
        if xim_img.properties['AcquisitionMode'] == 'Highres':
            files_to_analyze.append(file)

    # now load the pixel data only for the files we're interested in
    for file in files_to_analyze:
        xim_img = XIM(file)
        # image is available, do what you want
        xim_img.plot()

An XIM has all the utility methods other pylinac image do, so use this to your advantage:

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    # process
    xim_img.crop(pixels=30)
    xim_img.filter()
    xim_img.fliplr()
    ...

.. _export-xim:

Exporting images
^^^^^^^^^^^^^^^^

Exporting .xim images is easy. The PNG format is recommended because its ~1/2 the size of the xim image and will
also include the properties. PNG images can usually be viewed easily across many devices and OSs and also loads very fast.

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    xim_img.save_as('myxim.png')
    # saved to PNG!

.. _reading-exported-xim:

Reading exported images
^^^^^^^^^^^^^^^^^^^^^^^

To load the image in python you can use any library that reads PNG. Pillow is recommended.
Opening these files are usually very fast (order of milliseconds), so
if you plan on doing research or analysis of a large number of .xim images, it may be worth it
to export to PNG en masse and then perform the analysis.

.. code-block:: python

    import numpy as np
    import PIL.Image
    import matplotlib.pyplot as plt

    xim_img = PIL.Image.open('myxim.png')

    # numpy array of the pixels
    xim_array = np.asarray(xim_img)

    # plot it
    plt.imshow(xim_array)
    plt.show()

To read the properties of a xim file that was saved to PNG we may to have to load from strings.
PNG tags are all strings, and some xim properties are arrays or numbers. In order to
easily save it, we convert them all to strings. In order to get the native datatype
if it wasn't originally a string is to use ``json``:

.. code-block:: python

    import json
    import PIL.Image

    xim_img = PIL.Image.open('myxim.png')

    system_version = xim_img.info['AcquisitionSystemVersion']
    # "2.7.304.16" already a string so no change needed

    couch_lat = xim_img.info['CouchLat']
    # '100.39021332'  it's a string even though it looks like a number
    # convert to the original type:
    couch_lat_num = float(couch_lat)

    # MLCs are a list; we need json
    mlc_a_string = xim_img.info['MLCLeafsA']
    # '[20.6643, 20.6992, ...]'
    mlc_a_list = json.loads(mlc_a_string)
    # now it's a normal list: [20.6643, 20.6992, ...]

.. _contrast:

Contrast
--------

Contrast is used in the catphan and planar imaging modules. There are two contrasts that are evaluated: high contrast
and low contrast. High contrast is also called spatial resolution, and refers to the ability of the device to resolve
high contrast objects that are abutting. This is usually measured with line pairs or a high-contrast point. Low contrast
refers to the ability of the device to measure differences between two similarly-attenuating materials. The materials
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

.. _low_contrast_topic:

Low contrast
^^^^^^^^^^^^

Low contrast calculations are also performed by analyzing multiple ROIs, but each ROI has only one value: the median pixel value.
These pixel values are compared to a reference ROI. However, that comparison is different depending on who you ask.
Previously, pylinac gave only the Michelson contrast as the low contrast option. However, there are now multiple options available.

.. note:: The combination of low contrast and ROI size is handled in the next section. Do not confuse low contrast with visibility/perception.


For all below :math:`I` is the given ROI and :math:`R` is the reference ROI.

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
Typically, low contrast is evaluated irrespective of the size of the object. However, as a phantom like the Las Vegas or CatPhan 515 module shows,
a large-sized object with small contrast might be seen, but a small-sized object of the same contrast might not. This
is referred to as visibility. Visibility in pylinac is a derivation of the `Rose <https://www.osapublishing.org/josa/abstract.cfm?uri=josa-38-2-196>`_ model,
defined here as:

.. math:: Visibility(I) = Contrast(I) * \sqrt{Area(I) * DQE(I)} = Contrast(I) * \frac{\sqrt{\pi * I_{radius}^2}}{I_{std}}

where contrast is an option from the :ref:`low contrast methods <low_contrast_topic>` and :math:`\pi * I_{radius}^2` is the area of the ROI, which is assumed to be circular.

.. note::
     What is meant by "noise" is unclear in the literature. Technically, it was meant to be the detective quantum efficiency (DQE).
     For simplicity and ease of understanding, the standard deviation works.

.. note::
    Pylinac ROIs are smaller than that actual size of the contrast ROI on the phantom. Uncertainty in the phantom detection
    algorithm means that the ROIs must be smaller to allow a small localization tolerance in the algorithm. Thus, visibility is a very specific
    number that depends on the size of the **sampling** ROI.

Contrast-to-noise ratio
^^^^^^^^^^^^^^^^^^^^^^^

The contrast to noise ratio (CNR) is defined as follows:

.. math:: CNR(I) = \frac{Contrast(I)}{noise(I)} = \frac{Contrast(I)}{stdev(I)}

where contrast is an option from the low contrast methods.

.. _mtf_topic:

Modulation Transfer Function (MTF)
----------------------------------

The MTF is used in CBCT and planar imaging metrics to describe high-contrast characteristics of the imaging system.
An excellent introduction is here: https://www.edmundoptics.com/knowledge-center/application-notes/optics/introduction-to-modulation-transfer-function/
In pylinac, MTF is calculated using equation 3 of the above reference:

.. math:: contrast = \frac{I_{max} - I_{min}}{I_{max} + I_{min}}

Then, all the contrasts are normalized to the largest one, resulting in a normalized MTF or rMTF (relative).
Pylinac only reports rMTF values. This is the first of two inputs. The other is the line pair spacing. The spacing
is usually provided by the phantom manufacturer. The rMTF is the plotted against the line pair/mm values. Also from
this data the MTF at a certain percentage (e.g. 50%) can be determined in units of lp/mm.

However, it's important to know what :math:`I_{max}` and :math:`I_{min}` means here. For a line pair set, each bar and space-between
is one contrast value. Thus, one contrast value is calculated for each bar/space combo. For phantoms with areas of the
same spacing (e.g. the Leeds), all bars and spaces are the same and thus we can use an area-based ROI for the input to
the contrast equation.
