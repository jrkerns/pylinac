
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
