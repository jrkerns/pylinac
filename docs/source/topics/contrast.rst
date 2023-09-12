
.. _contrast:

Contrast
--------

Contrast is used in the catphan and planar imaging modules. There are two contrasts that are evaluated: high contrast
and low contrast. High contrast is also called spatial resolution, and refers to the ability of the device to resolve
high contrast objects that are abutting. This is usually measured with line pairs or a high-contrast point. Low contrast
refers to the ability of the device to measure differences between two similarly-attenuating materials. The materials
and regions need not be abutting as for high contrast.

Depending on who you ask/read, there are multiple definitions of contrast. For high contrast, this is less contentious than
low contrast. We describe here the equations used or offered in pylinac to calculate contrast. A summary of methods
can be read here: https://en.wikipedia.org/wiki/Display_contrast

High contrast
^^^^^^^^^^^^^

High contrast calculations are performed by analyzing multiple ROIs and calculating the maximum and minimum pixel value from each ROI.
An ROI is used for each high contrast region (e.g. each line pair region). The contrast is first calculated, then normalized.
The high contrast calculation uses the :ref:`Michelson <michelson>` algorithm.

.. math:: \frac{ \frac{I_{max} - I_{min}}{I_{max} + I_{min}}}{\max{\left( \frac{I_{max} - I_{min}}{I_{max} + I_{min}}\right)}}

where :math:`I = {1...n}` line pair ROIs. See also the :ref:`MTF <mtf_topic>` section.

.. _low_contrast_topic:

Low contrast
^^^^^^^^^^^^

Low contrast calculations are also performed by analyzing multiple ROIs, but each ROI has only one value: the median pixel value.
These pixel values are compared to a reference ROI. However, that comparison is different depending on who you ask.
Previously, pylinac gave only the Michelson contrast as the low contrast option. However, there are now multiple options available.

.. important:: The combination of low contrast and ROI size, aka visibility, is handled in the next section. Do not confuse low contrast with visibility/perception.


For all below :math:`I` is the given ROI and :math:`R` is the reference ROI.

.. _michelson:

Michelson
"""""""""

.. note:: For backwards compatibility with older pylinac versions, Michelson is the default contrast algorithm.

Michelson is a good algorithm for contrast and is the pylinac default for low contrast calculations.
It is the only algorithm used for high contrast evaluation. The Michelson contrast can range from 0 to 1.
The official definition is:

.. math:: \frac{I_{max} - I_{min}}{I_{max} + I_{min}}

This is how high-contrast ROIs are evaluated, where the ROI contains both the high and low values together.

When applied to low contrast, there are usually two ROIs, the ROI of the contrast value in question and some reference
value, usually the background. For low-contrast evaluation, the equation is the same, but is restated given the two
individual ROIs:

.. math:: \frac{I_{mean} - R_{mean}}{I_{mean} + R_{mean}}

where :math:`I` is the ROI of the contrast region in question and :math:`R` is the background ROI, usually
placed somewhere within the phantom area that is uniform.

It is primarily used in high-contrast, periodic pattern situations, although it can
be used in most situations.

An example calculation:

.. code-block:: python

    from pylinac.core import contrast

    my_roi = np.array((1.17, 1.31, 1.26, ...))
    c = contrast.michelson(my_roi)

.. seealso::

  `Wikipedia <https://en.wikipedia.org/wiki/Contrast_(vision)#Michelson_contrast>`__

Weber
"""""

The Weber algorithm is generally used when a large image has a small feature and the majority of the image is background.
However, how the "feature" and "background" values are calculated is ambiguous. The Weber contrast value ranges from -1 to infinity.

The official definition is:

.. math:: \frac{I - I_{background}}{I_{background}}

Within pylinac, this is interpreted to be the following:

.. math:: \frac{|I_{mean} - R_{mean}|}{R_{mean}}

where :math:`I` is the ROI of the contrast region in question and :math:`R` is the background ROI, usually
placed somewhere within the phantom area that is uniform.

.. important::

    For historical reasons, the numerator is the absolute difference. This means the range is from
    0 to infinity vs -1 to infinity. The repercussions is that contrast is symmetric. I.e. -1 and +1 both go
    to +1.

An example calculation:

.. code-block:: python

    from pylinac.core import contrast

    feature_value = np.max(my_array)
    background = np.median(my_array)
    c = contrast.weber(feature=feature_value, background=background)

.. seealso::

  `Wikipedia <https://en.wikipedia.org/wiki/Contrast_(vision)#Weber_contrast>`__.

Ratio
"""""

The ratio algorithm is simply the value of interest over the reference or background value.

.. math:: \frac{feature}{reference}

Within pylinac, this is interpreted as:

.. math:: \frac{I_{mean}}{R_{mean}}

where :math:`I` is the ROI of the contrast region in question and :math:`R` is the background ROI, usually
placed somewhere within the phantom area that is uniform.

An example calculation:

.. code-block:: python

    from pylinac.core import contrast

    feature_value = np.max(my_array)
    reference = np.min(my_array)
    c = contrast.ratio(feature=feature_value, reference=reference)

Difference
""""""""""

The difference algorithm is the absolute difference of the feature ROI and the reference or background ROI.
You might prefer this algorithm if you want to have a strictly normal definition of CNR like `this <https://en.wikipedia.org/wiki/Contrast-to-noise_ratio>`__.


.. math:: |feature - background|

.. note::

  The absolute difference is used; i.e. the difference algorithm is symmetric.

Within pylinac, this is interpreted as:

.. math:: I_{mean} - R_{mean}

where :math:`I` is the ROI of the contrast region in question and :math:`R` is the background/reference ROI, usually
placed somewhere within the phantom area that is uniform.

An example calculation:

.. code-block:: python

    from pylinac.core import contrast

    c = contrast.ratio(feature=10, background=5)

.. seealso::

    `Wikipedia <https://en.wikipedia.org/wiki/Contrast-to-noise_ratio>`__.

Root-mean-square
""""""""""""""""

The RMS algorithm is another good general algorithm for evaluating contrast in myriad situations. It is defined as:

.. math:: \sqrt{ \frac{1}{M*N} * \sum_{i=0}^{N-1}\sum_{j=0}^{M-1} (I_{i,j} - \bar{I})^2 }

where an image/array is of size :math:`M` by :math:`N`. :math:`\bar{I}` is the average intensity of the image. :math:`I_{i,j}`
is the element at the :math:`i`-th and :math:`j`-th position within the image array dimensions.

.. warning::

    RMS calculations require the input values to be within the range 0 and 1. You might need
    to normalize your image/array first.

An example calculation:

.. code-block:: python

    from pylinac.core import contrast

    my_roi = np.array((0.34, 0.67, 0.44, ...))
    c = contrast.rms(my_roi)

.. seealso::

    `Wikipedia <https://en.wikipedia.org/wiki/Contrast_(vision)#RMS_contrast>`__

.. _visibility:

Visibility
^^^^^^^^^^

Visibility is the ability for humans to detect signal against noise within a certain context. Visibility is a component of low contrast detectability.
Traditionally, low contrast is evaluated irrespective of the size of the object. However, as a phantom like the Las Vegas or CatPhan 515 module shows,
a large-sized object with small contrast might be seen, but a small-sized object of the same contrast might not. This
is referred to as visibility. Visibility following the `Rose <https://www.osapublishing.org/josa/abstract.cfm?uri=josa-38-2-196>`_ model is
defined as:

.. math:: Visibility \approx C * \sqrt{Area * N}

where :math:`C` is contrast and :math:`N` is the number of photons.

Within pylinac, this equation is interpreted as:

.. math:: Visibility(I, R) = Contrast(I, R) * \sqrt{Area(I) * DQE(I)} = Contrast(I, R) * \frac{\sqrt{\pi * radius^2}}{I_{std}}

where contrast is an option from the :ref:`low contrast methods <low_contrast_topic>` and :math:`\pi * radius^2` is the area of the ROI, which is assumed to be circular.
:math:`I` is the contrast image/array and :math:`R` is the reference image/array.

.. note::
     What is meant by "noise" is unclear in the literature. Technically, it was meant to be the detective quantum efficiency (DQE) which correlates
     to the number of photons counted.
     For simplicity and ease of understanding, the standard deviation, aka noise, of the ROI works as a simple inverse surrogate.
     I.e.

     .. math:: \sqrt{DQE(I)}\approx\sqrt{N}\approx\frac{1}{stdev_{I}}

.. note::
    Pylinac ROIs are smaller than that actual size of the contrast ROI on the phantom. Uncertainty in the phantom detection
    algorithm means that the ROIs must be smaller to allow a small localization tolerance in the algorithm. Thus, visibility is a very specific
    number that depends on the size of the **sampling** ROI.

Contrast-to-noise ratio
^^^^^^^^^^^^^^^^^^^^^^^

The contrast to noise ratio (CNR) is defined as follows:

.. math:: CNR(I) = \frac{Contrast(I)}{noise(I)} = \frac{Contrast(I)}{stdev(I)}

where contrast is an option from the low contrast methods.
If you prefer the `classic definition <https://en.wikipedia.org/wiki/Contrast-to-noise_ratio>`__ of CNR
then use the "Difference" algorithm.
