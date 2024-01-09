
.. _mtf_topic:

Modulation Transfer Function
----------------------------

Peak-Valley MTF
~~~~~~~~~~~~~~~

The modulation transfer function (MTF) is used in CBCT and planar imaging metrics to describe high-contrast characteristics of the imaging system.
An excellent introduction is `here <https://www.edmundoptics.com/knowledge-center/application-notes/optics/introduction-to-modulation-transfer-function/>`__.
In pylinac, MTF is calculated using equation 3 of the above reference, which is also the :ref:`Michelson <michelson>` contrast definition.

.. math:: contrast = \frac{I_{max} - I_{min}}{I_{max} + I_{min}}

Then, all the contrasts are normalized to the largest one, resulting in a normalized MTF or rMTF (relative).
Pylinac only reports rMTF values. This is the first of two inputs. The other is the line pair spacing. The spacing
is usually provided by the phantom manufacturer. The rMTF is the plotted against the line pair/mm values. Also from
this data the MTF at a certain percentage (e.g. 50%) can be determined in units of lp/mm.

However, it's important to know what :math:`I_{max}` and :math:`I_{min}` means here. For a line pair set, each bar and space-between
is one contrast value. Thus, one contrast value is calculated for each bar/space combo. For phantoms with areas of the
same spacing (e.g. the Leeds), all bars and spaces are the same and thus we can use an area-based ROI for the input to
the contrast equation.

Moments-based MTF
~~~~~~~~~~~~~~~~~

The MTF can also be calculated using the moments of the line pair spread function (LPSF).
This algorithm is based on the work of Hander et al [1]_. Specifically, equation 8:

.. math::

   MTF = \frac{\sqrt{2 * (\sigma^{2} - \mu)}}{\mu}

where :math:`\mu` is the mean pixel value of the ROI and :math:`\sigma` is the standard deviation of the ROI pixel values.


.. [1] `Hander et al. <https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.597928>`__ "Rapid objective measurement of gamma camera resolution using statistical moments" (1998).
