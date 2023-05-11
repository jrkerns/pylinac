
.. _mtf_topic:

Modulation Transfer Function
----------------------------

The modulation transfer function (MTF) is used in CBCT and planar imaging metrics to describe high-contrast characteristics of the imaging system.
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
