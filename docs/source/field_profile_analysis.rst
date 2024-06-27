.. _field-profile-analysis:

======================
Field Profile Analysis
======================

The modular field analysis is the newer implementation of the :ref:`field_analysis_module` concept.
Colloquially, it is referred to as "v2" of the field analysis module.
It leverages the new :ref:`profiles` framework to provide a more flexible way
to analyze fields.

.. note::

    If you are familiar with the old field analysis module and want to convert, you can read about the differences in
    :ref:`comparison-to-field-analysis-v1`.

Overview
========

The field profile analysis module (``pylinac.field_profile_analysis``) allows a physicist to analyze various metrics of an open radiation
field from extracted profiles. The user can specify the location and width of the profiles
from the image. They can also select various metrics to compute on the profiles as well
as write custom plugins to compute additional metrics as necessary for the specific application
of interest. The use cases and assumptions/constraints are the same as for the :ref:`profiles <profiles>` framework.

Typical Usage
=============

To analyze profiles from an open field, the user should import the class, pass the path to the open beam, and then analyze the profiles from the beam.

.. code-block:: python

  from pylinac import FieldProfileAnalysis, Centering, Normalization, Edge
  from pylinac.metrics.profile import (
      PenumbraLeftMetric,
      PenumbraRightMetric,
      SymmetryAreaMetric,
      FlatnessDifferenceMetric,
  )

  path = r"C:\path\to\open_field.dcm"
  field_analyzer = FieldProfileAnalysis(path)
  field_analyzer.analyze(
      centering=Centering.BEAM_CENTER,
      x_width=0.02,
      y_width=0.02,
      normalization=Normalization.BEAM_CENTER,
      edge_type=Edge.INFLECTION_DERIVATIVE,
      ground=True,
      metrics=(
          PenumbraLeftMetric(),
          PenumbraRightMetric(),
          SymmetryAreaMetric(),
          FlatnessDifferenceMetric(),
      ),
  )
  field_analyzer.plot_analyzed_images(show_grid=True, mirror="beam")

which will output 3 images: the image itself marked with where the profiles were taken
and the X and Y profile plots along with the metrics plotted.

.. plot::
  :include-source: false

  from pylinac import FieldProfileAnalysis, Centering, Normalization
  from pylinac.metrics.profile import PenumbraLeftMetric, PenumbraRightMetric, SymmetryAreaMetric, FlatnessDifferenceMetric

  field_analyzer = FieldProfileAnalysis.from_demo_image()
  field_analyzer.analyze(centering=Centering.BEAM_CENTER, x_width=0.02, y_width=0.02, normalization=Normalization.BEAM_CENTER, ground=True, metrics=(PenumbraLeftMetric(), PenumbraRightMetric(), SymmetryAreaMetric(), FlatnessDifferenceMetric()))
  field_analyzer.plot_analyzed_images(mirror='beam')

Analysis Options
================

Centering
^^^^^^^^^

There are 3 centering options: manual, beam center, and geometric center.

Manual
######

Manual centering means that you as the user specify the position of the image that the profiles are taken from.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Centering

    fa = FieldProfileAnalysis(...)
    fa.analyze(..., centering=Centering.MANUAL)  # default is the middle of the image.

    # or specify a custom location
    fa.analyze(..., centering=Centering.MANUAL, position=(0.3, 0.8))
    # take profile at 30% height and 80% width

Beam center
###########

This is the default centering option. It first looks for the field to find the approximate center along each axis.
Then it extracts the profiles at these coordinates. This is helpful if you always want to be near the center of the field, even for offset fields or wedges.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Centering

    fa = FieldProfileAnalysis(...)
    fa.analyze(...)  # nothing special needed as it's the default

    # Specifying a position here will be ignored
    fa.analyze(..., centering=Centering.BEAM_CENTER, position=(0.1, 0.4))
    # this is allowed but will result in the same result as above

Geometric center
################

The geometric center will always find the middle pixel and extract the profiles from there.
This is helpful if you always want to be at the center of the image.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Centering

    fa = FieldProfileAnalysis(...)
    fa.analyze(..., centering=Centering.GEOMETRIC_CENTER)

.. _edge-type:

Edge detection
^^^^^^^^^^^^^^

Edge detection is important for determining the field width and beam center (which is often used for symmetry).
There are 3 detection strategies: FWHM, inflection via derivative, and inflection via the Hill/sigmoid/4PNLR function.

FWHM
####

The full-width half-max strategy is traditional and works for flat beams. It can give poor values for FFF beams.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Edge

    fa = FieldProfileAnalysis(...)
    fa.analyze(..., edge_type=Edge.FWHM)


Inflection (derivative)
#######################

The inflection point via the derivative is useful for both flat and FFF beams, and is thus the default for ``FieldProfileAnalysis``.
The method will find the positions of the max and min derivative of the values. Using a 0-crossing of the 2nd derivative
can be tripped up by noise so it is not used.

.. note::

    This method is recommended for high spatial resolution images such as the EPID, where the derivative has several points to use at the beam edge.
    It is not recommended for 2D device arrays.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Edge

    fa = FieldProfileAnalysis(...)  # nothing special needed as it's the default

    # you may also specify the edge smoothing value. This is a gaussian filter applied to the derivative just for the purposes of finding the min/max derivative.
    # This is to ensure the derivative is not caught by some noise. It is usually not necessary to change this.
    fa = FieldProfileAnalysis(..., edge_smoothing_ratio=0.005)


Inflection (Hill)
#################

The inflection point via the Hill function is useful for both flat and FFF beams
The fitting of the function is best for low-resolution data.
The Hill function, the sigmoid function, and 4-point non-linear regression belong to a family of logistic equations to fit a dual-curved value.
Since these fit a function to the data the resolution problems are eliminated. Some examples can be
seen `here <https://en.wikipedia.org/wiki/Sigmoid_function#Examples>`__. The generalized logistic function has helpful visuals as well
`here <https://en.wikipedia.org/wiki/Generalised_logistic_function>`__.

The function used here is:

:math:`f(x) = A + \frac{B - A}{1 + \frac{C}{x}^D}`

where :math:`A` is the low asymptote value (~0 on the left edge of a field),
:math:`B` is the high asymptote value (~1 for a normalized beam on the left edge),
:math:`C` is the inflection point of the sigmoid curve,
and :math:`D` is the slope of the sigmoid.

The function is fitted to the edge data of the field on each side to return the function. From there, the inflection point, penumbra, and slope can be found.

.. note::

    This method is recommended for low spatial resolution images such as 2D device arrays, where there is very little data at the beam edges.
    While it can be used for EPID images as well, the fit can have small errors as compared to the direct data.
    The fit, however, is much better than a linear or even spline interpolation at low resolutions.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Edge

    fa = FieldProfileAnalysis(..., edge_type=Edge.INFLECTION_HILL)

    # you may also specify the Hill window. This is the size of the window (as a ratio) to use to fit the field edge to the Hill function.
    fa = FieldProfileAnalysis(..., edge_type=Edge.INFLECTION_HILL, hill_window_ratio=0.05)
    # i.e. use a 5% field width about the edges to fit the Hill function.


.. note::

    When using this method, the fitted Hill function will also be plotted on the image. Further, the exact field edge
    marker (green x) may not align with the Hill function fit. This is just a rounding issue due to the plotting mechanism.
    The field edge is really using the Hill fit under the hood.

Normalization
-------------

There are 4 options for interpolation: ``NONE``, ``GEOMETRIC_CENTER``, ``BEAM_CENTER``, and ``MAX``. These should be
self-explanatory, especially in light of the centering explanations.

.. code-block:: python

    from pylinac import FieldProfileAnalysis, Normalization

    fa = FieldProfileAnalysis(...)
    fa.analyze(..., normalization_method=Normalization.BEAM_CENTER)


.. _choosing-field-metrics:

Choosing Metrics
================

The metrics the user can use in the field/profile analysis that come out of the box
can be viewed in this list: :ref:`profile_builtin_plugins`. These metrics can then
be passed to the ``metrics`` parameter. The user can mix and match as desired.
It's also possible to write custom plugins; see :ref:`using-custom-field-analysis-plugins`.

.. code-block:: python

    from pylinac import FieldProfileAnalysis
    from pylinac.metrics.profile import PenumbraRightMetric, SlopeMetric, SymmetryAreaMetric

    p = Path(r"C:\path\to\open_field.dcm")
    field_analyzer = FieldProfileAnalysis(p)
    field_analyzer.analyze(
        metrics=[PenumbraRightMetric(), SlopeMetric(), SymmetryAreaMetric()]
    )
    field_analyzer.plot_analyzed_images()

Customizing Metrics
===================

To change the name or parameters of a metric, pass the relevant parameters to the metric's constructor.
E.g. to use 90/10 for the penumbra distance instead of 80/20:

.. code-block:: python

  from pylinac import FieldProfileAnalysis
  from pylinac.metrics.profile import PenumbraRightMetric

  fa = FieldProfileAnalysis(...)
  fa.analyze(metrics=[PenumbraRightMetric(lower=10, upper=90, color="r", ls="--")])
  # will use 90/10 instead of 80/20 and plot the line in red and dashed
  fa.plot_analyzed_images()

Each metric has its own parameters that can be customized. See the metric's documentation for more information.

Interpreting results
====================

Results from the analysis can be printed as a simple string using :meth:`~pylinac.field_profile_analysis.FieldProfileAnalysis.results`.
The ideal method is to access the results using the :meth:`~pylinac.field_profile_analysis.FieldProfileAnalysis.results_data` method, which
can return a JSON string, a dictionary, or a :class:`~pylinac.field_profile_analysis.FieldResult` that can be used to access the results.

See also :ref:`exporting-results`.


The outcome from analyzing the phantom available in RadMachine or from
``results_data`` is:

* ``normalization``: The normalization method used.
* ``centering``: The centering method used.
* ``edge_type``: The edge detection method used.
* ``center``: The center ROI of the field with the following items:

  * ``mean``: The mean pixel value of the center ROI.
  * ``stdev``: The standard deviation of the pixel values in the center ROI.
  * ``min``: The minimum pixel value in the center ROI.
  * ``max``: The maximum pixel value in the center ROI.

* ``x_metrics``: The metrics computed from the X/crossplane profile.
  The items included depend on the metrics given to the analysis. See: :ref:`profile_builtin_plugins`. The
  key will be the name of the metric and the value will be what's returned from
  the ``calculate`` method of the metric. In addition, the following items
  are always given:

  * ``Field Width (mm)``: The width of the field in mm.
  * ``values``: A list of pixel values of the profile that the metrics were calculated from.

* ``y_metrics``: The metrics computed from the Y/inplane profile.
  The items and logic are the same as for the ``x_metrics``.

.. _using-custom-field-analysis-plugins:

Using Custom Plugins
====================

To analyze custom metrics on the profiles, the user can write custom plugins.
First, :ref:`create the plugin itself <writing-profile-plugins>`. Then,
pass it to the ``metrics`` parameter per :ref:`choosing-field-metrics`.

.. code-block:: python

    from pylinac import FieldProfileAnalysis
    from pylinac.metrics.profile import PenumbraLeftMetric
    from my_custom_plugins import MyCustomMetric

    p = Path(r"C:\path\to\open_field.dcm")
    field_analyzer = FieldProfileAnalysis(p)
    field_analyzer.analyze(metrics=[MyCustomMetric(), PenumbraLeftMetric()])
    field_analyzer.plot_analyzed_images()

.. _comparison-to-field-analysis-v1:

Comparison to field analysis v1
===============================

In comparison to the original :ref:`field_analysis_module` module,
this updated module:

* Almost every metric will result in the same value. E.g. flatness, symmetry, etc. The only
  times that metrics showed differences was if interpolation was "None". In such a case,
  the v1 penumbra may be slightly larger. This is due to the field edge index "snapping" to the nearest pixel.
  In v2, the edge is always interpolated and never snaps.
* Is much easier and more logical to write new plugins for. The plugins are modular and are not tightly coupled
  to this module.
* Is loosely coupled to the metrics themselves. Since the metrics themselves are
  self-contained, can be customized directly, and don't rely on any specific data massaging or formatting, they can be used in other modules.
* Does not do any device-specific analysis (e.g. IC Profiler). The original module and intent
  was to have a generic device-to-profile framework. This ended up being a considerable amount
  of work and ultimately we felt it fell out of the scope of pylinac, which is to focus on image
  analysis. If you desire device-specific profile analysis, the conversion from the original file/format
  a 1D array of values is your responsibility.
* There is no explicit interpolation method. Values are always interpolated linearly under the hood.
* v1 uses language like "top, left, bottom, right". In v2 everything is left and right, including the
  y-profile (vertical). Rotating the y-profile 90 degrees counterclockwise will
  give the same understanding of language. I.e. the "top" is the y-left value.

  .. image:: images/field_analysis_v2.png
      :width: 750
      :align: center

* While not strictly a user benefit, the modularity of the metrics and the analysis itself
  makes the codebase much cleaner and easier to maintain (v1 is ~1500 LOC vs v2 at ~400).


API Documentation
=================


.. autoclass:: pylinac.field_profile_analysis.FieldProfileAnalysis
    :members:
    :inherited-members:

.. autopydantic_model:: pylinac.field_profile_analysis.FieldProfileResult
