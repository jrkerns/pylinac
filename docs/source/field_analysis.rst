.. _field_analysis_module:

==============
Field Analysis
==============

.. warning::

    This module will be deprecated in favor of the newer method of :ref:`field-profile-analysis` at some point in the future.
    This module will either be moved to a ``deprecated`` module in v4 or simply removed.

Overview
--------

The field analysis module (``pylinac.field_analysis``) allows a physicist to analyze metrics from an EPID
to measure penumbra, field width, etc. Additionally, protocols can be used which can calculate flatness & symmetry.
The module is very flexible, allowing users to choose different types of interpolation, normalization, centering, etc.
Users can also creat custom protocols to perform other types of field analysis within the main pylinac flow.

The module implements traditional analysis like FWHM as well as new methods as outlined in the pre-publication of the
`NCS-33 report <https://radiationdosimetry.org/files/Prepublication_-_NCS_Report_33_Beam_parameters_V2020-07-29.pdf>`__
which include edge fitting for FFF fields as well as a "top" calculation for the center position of FFF beams.

.. note::
    This is not a purely faithful implementation of NCS-33.
    There are a few differences w/r/t how stringent field widths are applied.
    E.g. the "top" calculation in NCS-33 is over the central 5cm.
    Pylinac simply uses a field width ratio which may or may not be 5cm.

The module's main class is :class:`~pylinac.field_analysis.FieldAnalysis` which is used for EPID images.

Running the Demo
----------------

To run the demo, import the main class and run the demo method:

.. plot::

    from pylinac import FieldAnalysis

    FieldAnalysis.run_demo()

Which will also result in the following output::

    Field Analysis Results
    ----------------------
    File: E:\OneDrive - F...demo_files\flatsym_demo.dcm
    Protocol: VARIAN
    Centering method: Beam center
    Normalization method: Beam center
    Interpolation: Linear
    Edge detection method: Inflection Derivative

    Penumbra width (20/80):
    Left: 2.7mm
    Right: 3.0mm
    Top: 3.9mm
    Bottom: 2.8mm

    Field Size:
    Horizontal: 140.9mm
    Vertical: 200.3mm

    CAX to edge distances:
    CAX -> Top edge: 99.8mm
    CAX -> Bottom edge: 100.5mm
    CAX -> Left edge: 60.4mm
    CAX -> Right edge: 80.5mm

    Top slope: -0.006%/mm
    Bottom slope: 0.044%/mm
    Left slope: 0.013%/mm
    Right slope: 0.014%/mm

    Protocol data:
    --------------
    Vertical symmetry: -2.631%
    Horizontal symmetry: -3.006%

    Vertical flatness: 1.700%
    Horizontal flatness: 1.857%

Typical Use
-----------

In most instances, a physicist is interested in quickly calculating the flatness, symmetry, or both of the
image in question. The ``field_analysis`` module allows you to do this easily and quickly.

To get started, import the :class:`~pylinac.field_analysis.FieldAnalysis` class:

.. code-block:: python

    from pylinac import FieldAnalysis

Loading images is easy and just like any other module:

.. code-block:: python

    # from a file
    my_file = r"C:/my/QA/folder/img.dcm"
    my_img = FieldAnalysis(path=my_file)

Alternatively, you can load data from a 2D device array:

.. code-block:: python

    from pylinac import DeviceFieldAnalysis

    # Profiler file
    my_file = r"C:/my/profiler/data.prm"
    my_img = DeviceFieldAnalysis(path=my_file)

If you don't have an image you can load the demo image:

.. code-block:: python

    my_img = FieldAnalysis.from_demo_image()

You can then calculate several field metrics with the :meth:`~pylinac.field_analysis.FieldAnalysis.analyze` method:

.. code-block:: python

    my_img.analyze()

After analysis, the results can be printed, plotted, or saved to a PDF:

.. code-block:: python

    print(my_img.results())  # print results as a string
    my_img.plot_analyzed_image()  # matplotlib image
    my_img.publish_pdf(filename="flatsym.pdf")  # create PDF and save to file
    my_img.results_data()  # dict of results


Analyze Options
---------------

The analysis algorithm allows the user to change numerous parameters about the analysis including automatic/manual
centering, profile extraction width, field width ratio, penumbra values, interpolation, normalization, and edge detection.
See :meth:`pylinac.field_analysis.FieldAnalysis.analyze` for details on each parameter

The protocol can also be specified here;
this is where both default and custom algorithms like flatness and symmetry can be used.
See :ref:`analysis_definitions` for the common flatness/symmetry algorithms provided out of the box.
For custom protocols, see :ref:`custom_protocols`.

.. code-block:: python

    from pylinac import Protocol, Centering, Edge, Normalization, Interpolation

    my_img.analyze(
        protocol=Protocol.ELEKTA,
        centering=Centering.BEAM_CENTER,
        in_field_ratio=0.8,
        is_FFF=True,
        interpolation=Interpolation.SPLINE,
        interpolation_resolution_mm=0.2,
        edge_detection_method=Edge.INFLECTION_HILL,
    )


.. _centering:

Centering
^^^^^^^^^

There are 3 centering options: manual, beam center, and geometric center.

Manual
######

Manual centering means that you as the user specify the position of the image that the profiles are taken from.

.. code-block:: python

    from pylinac import FieldAnalysis, Centering

    fa = FieldAnalysis(...)
    fa.analyze(..., centering=Centering.MANUAL)  # default is the middle of the image

    # or specify a custom location
    fa.analyze(..., centering=Centering.MANUAL, vert_position=0.3, horiz_position=0.8)
    # take profile at 30% width (i.e. left side) and 80% height

Beam center
###########

This is the default for EPID images/``FieldAnalysis``. It first looks for the field to find the approximate center along each axis.
Then it extracts the profiles and continues. This is helpful if you always want to be at the center of the field, even for offset fields or wedges.

.. code-block:: python

    from pylinac import FieldAnalysis, Centering

    fa = FieldAnalysis(...)
    fa.analyze(...)  # nothing special needed as it's the default

    # You CANNOT specify a position. These values will be ignored
    fa.analyze(..., centering=Centering.BEAM_CENTER, vert_position=0.3, horiz_position=0.8)
    # this is allowed but will result in the same result as above

Geometric center
################

This is the default for 2D device arrays/``DeviceFieldAnalysis``. It will always find the middle pixel and extract the profiles from there.
This is helpful if you always want to be at the center of the image.

.. code-block:: python

    from pylinac import FieldAnalysis, Centering

    fa = FieldAnalysis(...)
    fa.analyze(...)  # nothing special needed as it's the default

    # You CANNOT specify a position. These values will be ignored
    fa.analyze(
        ..., centering=Centering.GEOMETRIC_CENTER, vert_position=0.3, horiz_position=0.8
    )
    # this is allowed but will result in the same result as above

.. _edge:

Edge detection
^^^^^^^^^^^^^^

Edge detection is important for determining the field width and beam center (which is often used for symmetry).
There are 3 detection strategies: FWHM, inflection via derivative, and inflection via the Hill/sigmoid/4PNLR function.

FWHM
####

The full-width half-max strategy is traditional and works for flat beams. It can give poor values for FFF beams.

.. code-block:: python

    from pylinac import FieldAnalysis, Edge

    fa = FieldAnalysis(...)
    fa.analyze(..., edge_detection_method=Edge.FWHM)


Inflection (derivative)
#######################

The inflection point via the derivative is useful for both flat and FFF beams, and is thus the default for ``FieldAnalysis``.
The method will find the positions of the max and min derivative of the values. Using a 0-crossing of the 2nd derivative
can be tripped up by noise so it is not used.

.. note::

    This method is recommended for high spatial resolution images such as the EPID, where the derivative has several points to use at the beam edge.
    It is not recommended for 2D device arrays.

.. code-block:: python

    from pylinac import FieldAnalysis, Edge

    fa = FieldAnalysis(...)  # nothing special needed as it's the default

    # you may also specify the edge smoothing value. This is a gaussian filter applied to the derivative just for the purposes of finding the min/max derivative.
    # This is to ensure the derivative is not caught by some noise. It is usually not necessary to change this.
    fa = FieldAnalysis(..., edge_smoothing_ratio=0.005)


Inflection (Hill)
#################

.. warning::

    The Hill method is prone to errors if the penumbra is not expressed well on both edges. E.g.
    if using a device array with a full, open field where the penumbra sharply falls to 0 but does not
    have a long "tail" of low values, the Hill function will not fit well as a tail is expected on both
    edges.

The inflection point via the Hill function is useful for both flat and FFF beams
The fitting of the function is best for low-resolution data, and is thus the default for ``DeviceFieldAnalysis``.
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
    The fit, however, can be better than a linear or even spline interpolation at low resolutions. However,
    this should be combined with a resolution method like linear interpolation. While Hill can sometimes
    find a better fit for low-resolution data, the fit itself should have many data points to work with. I.e.
    the call should look similar to:

    .. code-block:: python

      fa = DeviceFieldAnalysis(...)
      fa.analyze(
          ..., interpolation=Interpolation.LINEAR, edge_detection_method=Edge.INFLECTION_HILL
      )


.. code-block:: python

    from pylinac import FieldAnalysis, Edge

    fa = FieldAnalysis(..., edge_detection_method=Edge.INFLECTION_HILL)

    # you may also specify the hill window. This is the size of the window (as a ratio) to use to fit the field edge to the Hill function.
    fa = FieldAnalysis(
        ..., edge_detection_method=Edge.INFLECTION_HILL, hill_window_ratio=0.05
    )
    # i.e. use a 5% field width about the edges to fit the Hill function.


.. note::

    When using this method, the fitted Hill function will also be plotted on the image. Further, the exact field edge
    marker (green x) may not align with the Hill function fit. This is just a rounding issue due to the plotting mechanism.
    The field edge is really using the Hill fit under the hood.

.. _normalization:

Normalization
-------------

There are 4 options for interpolation: ``None``, ``GEOMETRIC_CENTER``, ``BEAM_CENTER``, and ``MAX``. These should be
self-explanatory, especially in light of the centering explanations.

.. code-block:: python

    from pylinac import FieldAnalysis, Normalization

    fa = FieldAnalysis(...)
    fa.analyze(..., normalization_method=Normalization.BEAM_CENTER)

.. _interpolation:

Interpolation
-------------

There are 3 options for interpolation: ``NONE``, ``LINEAR``, and ``SPLINE``.

None
^^^^

A method of ``NONE`` will obviously apply no interpolation. Other interpolation parameters (see below) are ignored.
This is the default method for ``DeviceFieldAnalysis``

.. note::

    When plotting the data, if interpolation is None and the data is from a device, the data will be plotted as individual markers (``+``).
    If interpolation is applied to device data or it is a DICOM/EPID image, the data is plotted as a line.

Linear
^^^^^^

This will apply a linear interpolation to the original data. Along with this, the parameter interpolation_resolution_mm determine the
amount of interpolation. E.g. a value of 0.1 will resample the data to get data points 0.1mm apart.
This is the default method for ``FieldAnalysis``.

.. code-block:: python

    from pylinac import FieldAnalysis, Interpolation

    fa = FieldAnalysis(...)
    fa.analyze(..., interpolation=Interpolation.LINEAR, interpolation_resolution_mm=0.1)

Spline
^^^^^^

This will apply a cubic spline interpolation to the original data. Along with this, the parameter interpolation_resolution_mm determine the
amount of interpolation. E.g. a value of 0.1 will resample the data to get data points 0.1mm apart.

.. code-block:: python

    from pylinac import FieldAnalysis, Interpolation

    fa = FieldAnalysis(...)
    fa.analyze(..., interpolation=Interpolation.SPLINE, interpolation_resolution_mm=0.1)

.. _analysis_definitions:

Protocol Definitions
--------------------

There are multiple definitions for both flatness and symmetry. Your machine vendor uses certain equations,
or your clinic may use a specific definition. Pylinac has a number of built-in definitions which you can use.
Know also that you can create your own if you don't like/want to extend these :ref:`custom_protocols`.

None
^^^^

Technically, you are allowed a "None" protocol (``Protocol.NONE``), which just means nothing beyond the basic field analysis
is performed. If you just want the penumbra, distances to CAX, etc, without flatness/symmetry/custom algos then this is for you.

.. _varian_protocol:

Varian
^^^^^^

This is the default protocol if you don't specify one (``Protocol.VARIAN``). Two metrics are included, flatness & symmetry.

.. code-block:: python

    from pylinac import FieldAnalysis, Protocol

    fa = FieldAnalysis(...)
    fa.analyze(protocol=Protocol.VARIAN, ...)
    ...


Flatness
########

Flatness is defined by the variation (difference) across the field values within the field width.

:math:`flatness = 100 * |D_{max} - D_{min}| / (D_{max} + D_{min})`

If the field width is set to, e.g. 80%, then the flatness is calculated over all the values within that 80%.
Flatness is a scalar and always positive.

Symmetry
########

Symmetry is defined as the Point Difference:

:math:`symmetry = 100 * max(|L_{pt} - R_{pt}|)/ D_{CAX}`

where :math:`L_{pt}` and :math:`R_{pt}` are equidistant from the beam center.

Symmetry is calculated over the specified field width (e.g. 80%) as set in by ``analyze()``.
Symmetry can be positive or negative. A negative value means the right side is higher.
A positive value means the left side is higher.

.. _elekta_protocol:

Elekta
^^^^^^

This is specified by passing ``protocol=Protocol.ELEKTA`` to ``analyze``.

.. code-block:: python

    from pylinac import FieldAnalysis, Protocol

    fa = FieldAnalysis(...)
    fa.analyze(protocol=Protocol.ELEKTA, ...)
    ...

Flatness
########

Flatness is defined by the ratio of max/min across the field values within the field width.

:math:`flatness = 100 * D_{max}/D_{min}`

If the field width is set to, e.g. 80%, then the flatness is calculated over all the values within that 80%.
Flatness is a scalar and always positive.

Symmetry
########

Symmetry is defined as the Point Difference Quotient (aka IEC):

:math:`symmetry = 100 * max(|L_{pt}/R_{pt}|, |R_{pt}/L_{pt}|)`

where :math:`L_{pt}` and :math:`R_{pt}` are equidistant from the beam center.

Symmetry is calculated over the specified field width (e.g. 80%) as set in by ``analyze()``.
Symmetry can be positive or negative. A negative value means the right side is higher.
A positive value means the left side is higher.

.. _siemens_protocol:

Siemens
^^^^^^^

This is specified by passing ``protocol=Protocol.SIEMENS`` to ``analyze``.

.. code-block:: python

    from pylinac import FieldAnalysis, Protocol

    fa = FieldAnalysis(...)
    fa.analyze(protocol=Protocol.SIEMENS, ...)
    ...

Flatness
########

Flatness is defined by the variation (difference) across the field values within the field width.

:math:`flatness = 100 * |D_{max} - D_{min}| / (D_{max} + D_{min})`

If the field width is set to, e.g. 80%, then the flatness is calculated over all the values within that 80%.
Flatness is a scalar and always positive.

Symmetry
########

Symmetry is defined as the ratio of area on each side about the CAX:

:math:`symmetry = 100 * (A_{left} - A_{right}) / (A_{left} + A_{right})`

Symmetry is calculated over the specified field width (e.g. 80%) as set in by ``analyze()``.
Symmetry can be positive or negative. A negative value means the right side is higher.
A positive value means the left side is higher.


.. _custom_protocols:

Creating & Using Custom Protocols
---------------------------------

Protocols allow the user to perform specific image metric algorithms. This includes things like flatness & symmetry.
Depending on the protocol, different methods of determining the flatness/symmetry/whatever exist. Pylinac provides a handful
of protocols out of the box, but it is easy to add your own custom algorithms.

To create a custom protocol you must 1) create custom algorithm functions, 2) create a protocol class that 3) inherits from Enum and 4) defines a dictionary
with a ``calc``, ``unit``, and ``plot`` key/value pair. The ``plot`` key is optional; it allows you to plot something if you
also want to see your special algorithm (e.g. if it used a fitting function and you want to plot the fitted values).

* ``calc`` should be a function to calculate a specific, singular value such as flatness.
* ``unit`` should be a string that specifies the unit of ``calc``. If it is unitless leave it as an empty string (``''``)
* ``plot`` is **OPTIONAL** and is a function that can plot something to the profile views (e.g. a fitting function)

The ``calc`` and ``plot`` values should be functions with
a specific signature as shown in the example below:

.. code-block:: python

    import enum


    # create the custom algorithm functions
    # the ``calc`` function must have the following signature
    def my_special_flatness(
        profile: SingleProfile, in_field_ratio: float, **kwargs
    ) -> float:
        # do whatever. Must return a float. ``profile`` will be called twice, once for the vertical profile and horizontal profile.
        # the kwargs are passed to ``analyze`` and can be used here for your own purposes (e.g. fitting parameters)
        my_special_value = kwargs.get("funkilicious")
        flatness = profile(...)
        return flatness


    # custom plot function for the above flatness function
    # This is OPTIONAL
    # If you do implement this, it must have the following signature
    def my_special_flatness_plot(instance, profile: SingleProfile, axis: plt.Axes) -> None:
        # instance is the FieldAnalysis instance; i.e. it's basically ``self``.
        # do whatever; typically, you will do an axis.plot()
        axis.plot(...)


    # custom protocol MUST inherit from Enum
    class MySpecialProtocols(enum.Enum):
        # note you can specify several protocols if you wish
        PROTOCOL_1 = {
            # for each protocol, you can specify any number of metrics to calculate. E.g. 2 symmetry calculations
            "my flatness": {
                "calc": my_special_flatness,
                "unit": "%",
                "plot": my_special_flatness_plot,
            },
            "my symmetry": ...,
            "my other flatness metric": ...,
        }
        PROTOCOL_2 = ...


    # proceed as normal
    fa = FieldAnalysis(...)
    fa.analyze(protocol=MySpecialProtocols.PROTOCOL_1, ...)
    ...

Passing in custom parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may pass custom parameters to these custom algorithms via the analyze method as simple keyword arguments:

.. code-block:: python

    fa = FieldAnalysis(...)
    fa.analyze(..., my_special_variable=42)

The parameter will then be passed to the custom functions:

.. code-block:: python

    def my_special_flatness(
        profile: SingleProfile, in_field_ratio: float, **kwargs
    ) -> float:
        my_special_value = kwargs.get("my_special_variable")  # 42
        flatness = profile(...)
        return flatness

.. note:: The :class:`~pylinac.core.profile.SingleProfile` passed to the functions is very powerful and can calculate
          numerous helpful data for you such as the field edges, minimum/maximum value within the field, and much more.
          Read the SingleProfile documentation before creating a custom algorithm.

.. _fff_fields:

FFF fields
----------

The field analysis module can handle FFF beams, or more specifically, calculating extra metrics associated with FFF fields.
These metrics are largely from the NCS-33 pre-publication and include the "top" position, and the slopes of the field on each side.

These metrics are always calculated (even for flat beams), but will be shown in the ``results()`` output and also on the
plotted image of ``plot_analyzed_image()`` if the ``is_FFF`` flag is true.

The easiest way to demonstrate this is through the Device demo, which is an FFF field:

.. plot::

    from pylinac import DeviceFieldAnalysis, Protocol

    fa = DeviceFieldAnalysis.from_demo_image()
    fa.analyze(protocol=Protocol.VARIAN, is_FFF=True)
    fa.plot_analyzed_image()

.. _field_analysis_top_metric:

"Top" metric
^^^^^^^^^^^^

The "top" metric is the fitted position of the peak of a FFF beam. It uses the central region of the field as
specified by the ``slope_exclusion_ratio``. E.g. if the value is 0.3 it will use the central 30% field width.

The central region is fitted to a 2nd order polynomial and then the max of the polynomial is found. That value is the
"top" position. This helps to account for noise in the profile.

When printing results for an FFF beam there will be a section like so::

    'Top' vertical distance from CAX: 1.3mm
    'Top' horizontal distance from CAX: 0.6mm
    'Top' vertical distance from beam center: 1.7mm
    'Top' horizontal distance from beam center: 0.3mm

.. _field_analysis_slope_metric:

Field slope
^^^^^^^^^^^

For FFF beams, an additional metric is calculated: the slope of each side of the field. Since traditional flatness
algorithms aren't tuned for FFF beams they can be noisy or non-sensible. By measuring the slope of each side of the field
the flatness can be measured more accurately (as a slope) for trending and consistency purposes.

The slope is calculated in the regions between the field width edges and the slope exclusion ratio. E.g. a field width ratio
of 0.8 and a slope exclusion ratio of 0.4 will mean that the regions between +/-0.4 (0.8/2) from the CAX to +/-0.2 (0.4/2)
will be used to fit linear regressions.

When printing results for an FFF beam there will be a section like so::

    Top slope: 0.292%/mm
    Bottom slope: -0.291%/mm
    Left slope: 0.295%/mm
    Right slope: -0.296%/mm

Accessing data
--------------

.. versionchanged:: 3.0

Using the module in your own scripts? While the analysis results can be printed out,
if you intend on using them elsewhere (e.g. in an API), they can be accessed the easiest by using the :meth:`~pylinac.field_analysis.FieldAnalysis.results_data` method
which returns a :class:`~pylinac.field_analysis.FieldResult` instance.

.. note::
    While the pylinac tooling may change under the hood, this object should remain largely the same and/or expand.
    Thus, using this is more stable than accessing attrs directly.

You can access most data you get from ``results()``:

.. code-block:: python

    fa = FieldAnalysis(...)
    fa.analyze(...)
    data = fa.results_data()

    data.top_penumbra_mm
    data.beam_center_to_left_mm

You may also access protocol data in the ``protocol_results`` dictionary. These results must be in a dictionary because
the protocol names and fields are dynamic and not known a priori.

.. code-block:: python

    data.protocol_results["flatness_vertical"]
    data.protocol_results["symmetry_horizontal"]

The keys of this dict are defined by the protocol names. Using the example from the :ref:`custom_protocols` section,
we would access that custom protocol data as:

.. code-block:: python

    data.protocol_results["my flatness_vertical"]
    data.protocol_results["my flatness_horizontal"]

because the protocol name was ``my flatness``.

Analysis Parameters
-------------------

.. tab-set::
   :sync-group: usage

   .. tab-item:: pylinac
      :sync: pylinac

      See :meth:`~pylinac.field_analysis.FieldAnalysis.analyze` for details.

   .. tab-item:: RadMachine
      :sync: radmachine

      * **Protocol**: The protocol to use for the analysis. See :ref:`analysis_definitions`.
      * **Centering method**: The method used to determine the center of the field. See :ref:`centering`.
      * **is FFF beam**: Whether the beam is a Flattening Filter Free beam. If so, the "Top" metric is calculated. See :ref:`fff_fields`.
      * **Normalized vertical position**: The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
        in the middle of the image. 0.0 is at the left edge of the image and 1.0 is at the right edge of the image.

        .. note::

            This is only used when centering is set to "MANUAL".

      * **Normalized horizontal position**: The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
        in the middle of the image. 0.0 is at the top edge of the image and 1.0 is at the bottom edge of the image.

        .. note::

            This is only used when centering is set to "MANUAL".

      * **Normalized vertical width**: The width of the vertical profile to take. 0 is a single pixel, 1 is the entire image. The resulting
        profile will be the mean of the pixels within the specified width.
      * **Normalized horizontal width**: The width of the horizontal profile to take. 0 is a single pixel, 1 is the entire image. The resulting
        profile will be the mean of the pixels within the specified width.
      * **In-field ratio**: The ratio of the field width to use for protocol analyses. E.g. 0.8 means the central 80% of the field is used for flatness and symmetry.
      * **Slope exclusion ratio**: This is the ratio of the field to use to 1) calculate the "top" metric for an FFF beam (:ref:`field_analysis_top_metric`) and 2) exclude from the
        "slope" calculation of each side of the field (:ref:`field_analysis_slope_metric`).
      * **Lower penumbra**: The lower penumbra value in % to use for the field edge.
      * **Upper penumbra**: The upper penumbra value in % to use for the field edge.
      * **Interpolation**: Interpolation technique to use. See :ref:`Interpolation`. Interpolation is usually only helpful
        for low-resolution data.
      * **Interpolation resolution**: The resolution that the interpolation will scale to.
        E.g. if the native dpmm is 2 and the resolution is set to 0.1mm the data will be interpolated to have a new dpmm of 10 (1/0.1).
      * **Ground the profile**: If True, the lowest value of the profile will be set to 0.
      * **Normalization method**: The method used to normalize the field. See :ref:`normalization`.
      * **Edge detection method**: The method used to detect the field edges. See :ref:`edge`.
      * **Edge smoothing ratio**:

        .. note:: Only applies to Inflection Derivative and Inflection Hill Edge detection methods.

        The ratio of the length of the values to use as the sigma for a Gaussian filter applied before searching for
        the inflection. E.g. 0.005 with a profile of 1000 points will result in a sigma of 5.
        This helps make the inflection point detection more robust to noise. Increase for noisy data.
      * **Hill window ratio**:

        .. note:: Only applies to Inflection Hill Edge detection methods.

        The ratio of the field size to use as the window to fit the Hill function. E.g. 0.2 will using a window
        centered about each edge with a width of 20% the size of the field width.


Algorithm
---------

There is little of a true "algorithm" in ``field_analysis`` other than analyzing profiles. Thus, this section is more terminology and
notekeeping.

**Allowances**

* The image can be any size.
* The image can be EPID (actually just DICOM) or a 2D device array file.
* The image can be either inversion (Radiation is dark or light).
* The image can be off-center.

**Restrictions**

* The module is only meant for photon analysis at the moment (there are sometimes different equations for electrons for the same
  definition name).
* Analysis is limited to normal/parallel directions. Thus if the image is rotated there is no way to account for it other than
  rotating the image before analysis.

**Analysis**

* **Extract profiles** - With the positions given, profiles are extracted and analyzed according to the method specified (see
  :ref:`analysis_definitions`). For symmetry calculations that operate around the CAX, the CAX must first be determined, which is
  the center of the FWHM of the profile.

Interpreting Results
--------------------

The field analysis result is slightly different if analyzing a device array vs
a 2D image (EPID/DICOM). In both cases the following is given:

* ``protocol``: The protocol used for the analysis. See: :ref:`analysis_definitions`.
* ``protocol_results``: This is a dictionary that contains the results of the protocol calculations. The keys are the protocol names.
  Generally, they are ``symmetry_horizontal``, ``symmetry_vertical``, ``flatness_horizontal``, and ``flatness_vertical``.
  The values themselves are the calculated values for that specific protocol equation.
* ``centering_method``: The method used to determine the center of the field. See :ref:`centering`.
* ``normalization_method``: The method used to normalize the field. See :ref:`normalization`.
* ``interpolation_method``: The method used to interpolate the field. See :ref:`interpolation`.
* ``edge_detection_method``: The method used to detect the field edges. See :ref:`edge`.
* ``top_penumbra_mm``: The penumbra width at the top of the field.
* ``bottom_penumbra_mm``: The penumbra width at the bottom of the field.
* ``left_penumbra_mm``: The penumbra width at the left of the field.
* ``right_penumbra_mm``: The penumbra width at the right of the field.
* ``geometric_center_index_x_y``: The geometric center of the field in pixel coordinates.
* ``beam_center_index_x_y``: The beam center of the field in pixel coordinates.
* ``field_size_vertical_mm``: The vertical field size in mm.
* ``field_size_horizontal_mm``: The horizontal field size in mm.
* ``beam_center_to_top_mm``: The distance from the beam center to the top edge of the field.
* ``beam_center_to_bottom_mm``: The distance from the beam center to the bottom edge of the field.
* ``beam_center_to_left_mm``: The distance from the beam center to the left edge of the field.
* ``beam_center_to_right_mm``: The distance from the beam center to the right edge of the field.
* ``cax_to_top_mm``: The distance from the CAX to the top edge of the field.
* ``cax_to_bottom_mm``: The distance from the CAX to the bottom edge of the field.
* ``cax_to_left_mm``: The distance from the CAX to the left edge of the field.
* ``cax_to_right_mm``: The distance from the CAX to the right edge of the field.
* ``top_position_index_x_y``: The top position of the field in pixel coordinates.
* ``top_horizontal_distance_from_cax_mm``: The horizontal distance from the top position to the image center.

  .. note::
    This says CAX, and by this we mean the CAX of the machine, regardless of field size/position.

* ``top_vertical_distance_from_cax_mm``: The vertical distance from the top position to the image center.
* ``top_horizontal_distance_from_beam_center_mm``: The horizontal distance from the top position to the beam center.
* ``top_vertical_distance_from_beam_center_mm``: The vertical distance from the top position to the beam center.
* ``left_slope_percent_mm``: The slope of the left side of the in-field region in percent per mm.
* ``right_slope_percent_mm``: The slope of the right side of the in-field region in percent per mm.
* ``top_slope_percent_mm``: The slope of the top side of the in-field region in percent per mm.
* ``bottom_slope_percent_mm``: The slope of the bottom side of the in-field region in percent per mm.
* ``top_penumbra_percent_mm``: The penumbra width at the top of the field in percent per mm.
* ``bottom_penumbra_percent_mm``: The penumbra width at the bottom of the field in percent per mm.
* ``left_penumbra_percent_mm``: The penumbra width at the left of the field in percent per mm.
* ``right_penumbra_percent_mm``: The penumbra width at the right of the field in percent per mm.

If analyzing a 2D image (EPID/DICOM) the following is also given.
The ROI in this case is where the finite profile areas overlap. E.g.
if ``vert_width`` and ``horiz_width`` are 0.1, there is a central 10% of the field that is used for the following ROI:

* ``central_roi_mean``: The mean pixel value of the central region of interest.
* ``central_roi_std``: The standard deviation of the pixel values of the central region of interest.
* ``central_roi_max``: The maximum pixel value of the central region of interest.
* ``central_roi_min``: The minimum pixel value of the central region of interest.

API Documentation
-----------------

Main classes
^^^^^^^^^^^^

These are the classes a typical user may interface with.

.. autoclass:: pylinac.field_analysis.FieldAnalysis
    :members:

.. autopydantic_model:: pylinac.field_analysis.FieldResult

.. autopydantic_model:: pylinac.field_analysis.DeviceResult

.. autoclass:: pylinac.field_analysis.Device
    :members:

.. autoclass:: pylinac.field_analysis.Protocol
    :members:

.. autoclass:: pylinac.field_analysis.Centering
    :members:

.. autoclass:: pylinac.field_analysis.Interpolation
    :members:

.. autoclass:: pylinac.field_analysis.Edge
    :members:

Supporting Classes
^^^^^^^^^^^^^^^^^^

You generally won't have to interface with these unless you're doing advanced behavior.

.. autofunction:: pylinac.field_analysis.flatness_dose_difference

.. autofunction:: pylinac.field_analysis.flatness_dose_ratio

.. autofunction:: pylinac.field_analysis.symmetry_point_difference

.. autofunction:: pylinac.field_analysis.symmetry_area

.. autofunction:: pylinac.field_analysis.symmetry_pdq_iec
