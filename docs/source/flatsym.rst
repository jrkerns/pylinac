.. _flatsym_module:

======================================
Flatness/Symmetry module documentation
======================================

Overview
--------

The Flatness & Symmetry module (``pylinac.flatsym``) allows a physicist to check their linac's beam profile
against well-known flatness/symmetry calculation standards, reporting absolute values. Film or EPID images
can be loaded in and analyzed. There is one main class :class:`~pylinac.flatsym.FlatSym`. If you will be using
custom algorithms there are two module dictionaries you will need: ``FLATNESS_EQUATIONS`` and ``SYMMETRY_EQUATIONS``.

Running the Demo
----------------

To run the demo, import the main class and run the demo method:

.. plot::

    from pylinac import FlatSym

    FlatSym.run_demo()

Which will also result in the following output::

    Flatness & Symmetry
    File: C:\Users\jkern\...inac\pylinac\demo_files\flatsym_demo.dcm

    Flatness method: Varian
    Vertical flatness: 1.931%
    Horizontal flatness: 1.857%
    Symmetry method: Varian
    Vertical symmetry: 2.462%
    Horizontal symmetry: 2.989%

    Penumbra (80/20):
    Horizontal: 2.6mm
    Vertical: 2.8mm

    Field Size:
    Horizontal: 140.9mm
    Vertical: 200.2mm

    CAX to edge distances:
    CAX -> Upper edge: 99.9mm
    CAX -> Lower edge: 100.4mm
    CAX -> Left edge: 60.4mm
    CAX -> Right edge: 80.6mm

Typical Use
-----------

In most instances, a physicist is interested in quickly calculating the flatness, symmetry, or both of the
image in question. The ``flatsym`` module allows you to do this easily, using any of multiple definitions of flatness
or symmetry.

To get started, import the :class:`~pylinac.flatsym.FlatSym` class:

.. code-block:: python

    from pylinac import FlatSym

Loading images is easy and just like any other module:

.. code-block:: python

    # from a file
    my_file = r"C:/my/QA/folder/img.dcm"
    my_img = FlatSym(path=my_file)

If you don't have an image you can load the demo image:

.. code-block:: python

    my_img = FlatSym.from_demo_image()

You can then calculate the flatness and symmetry with the :meth:`~pylinac.flatsym.FlatSym.analyze` method:

.. code-block:: python

    my_img.analyze(flatness_method='varian', symmetry_method='varian', vert_position=0.5, horiz_position=0.5)

After analysis, the results can be printed, plotted, or saved to a PDF:

.. code-block:: python

    print(my_img.results())  # print results
    my_img.plot_analyzed_image()  # matplotlib image
    my_img.publish_pdf(filename="flatsym.pdf")  # create PDF and save to file

Raw Data
--------

The raw data values of analysis are also available within the public attributes of the class:

.. code-block:: python

    my_img = FlatSym.from_demo_image()
    my_img.analyze(flatness_method='varian', symmetry_method='varian')
    my_img.symmetry['horizontal']['value']  # the actual symmetry value
    my_img.flatness['vertical']['value']

Analysis Options
----------------

The flatness/symmetry algorithms can be specified as well as the position of the analysis within the image and the width of the profile.
See :ref:`analysis_definitions` for the common algorithms.

.. code-block:: python

    my_img.analyze(flatness_method='elekta', symmetry_method='point difference', vert_position=0.4, horiz_position=0.6,
    vert_width=0.05, horiz_width=0.05)

You can also create your own algorithms.

.. _analysis_definitions:

Analysis Definitions
--------------------

.. warning:: The following definitions are for **photons** only.

There are multiple definitions for both flatness and symmetry. Your machine vendor uses certain equations,
or your clinic may use a specific definition. Pylinac has a number of built-in definitions which you can use.

Symmetry:

================================ ================================= ========== ================================================================================================================================
Name                             Parameter values                  Vendors    Equation
================================ ================================= ========== ================================================================================================================================
Point Difference                 ``varian``, ``point difference``  Varian     :math:`100 * max(|L_{pt} - R_{pt}|)/ D_{CAX}` over 80%FW, where :math:`L_{pt}` and :math:`R_{pt}` are equidistant from the CAX.
Point Difference Quotient (IEC)  ``elekta``, ``pdq iec``           Elekta     :math:`100 * max(|L_{pt}/R_{pt}|, |R_{pt}/L_{pt}|)` over 80%FW if 10<FW<30cm [#elekta]_
================================ ================================= ========== ================================================================================================================================



* -- **Parameter value(s), Name, Vendors that use it -- Equation**
* -- ``varian``, ``point difference``, **Point Difference, Varian** -- :math:`100 * max(|L_{pt} - R_{pt}|)/ D_{CAX}` over 80%FW, where :math:`L_{pt}` and :math:`R_{pt}` are
  equidistant from the CAX.
* -- ``elekta``, ``pdq iec``, **Point Difference Quotient (IEC), Elekta** -- :math:`100 * max(|L_{pt}/R_{pt}|, |R_{pt}/L_{pt}|)` over 80%FW if 10<FW<30cm [#elekta]_.

Flatness:

* -- **Parameter value(s), Name, Vendors that use it -- Equation**
* -- ``varian``, ``vom80``, ``siemens``, **Variation over mean (80%), Varian** -- :math:`100 * |D_{max} - D_{min}| / (D_{max} + D_{min})` within 80%FW.
* -- ``elekta``, ``iec``, **Dmax/Dmin (IEC), Elekta** -- :math:`100 * D_{max}/D_{min}` within 80%FW for 10<FW<30cm [#elekta]_.

.. note:: Siemens and other definitions (e.g. Area, Area/2) will be added if the community `asks <https://github.com/jrkerns/pylinac/issues>`_ for it.

.. [#elekta]
    The region calculated over actually varies by the following: for 5<FW<10cm, FW - 2*1cm; for 10<FW<30cm, FW - 2*0.1*FW (i.e. 80%FW);
    for 30cm<FW, FW - 2*6cm. Pylinac currently only uses the 80%FW no matter the FW, but accounting for the FW will come in a future
    version.

Creating & Using Custom Algorithms
----------------------------------

Custom Algorithm Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

The flatness/symmetry algorithms can easily be extended. All algorithms are in a module dictionary, so the required
steps are to 1) create the custom analysis function and 2) add it to the equation dictionary.

The custom algorithms must be functions and follow the below structure. All names can be changed, but the structure of
one input and correct number and order of outputs is fixed:

.. code-block:: python

    def custom_flatness(profile: SingleProfile) -> Sequence[float]:
        ...
        return flatness_value, max_value, min_value, left_edge_index, right_edge_index

    def custom_symmetry(profile: SingleProfile) -> Tuple[float, Sequence[float], float, float]:
        ...
        return symmetry_value, symmetry_array, left_edge, right_edge

While most values are easily understood the ``symmetry_array`` should be a list of floats or numpy array of the symmetry value at each point,
assuming to start at the CAX and move outward.

.. note:: When creating custom algorithms, use the above for guidance. It 1) must be a function and
          2) must have one input argument, which will be a :class:`~pylinac.core.profile.SingleProfile`.
          Depending on whether it is a flatness or symmetry calculation the return values must follow the above
          structure. For further examples, see the source code for the built-in equations;
          e.g. :func:`~pylinac.flatsym.flatness_varian`.

.. note:: The :class:`~pylinac.core.profile.SingleProfile` given to the function is very powerful and can calculate
          numerous helpful data for you such as the field edges, minimum/maximum value within the field, and much more.
          Read the documentation before creating a custom algorithm.

.. note:: For flatness equations, if there is no ``max_value`` or ``min_value`` (e.g. area calculations) then just set
          those values to ``0``. Additionally, if you do not want or have a ``symmetry_array`` to plot, then pass ``0`` and
          it will not be plotted.

Using Custom Algorithms
^^^^^^^^^^^^^^^^^^^^^^^

To use the custom algorithm we must now add it to the equation dictionary using a custom name.
The name can be anything you like that is a valid dictionary key:

.. code-block:: python

    from pylinac import FlatSym
    from pylinac.flatsym import SYMMETRY_EQUATIONS, FLATNESS_EQUATIONS

    def custom_flatness(profile: SingleProfile):
        ...

    FLATNESS_EQUATIONS['my custom flatness'] = custom_flatness  # no parentheses; we don't want to call it, just assign it

    my_img = FlatSym.from_demo_image()
    my_img.analyze(flatness_method='my custom flatness', ...)  # use whatever key you defined in the equation dict
    ...

Algorithm
---------

There is little of a true "algorithm" in ``flatsym`` other analyzing profiles. Thus, this section is more terminology and
notekeeping.

**Allowances**

* The image can be any size.
* The image can be digitized film or EPID (most image formats and DICOM).
* The image can be either inversion (Radiation is dark or light).

**Restrictions**

* The module is only meant for photon analysis at the moment (there are sometimes different equations for electrons for the same
  definition name).
* Analysis is limited to normal/parallel directions. Thus if the image is rotated there is no way to account for it other than
  rotating the image before analysis.

**Analysis**

* *Extract profiles* - With the positions given, profiles are extracted and analyzed according to the method specified (see
  :ref:`analysis_definitions`). For symmetry calculations that operate around the CAX, the CAX must first be determined, which is
  the center of the FWHM of the profile.

API Documentation
-----------------

.. autoclass:: pylinac.flatsym.FlatSym
    :no-show-inheritance:

.. autofunction:: pylinac.flatsym.flatness_varian

.. autofunction:: pylinac.flatsym.flatness_elekta

.. autofunction:: pylinac.flatsym.symmetry_point_difference

.. autofunction:: pylinac.flatsym.symmetry_pdq_iec
