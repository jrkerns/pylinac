
.. _scale:

Machine Scale
-------------

A machine scale or coordinate system is a specific reference framework
used in the context of machines or equipment.
It provides a standardized way to define positions, orientations, and measurements within the machine or equipment.

Other examples of coordinate systems are Cartesian, Polar, and Spherical.

Within the context of pylinac and machines, we're generally referring to IEC61217 or other systems defined by a
manufacturer.

It is sometimes required to convert coordinate systems. An example is Dan Low's Winston-Lutz paper, where the
coordinate system must be corrected for the equations to work correctly. See :ref:`here <wl-algorithm>`.

Converting scales
^^^^^^^^^^^^^^^^^

To convert the values from one machine scale to another, use the :func:`~pylinac.core.scale.convert` function.

.. code-block:: python

    from pylinac.core.scale import convert, MachineScale

    gantry = 0
    coll = 90
    couch = 45

    new_gantry, new_coll, new_couch = convert(
        input_scale=MachineScale.Varian,
        output_scale=MachineScale.IEC61217,
        gantry=gantry,
        collimator=coll,
        rotation=couch,
    )

.. autoclass:: pylinac.core.scale.MachineScale
    :members:
    :inherited-members:

.. autofunction:: pylinac.core.scale.convert
