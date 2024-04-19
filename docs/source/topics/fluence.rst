
.. _plan_fluence:

=====================
Plotting Plan Fluence
=====================

We can plot fluence from an RT Plan using the ``plan_generator`` module.
This is a simplistic approach.

.. warning::

    This should not be used for evaluating patient treatments or fluences. It is for QA/research purposes only.

.. note::

    This will generate the fluence of the MLCs only. Collimator rotation is not accounted for.
    Further, jaws are shown as a rectangle outline, but the fluence calculation does not
    consider them. This is on purpose; accounting for jaw position in the fluence itself will not be supported.

First, we need to load a plan and then call ``plot_fluences``.

.. code-block:: python

    import pydicom
    from pylinac.plan_generator.fluence import plot_fluences

    plan = pydicom.dcmread("path/to/plan.dcm")
    plot_fluences(plan, width_mm=200, resolution_mm=0.5)

We can also generate the fluence as a 3D numpy array using the ``generate_fluences`` function.

.. code-block:: python

    from pylinac.plan_generator.fluence import generate_fluences

    fluences = generate_fluences(plan, width_mm=200, resolution_mm=0.5)

The fluence is of shape (num_beams, height, width).

API Documentation
-----------------

.. autofunction:: pylinac.plan_generator.fluence.plot_fluences

.. autofunction:: pylinac.plan_generator.fluence.generate_fluences
