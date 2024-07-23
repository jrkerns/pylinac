.. _v4-changes:

==========
V4 Changes
==========

This page will track changes that will be made in version 4.x of pylinac.

.. important::

    As of right now, there is no timeline for the release of version 4.x.

.. warning::

    The changes listed here are subject to change. Some changes may be incorporated
    into 3.x if a deprecation policy and migration path that is not overbearing to
    users is determined.


The point is to give users a heads up on what to expect in the next major release so that
code can be changed or adjusted with plenty of preparation time.

Philosophy
----------

* Generally speaking and with some exceptions, pylinac will no longer be making pass/fail evaluations.
  Now that pylinac is a part of the RadMachine product, the ideal workflow is for a user to make
  pass/warn/fail evaluations within the RadMachine ecosystem. Obviously, many people that use
  pylinac don't use RadMachine, but this is the philosophical line that we're drawing: pylinac
  is for analysis and not evaluation. Subtly though, it is still likely that pylinac will
  allow *inputting* of tolerances, but the output will be more raw data oriented
  (e.g. the difference between nominal and measured) and not pass/fail. We do understand that
  seeing red/green in analyzed images is helpful, which is why we will continue to support
  inputting tolerances. To state this another way, results won't have things like ``passed``
  but may have things like ``deviation``.
* The ``results_data`` method is the primary method for getting analyzed data out. This
  is already true but we will be pressing into this even more. This may mean new sections
  that includes warnings or deprecations for certain attributes or methods. If you are
  looking up attributes directly from the class to retrieve analyzed data you should file
  a ticket so it can be added to the ``results_data`` method.


Core
----

* Several utility functions and methods will be removed to trim some "fat" from the package.

  * ``load_url``
  * ``.from_url(...)`` class methods. I've not seen a single use case for loading from a URL directly.

* The plotting of analyzed images will migrate from ``matplotlib`` to ``plotly``.
  This is mainly related to being able to zoom in and otherwise examine an analyzed image
  more closely. Once plotly is supported, there is no reason to keep both plotting methods
  around. To reduce maintenance burden, the matplotlib plotting will be removed.



Open-Field Analysis
-------------------

* The ``FieldAnalysis`` class and :ref:`module <field_analysis_module>` will be deprecated.
  Users should migrate to the :ref:`new module <field-profile-analysis>`. This new module
  uses a composable approach and is far more flexible for adding and editing metrics.
* The ``SNCProfiler`` class will be removed. The original aim was to have a module for
  parsing different file types from various vendors. Pylinac's scope has grown and evolved
  and creating a module for vendor file parsing is no longer a priority. If it
  does become a priority, it will be a separate package.

Winston-Lutz
------------

* The ``image_details`` attribute will be removed from the ``results_data`` call.
  The ``keyed_image_details`` attribute already exists and is the same data but
  is a dictionary keyed with the gantry, coll, couch angles. ``image_details``
  is simply a list and is not ordered in any meaningful way.

Cheese
------

* The ``roi_<n>`` keys from the ``TomoCheese`` class will be removed. Instead, the
  ``rois`` attribute already exists and is a dictionary of ROI objects. It is the
  same data but is more Pythonic and extensible. Other and future "cheese" phantoms
  will and only have the ``rois`` attribute and no ``roi_<n>`` keys.

Picket Fence
------------

* The ``mlc_positions_by_leaf`` attribute of the picket fence ``results_data`` will be
  **changed** from being from the top/left of the image to being relative to the CAX.
