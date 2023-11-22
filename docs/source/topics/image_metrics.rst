.. _image-metrics:

Custom Image Metrics
====================

.. versionadded:: 3.16

Pylinac images can now have arbitrary metrics calculated on them, similar to profiles.
This can be useful for calculating and finding values and regions of interest in images.
The system is quite flexible and allows for any number of metrics to be calculated on an image.
Furthermore, this allows for re-usability of metrics, as they can be applied to any image.


Use Cases
---------

* Calculate the mean pixel value of an area of an image.
* Finding an object in the image.
* Calculating the distance between two objects in an image.

Basic Usage
-----------

To calculate metrics on an image, simply pass the metric(s) to the ``compute`` method of the image:

.. code-block:: python

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  metric = img.compute(
      metrics=DiskLocator(
          expected_position=(100, 100),
          search_window=(30, 30),
          radius=10,
          radius_tolerance=2,
      )
  )
  print(metric)

You may compute multiple metrics by passing a list of metrics:

.. code-block:: python

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  metrics = img.compute(
      metrics=[
          # disk 1
          DiskLocator(
              expected_position=(100, 100),
              search_window=(30, 30),
              radius=10,
              radius_tolerance=2,
          ),
          # disk 2
          DiskLocator(
              expected_position=(200, 200),
              search_window=(30, 30),
              radius=10,
              radius_tolerance=2,
          ),
      ]
  )
  print(metrics)

Metrics might have something to plot on the image. If so, the ``plot`` method of the image will plot the metric(s) on the image:

.. code-block:: python

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  metrics = img.compute(
      metrics=[
          # disk 1
          DiskLocator(
              expected_position=(100, 100),
              search_window=(30, 30),
              radius=10,
              radius_tolerance=2,
          ),
          # disk 2
          DiskLocator(
              expected_position=(200, 200),
              search_window=(30, 30),
              radius=10,
              radius_tolerance=2,
          ),
      ]
  )
  img.plot()  # plots the image with the BB positions overlaid

Built-in Metrics
----------------

Out of the box, three metrics currently exist:
:class:`~pylinac.core.metrics.DiskLocator`, :class:`~pylinac.core.metrics.DiskRegion` and
:class:`~pylinac.core.metrics.GlobalDiskLocator`.

These metrics will find disks, usually BBs, in an image and then return the location or region properties.

Single Disk Locators
^^^^^^^^^^^^^^^^^^^^

.. note::

  The values provided below are in pixels. The following sections show how variants of how to use the metrics
  using physical units and relative to the center of the image.

Here's an example of using the ``DiskLocator``:

.. code-block:: python
  :caption: Search for a disk 100 pixels right and 100 pixels down from the top left of the image

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  img.compute(
      metrics=[
          DiskLocator(
              expected_position=(100, 100),
              search_window=(30, 30),
              radius=10,
              radius_tolerance=2,
          )
      ]
  )
  img.plot()

This will search for a disk (BB) in the image at the expected position and window size for a disk of a given radius and tolerance.
If the disk is found, the location will be returned as a :class:`~pylinac.core.geometry.Point` object.
If the disk is not found, a ``ValueError`` will be raised.

The :class:`~pylinac.core.metrics.DiskRegion` metric is similar, but instead of returning the location, it returns a
`scikit-image regionprops <https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.regionprops>`__ object that is the region of the disk.
This allows one to then calculate things like the weighted centroid, area, etc.

Using physical units
####################

While pixels are useful, it is sometimes easier to use physical units.

To perform the same Disk/BB location using mm instead of pixels:

.. code-block:: python
  :caption: Search for a disk 30mm right and 30mm down from the top left of the image

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  img.compute(
      metrics=[
          # these are all in mm
          DiskLocator.from_physical(
              expected_position_mm=(30, 30),
              search_window_mm=(10, 10),
              radius_mm=4,
              radius_tolerance_mm=2,
          )
      ]
  )
  img.plot()


Relative to center
##################

We can also specify the expected position relative to the center of the image.

.. important::

  We can do this using pixels OR physical units.

This will look for the disk/BB 30 pixels right and 30 pixels down from the center of the image:

.. code-block:: python
  :caption: Relative to center using pixels

  from pylinac.core.image import DicomImage
  from pylinac.core.metrics import DiskLocator, DiskRegion

  img = DicomImage("my_image.dcm")
  img.compute(
      metrics=[
          # these are all in pixels
          DiskLocator.from_center(
              expected_position=(30, 30),
              search_window=(10, 10),
              radius=4,
              radius_tolerance=2,
          )
      ]
  )
  img.plot()

This will look for the disk/BB 30mm right and 30mm down from the center of the image:

.. code-block:: python
  :caption: Relative to center using physical units

  img.compute(
      metrics=[
          # these are all in mm
          DiskLocator.from_center_physical(
              expected_position_mm=(30, 30),
              search_window_mm=(10, 10),
              radius_mm=4,
              radius_tolerance_mm=2,
          )
      ]
  )
  img.plot()

Global Disk Locator
^^^^^^^^^^^^^^^^^^^

.. versionadded:: 3.17

The :class:`~pylinac.core.metrics.GlobalDiskLocator` metric is similar to the :class:`~pylinac.core.metrics.DiskLocator` metric
except that it searches the entire image for disks/BB, not just a small window. This is useful for finding the BB in images
where the BB is not in the expected location or unknown. This is also efficient for finding BBs in images,
even if the locations are known.

For example, here is an example analysis of an MPC image:

.. code-block:: python

  from pylinac.core.image import XIM
  from pylinac.core.metrics import GlobalDiskLocator

  img = XIM("my_image.xim")
  bbs = img.compute(
      metrics=GlobalDiskLocator(
          radius_mm=3.5,
          radius_tolerance_mm=1.5,
          min_number=10,
      )
  )
  img.plot()

This will result in an image like so:

.. image:: ../images/global_disk_locator.png
  :width: 600
  :align: center

.. _global_sized_field_locator:

Global Sized Field Locator
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 3.17

The :class:`~pylinac.core.metrics.GlobalSizedFieldLocator` metric is similar to the :class:`~pylinac.core.metrics.GlobalDiskLocator` metric
except that it searches the entire image for fields of a given size. This is useful for finding one or more fields in images
where the field is not in the expected location or unknown. This is also efficient when multiple fields are present in the image.

The locator will find the weighted center of the field(s) and return the location(s) as a :class:`~pylinac.core.geometry.Point` objects.
The boundary of the detected field(s) will be plotted on the image in addition to the center.

The locator will use pixels by default, but also has a ``from_physical`` class method to use physical units.

An example plot of finding multiple fields can be seen below:

.. image:: ../images/global_sized_field_locator.png
  :width: 600
  :align: center

For example:

.. code-block:: python
   :caption: Search for at least 2 fields of size 30x30 pixels with a tolerance of 4 pixels & plot

   img = DicomImage("my_image.dcm")
   img.compute(
       metrics=GlobalSizedFieldLocator(
           field_width_px=30, field_height_px=30, field_tolerance_px=4, max_number=2
       )
   )
   img.plot()  # this will plot the image with the fields overlaid

Using physical units
####################

To perform a similar field location using mm instead of pixels:

.. code-block:: python
   :caption: Search for at least 2 fields of size 30x30mm with a tolerance of 4mm

   img = DicomImage("my_image.dcm")
   img.compute(
       metrics=GlobalSizedFieldLocator.from_physical(
           field_width_mm=30, field_height_mm=30, field_tolerance_mm=4, max_number=2
       )
   )

Usage tips
##########

* Whenever possible, set the ``max_number`` parameter. This can **greatly** speed up the computation for several reasons.
  First, it will stop searching once the number of fields is found. Second, the thresholding algorithm will have a much
  better initial guess and also a better step size. This is because the approximate area of the field is known relative
  to the total image size.
* The ``field_tolerance_<mm|px>`` parameter can be relatively tight if the ``max_number`` parameter is set. Without a
  ``max_number`` parameter, you may have to increase the field tolerance to find all fields.

Writing Custom Plugins
----------------------

The power of the plugin architecture is that you can write your own metrics and use them on any image
as well as reuse them where needed.

To write a custom plugin, you must

* Inherit from the :class:`~pylinac.core.metrics.MetricBase` class
* Specify a ``name`` attribute.
* Implement the ``calculate`` method.
* (Optional) Implement the ``plot`` method if you want the metric to plot on the image.

.. warning::

    Do not modify the image in the ``calculate`` method as this will affect the image for other metrics and/or plotting.

For example, let's built a simple plugin that finds and plots an "X" at the center of the image:

.. plot::

    from pylinac.core.image_generator import AS1000Image, FilteredFieldLayer, GaussianFilterLayer
    from pylinac.core.image import DicomImage
    from pylinac.core.metrics import MetricBase

    class ImageCenterMetric(MetricBase):
        name = "Image Center"

        def calculate(self):
            return self.image.center

        def plot(self, axis: plt.Axes):
            axis.plot(self.image.center.x, self.image.center.y, 'rx', markersize=10)

    # now we create an image to compute over
    as1000 = AS1000Image(sid=1000)  # this will set the pixel size and shape automatically
    as1000.add_layer(
        FilteredFieldLayer(field_size_mm=(100, 100))
    )  # create a 100x100mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=2)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    ds = as1000.as_dicom()

    # now we can compute the metric on the image
    img = DicomImage.from_dataset(ds)
    center = img.compute(metrics=ImageCenterMetric())
    print(center)
    img.plot()

API
---

.. autoclass:: pylinac.core.metrics.MetricBase
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.metrics.DiskLocator
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.metrics.DiskRegion
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.metrics.GlobalDiskLocator
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.metrics.GlobalSizedFieldLocator
    :inherited-members:
    :members:
