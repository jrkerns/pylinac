.. _image_conversion_guide:

======================
Image Conversion Guide
======================

This guide will provide a starting point based on the data in hand and desired data formats.
The point is to provide a quick reference for the user to convert their data to a format that Pylinac can use.

.. note::

    This does not address the "direct" use of files in Pylinac such as picket fence, starshot, etc.
    If you are just using Pylinac to analyze images, see the relevant module documentation. The point
    here is to apply intermediate manipulations, convert from one format to the other or to create new
    images.

+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| I have a...                                        |  Load with                                  | Export as          | Export call                                     |
+====================================================+=============================================+====================+=================================================+
| DICOM dataset                                      | ``DICOMImage.from_dataset(...)``            | DICOM file         | ``<>.as_dicom(...)``                            |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| TIFF/JPG image                                     | ``FileImage(path=...)``                     | DICOM Dataset      | ``<>.as_dicom()`` or ``tiff_to_dicom(...)``     |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| DICOM image                                        |  ``DICOMImage(path=...)``                   | DICOM Dataset      | ``<>.as_dicom()``                               |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
|                                                    |                                             | DICOM file         | ``<>.save(...)                                  |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| Multiple DICOM files representing one final image  | ``DICOMImage.from_multiples(filelist=...)`` |                    |                                                 |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| XIM Image                                          | ``XIM(...)``                                | PNG, JPG           | ``<>.save_as(...)``                             |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
|                                                    |                                             | DICOM Dataset      | ``<>.as_dicom(...)``                            |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| 2D Numpy array representing an image               | ``ArrayImage(array=...)``                   | DICOM Dataset      | ``<>.as_dicom(...)`` or ``array_to_dicom(...)`` |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| 3D Numpy array representing a volumetric scan (CT) |                                             | Set of DICOM files | ``create_dicom_files_from_3d_array(...)``       |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| Cyberknife RAW image                               | ``load_raw_cyberknife(path=...)``           |                    |                                                 |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+
| VisionRT RAW image                                 | ``load_raw_visionrt(path=...)``             |                    |                                                 |
+----------------------------------------------------+---------------------------------------------+--------------------+-------------------------------------------------+

Converting to DICOM
-------------------

As the above table shows, there are multiple ways to convert from a given format to DICOM.
DICOM is the preferred format of Pylinac and is obviously useful in other contexts as well.

When it comes to converting to DICOM in pylinac, the intent is create datasets and files
that can be used by Pylinac. While required tags are added, some keys may or may not be
present in the final DICOM file. Here we will address issues and considerations when using
the above flows.

* **DICOM images are written as RT Images** - Pylinac deals with RT Images. Written files
  will be RT Images unless the tag is overridden.
* **Pylinac can convert floating point arrays** - Pylinac will convert floating point arrays
  to uint16 arrays using :func:`~pylinac.core.array_utils.convert_to_dtype`.
* **Pylinac will not touch uint data types** - If you have a uint8 or uint16 array, Pylinac
  will not change the data type. And, although the DICOM standard does not officially support it
  uint32 arrays are also allowed.

Overriding Tags
---------------

When saving images to DICOM, you may override tags as needed. This is done by passing a dictionary
of the tags and their values.

.. warning::

    Passing a tag here will override the tag regardless of Pylinac's default or computed value might be.
    E.g. passing the ``PixelData`` tag will overwrite the image data.

.. code-block:: pylinac

    img = FileImage(path='my_image.tif')
    img.as_dicom(extra_tags={'PatientID': '12345', 'PatientName': 'John Doe'})
