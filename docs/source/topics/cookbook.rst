.. _cookbook:

========
Cookbook
========

This is a collection of tips and tricks that can be useful when working with pylinac.

Combining DICOM files into a single file
----------------------------------------

There are cases when multiple DICOM files are related and should be treated as a single file.

.. code-block:: python

   from pylinac.image import load_multiples, load

   files = ["file1.dcm", "file2.dcm", "file3.dcm"]
   combined_img = load_multiples(files, method="mean", stretch_each=True, loader=load)
   # save to a new DICOM file
   combined_img.save("combined.dcm")
   # plot, etc
   combined_img.plot()

See :func:`~pylinac.core.image.load_multiples` for more information.

Filtering DICOM files from a directory
--------------------------------------

This will return only DICOM files from a directory, assuming there are other files (PDFs, text files, etc) in the directory.

.. code-block:: python
   :caption: Looking for DICOM files **only** at the root directory

   from pathlib import Path
   from pylinac.core.io import is_dicom

   dicom_files = [f for f in Path("my_directory").iterdir() if f.is_file() and is_dicom(f)]

.. code-block:: python
   :caption: Looking for DICOM files in this and all subdirectories

    import os
    from pathlib import Path
    from pylinac.core.io import is_dicom

    dicoms = []
    for dirpath, _, filenames in os.walk("my_directory"):
        for f in filenames:
            full_path = Path(dirpath) / f
            if is_dicom(full_path):
                dicoms.append(full_path)

See :func:`~pylinac.core.io.is_dicom` for more information.

Filtering DICOM images only
---------------------------

This will find DICOM files that are images specifically. E.g. pulling out RT images from amongst RT plans, etc.

.. code-block:: python

   from pathlib import Path
   from pylinac.core.io import is_dicom_image

   dicom_image_files = [
       f for f in Path("my_directory").iterdir() if f.is_file() and is_dicom_image(f)
   ]

See :func:`~pylinac.core.io.is_dicom_image` for more information.

Temporarily Unzip Files
-----------------------

It's quite common that you may want to analyze files inside a zip archive. Assuming that no modification
is required and that you want to leave the zip archive as you found it, you can use the following code to
perform actions in a temporary directory:

.. code-block:: python

   from pathlib import Path
   from pylinac.core.io import TemporaryZipDirectory

   zip_file = "my_stuff.zip"
   # contains a.dcm and b.txt
   with TemporaryZipDirectory(zip_file) as zip_dir:
       pf = PicketFence(Path(zip_dir / "a.dcm"))
       ...


Converting TIFF to DICOM
------------------------

See :ref:`tiff-to-dicom`.


Manipulating a DICOM file and saving back
-----------------------------------------

This use case is to modify a DICOM file and save it back to disk. This might be to apply known corrections before
analysis, etc.

.. code-block:: python

   from pylinac.image import DicomImage

   img = DicomImage("my_dicom.dcm")
   img.rotate(45)
   img.normalize()
   img.filter()
   img.save("my_modified_dicom.dcm")
   # my_modified_dicom.dcm can be treated like any other DICOM file

See :ref:`how-image-data-is-loaded` and :class:`~pylinac.core.image.DicomImage` for available manipulations.
