
.. _xim-images:

XIM images
----------

Images ending in ``.xim`` are generally produced by a Varian TrueBeam or newer linac. They are images with additional
tags. Unfortunately, they are written in binary into a custom format so using a typical image library will not work.

The binary file specification appears to be unofficial, but it does work. You can find the spec `here <https://bitbucket.org/dmoderesearchtools/ximreader/raw/4900d324d5f28f8b6b57752cfbf4282b778a4508/XimReader/xim_readme.pdf>`__
which comes from this repo: https://bitbucket.org/dmoderesearchtools/ximreader/src/master/

.. warning::

    Rant ahead.

The XIM images used a custom compression format. Why they chose to use a custom format is beyond me. Moreso, the
format they chose was that of a PNG algorithm, but not as good. So, XIM images are just PNG images but with a custom lookup table
and property tags. Everyday PNG format would've worked just as well. It's possible this is security by obscurity or NIH syndrome.

Loading an XIM image
^^^^^^^^^^^^^^^^^^^^

To load an XIM images use the :class:`~pylinac.core.image.XIM` class:

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    # plot the image
    xim_img.plot()

    # see the XIM properties
    print(xim_img.properties)

Reconstructing the image pixels is relatively slow (~1s for AS1200 image) thanks to the custom compression format,
so if you are only searching through the properties you can skip reconstructing the pixels. Skipping the
pixels and only reading the properties is relatively fast (order of milliseconds):

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_files = [r"C:\TDS\H12345\QA\image.xim", ...]
    files_to_analyze = []
    for file in my_xim_files:
        # will load relatively fast
        xim_img = XIM(file, read_pixels=False)
        if xim_img.properties["AcquisitionMode"] == "Highres":
            files_to_analyze.append(file)

    # now load the pixel data only for the files we're interested in
    for file in files_to_analyze:
        xim_img = XIM(file)
        # image is available, do what you want
        xim_img.plot()

An XIM image has all the utility methods other pylinac images do, so use this to your advantage:

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    # process
    xim_img.crop(pixels=30)
    xim_img.filter()
    xim_img.fliplr()
    ...

.. _export-xim:

Exporting images
^^^^^^^^^^^^^^^^

Exporting ``*.xim`` images is easy. The PNG format is recommended because its ~1/2 the size of the original xim image and will
also include the properties. PNG is also lossless, so all information is retained.
PNG images can usually be viewed easily across many devices and OSs and also loads very fast.

.. code-block:: python

    from pylinac.core.image import XIM


    my_xim_file = r"C:\TDS\H12345\QA\image.xim"
    xim_img = XIM(my_xim_file)

    xim_img.save_as("myxim.png")
    # saved to PNG!

.. _reading-exported-xim:

Reading exported images
^^^^^^^^^^^^^^^^^^^^^^^

To load the image in python you can use any library that reads PNG. Pillow is recommended.
Opening these files are usually very fast (order of milliseconds), so
if you plan on doing research or analysis of a large number of .xim images, it may be worth it
to export to PNG en masse and then perform the analysis.

.. code-block:: python

    import numpy as np
    import PIL.Image
    import matplotlib.pyplot as plt

    xim_img = PIL.Image.open("myxim.png")

    # numpy array of the pixels
    xim_array = np.asarray(xim_img)

    # plot it
    plt.imshow(xim_array)
    plt.show()

To read the properties of an XIM file that was saved to PNG we may to have to load from strings.
PNG tags are all strings, and some xim properties are arrays or numbers. In order to
easily save it, we convert them all to strings. In order to get the native datatype
for non-string types we cast to the inferred type. For numbers, use ``float`` and for lists use ``json``:

.. code-block:: python

    import json
    import PIL.Image

    xim_img = PIL.Image.open("myxim.png")

    system_version = xim_img.info["AcquisitionSystemVersion"]
    # "2.7.304.16" already a string so no change needed

    couch_lat = xim_img.info["CouchLat"]
    # '100.39021332'  it's a string even though it looks like a number
    # convert to the original type:
    couch_lat_num = float(couch_lat)

    # MLCs are a list; we need json
    mlc_a_string = xim_img.info["MLCLeafsA"]
    # '[20.6643, 20.6992, ...]'
    mlc_a_list = json.loads(mlc_a_string)
    # now it's a normal list: [20.6643, 20.6992, ...]
