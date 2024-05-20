.. _exporting-results:

=================
Exporting Results
=================

After analyzing a dataset in pylinac, you may want to export the data to a file for further analysis or for record-keeping.
Pylinac has built-in methods for **most** analyses to export to PDF, save figures, and dump data to string or JSON.


Exporting to PDF
----------------

For most analyses, the ``publish_pdf`` method is available. This method will create a PDF report of the analysis. The PDF will contain a summary of the analysis, the figure(s) generated, and any other relevant information.

.. code-block:: python

    from pylinac import WinstonLutz

    wl = WinstonLutz.from_zip("myzip.zip")
    wl.analyze()
    wl.publish_pdf("my_report.pdf")  # will save the report to a PDF file

You can even open this upon saving by setting the ``open_file`` parameter to ``True``.

.. code-block:: python

    ...
    wl.publish_pdf(
        "my_report.pdf", open_file=True
    )  # will save the report to a PDF file and open it

Exporting Figures
-----------------

If you want to save the figures generated during an analysis, you can use the ``save_analyzed_image`` method. This method will save the figure to a file.

.. code-block:: python

    from pylinac import Starshot

    star = Starshot.from_zip("myzip.zip")
    star.analyze()
    star.save_analyzed_image("my_starshot.png")  # will save the figure to a PNG file

Exporting to JSON
-----------------

If you want to save the data to a JSON file, you can use the ``results_data(to_json=True)`` method. This method will return a JSON string of the data.
This string can be saved to a file or otherwise passed around.

.. code-block:: python

    from pylinac import Starshot

    star = Starshot.from_zip("myzip.zip")
    star.analyze()
    json_str = star.results_data(
        to_json=True
    )  # dumps to a string that can be loaded by any JSON reader
    with open("my_starshot.json", "w") as f:
        f.write(json_str)

Exporting to QuAAC
------------------

.. note::

   `QuAAC <https://quaac.readthedocs.io/en/latest/index.html>`__ is a new standard we designed to allow for interoperable data exchange between vendors and clinics.
   It is in in the early stages of development and is not yet widely adopted.

To save results from an analysis to the QuAAC format, use the ``to_quaac()`` method. This will save the
data to a file in the QuAAC format.

.. code-block:: python

    from pylinac import Starshot
    from quaac import User, Equipment

    star = Starshot.from_zip("myzip.zip")
    star.analyze()
    performer = User(name='Michael Bluth', email='mbluth@sitwell.com')
    equipment = Equipment(name='Halcyon 1', type='Linac', serial_number='1234', manufacturer="Varian", model='Halcyon')
    star.to_quaac('star_quaac.yaml', format='yaml', performer=performer, primary_equipment=equipment)

It is important to know that QuAAC can store binary data like plotted images and PDFs that pylinac might produce.
This is not done by default due to the many options of publishing the PDF and plotting the images.
However, this can be done by separately saving the images and PDFs and then adding them to the QuAAC file.
as attachments.

.. code-block:: python

    from pylinac import Starshot
    from quaac import User, Equipment, Attachment

    star = Starshot.from_zip("myzip.zip")
    star.analyze()
    # save the images and PDFs
    star.save_analyzed_image('starshot.png')
    star.publish_pdf('starshot_report.pdf')
    # now dump to QuAAC
    performer = User(name='Michael Bluth', email='mbluth@sitwell.com')
    equipment = Equipment(name='Halcyon 1', type='Linac', serial_number='1234', manufacturer="Varian", model='Halcyon')
    png_attachment = Attachment.from_file(path='starshot.png', name='Starshot Image')
    pdf_attachment = Attachment.from_file(path='starshot_report.pdf', name='Starshot Report')
    original_image = Attachment.from_file(path='myzip.zip', name='Raw Image')
    star.to_quaac('star_quaac.yaml', format='yaml', performer=performer, primary_equipment=equipment, attachments=, attachments=[png_attachment, pdf_attachment, original_image])

See more `here <https://quaac.readthedocs.io/en/latest/writing_quaac.html#writing-quaac-files>`__ for writing QuAAC files.
