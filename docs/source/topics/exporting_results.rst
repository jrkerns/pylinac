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
