
=========================
CBCT module documentation
=========================

Overview
--------

.. automodule:: pylinac.cbctqa.cbct

Running the Demo
----------------

To run the CBCT demo, run cbct.py, or create a script and run::

    from pylinac import CBCT
    CBCT().run_demo_head()
    # You can run the other two demos just as easily as well
    # CBCT().run_demo_thorax()
    # CBCT().run_demo_pelvis()

Results will be printed to the console and a figure showing the slices analyzed will pop up::

    - CBCT QA Test -
    HU Regions:  {'PMP': -199.748046875, 'Teflon': 997.01953125, 'Air': -997.205078125, 'LDPE': -103.24596774193549, 'Acrylic': 116.6078431372549, 'Delrin': 340.61328125, 'Poly': -45.062745098039215}
    HU Passed?:  True
    Uniformity:  {'Top': 5.9677290836653389, 'Left': 10.206608280254777, 'Center': 21.318362480127185, 'Right': -0.3411624203821656, 'Bottom': 2.9124203821656049}
    Uniformity Passed?:  True
    MTF 50% (lp/mm):  1.11
    Geometric distances:  {'Top-Horiz': 49.695505782164226, 'Bottom-Horiz': 49.689076691833307, 'Right-Vert': 49.668186459198381, 'Left-Vert': 49.692334275563923}
    Geometry Passed?:  True

.. image:: /images/cbct_head.png
   :height: 550px
   :width: 550px


Typical Use
-----------

CBCT analysis as done by this module closely follows what is specified in the CatPhan manuals, replacing the need for hand measurements.
Assuming you've made a CBCT object as follows::

    from pylinac import CBCT
    mycbct = CBCT()

The minimum needed to get going is to:

* **Load images** -- Loading the DICOM images into your CBCT object can be done by passing the folder the images are located in.
 This can be done directly, or by using a UI. The code might look like either of the following::

    # set the folder path
    cbct_folder = r"C:/QA Folder/CBCT/June monthly"  # use of 'r' is for raw string; otherwise spaces and backslashes aren't interpreted properly
    # load the images from the file path
    mycbct.load_folder(cbct_folder)

    # *OR*

    # Identify the folder using a UI
    mycbct.load_folder_UI()


* **Analyze the images** -- Once the folder/images are loaded, tell CBCT to start analyzing the images. See the
  Algorithm section for details on how this is done::

    mycbct.analyze()

* **View the results** -- The CBCT module can print out the summary of results to the console as well as draw a matplotlib image to show where the
  samples were taken and their values::

      # print results to the console
      mycbct.return_results()
      # view analyzed images
      mycbct.plot_analyzed_image()

Algorithm
---------

The CBCT module is based on the tests and values given in the CatPhan 504 Manual. The algorithm works like such:

**Allowances**

* The images can be any size.
* For Varian machines, the images can be acquired with any protocol (Pelvis, Head, etc).
* The phantom can have significant translation in the Left-Right, Up-Down direction.
* The phantom can have significant rotation in any direction (yaw, pitch, and roll).

**Restrictions**

* The phantom used must be an unmodified CatPhan 504, as endorsed and supplied by Varian.
* The phantom must have <1cm offset in the In-Out direction (work to remove this is in the plans).

**Pre-Analysis**

* **Determine image properties** -- Upon load, the image set is analyzed for its DICOM properties to determine mm/pixel
  spacing, rescale intercept and slope, and manufacturer. All the images are resized to 512x512 pixels for streamlined
  analysis; this doesn't affect the physical spacing measurements.
* **Convert to HU** -- The entire image set is converted from its raw values to HU by applying the rescale intercept
  and slope which is contained in the DICOM properties.

**Analysis**

* **Determine phantom roll** -- Precise knowledge of the ROIs to analyze is important, and small changes in rotation
  could invalidate automatic results. The roll of the phantom is determined by examining the HU module and converting to
  binary. The air holes are then located and the angle of the two holes determines the phantom roll.

.. note::
    For each step below, the "module" analyzed is actually the mean, median, or maximum of 3 slices (+/-1 slice around and
    including the nominal slice) to ensure robust measurements. Also, for each step/module, the phantom center is
    determined, which corrects for the phantom pitch and yaw.

* **Determine HU linearity** -- The HU module (CTP404) contains several materials with different HU values. Using
  hardcoded angles (corrected for roll) and radius from the center of the phantom, circular ROIs are sampled which
  correspond to the HU material regions. The mean pixel values of these ROIs are the HU values.
* **Determine HU uniformity** -- HU uniformity (CTP486) is calculated in a very similar manner to HU linearity, but
  within the CTP486 module/slice.
* **Calculate Geometry/Scaling** -- The HU module (CTP404), besides HU materials, also contains several "nodes" which
  have an accurate spacing (5cm apart). Again, using hardcoded but corrected angles, the area around the 4 nodes are
  sampled and then a threshold is applied which identifies the node within the ROI sample. The center of this node is
  determined and then the space between nodes is calculated.
* **Calculate Spatial Resolution/MTF** -- The Spatial Resolution module (CTP528) contains 21 pairs of aluminum bars
  having varying thickness, which also corresponds to the thickness between the bars. One unique advantage of these
  bars is that they are all focused on and equally distant to the phantom center. Thus, several circular profiles are
  taken that go along all the line pairs. Starting on the proximal side of the line pairs, five profiles, each 1 pixel
  apart, sample an entire circle that cuts through the line pairs. The five profiles are concatenated and the median
  value extracted to effectively form one median profile. The peaks and valleys of the profile are located; peaks and
  valleys of the same line pair are averages. The relative MTF (i.e. normalized to the first line pair) is then
  calculated from these values.
* **Calculate Low Contrast Resolution** -- Not yet implemented, but it's being examined.

**Post-Analysis**

* **Test if values are within tolerance** -- For each module, the determined values are compared with the nominal values.
  If the difference between the two is below the specified tolerance then the module passes.


API Documentation
-----------------

The CBCT class is derived from a Base Class; inherited members are included in this list. More at :ref:`baseclass_api_doc`.

.. autoclass:: pylinac.cbctqa.cbct.CBCT
    :members:
    :inherited-members:

.. autoclass:: pylinac.cbctqa.cbct.Slice
