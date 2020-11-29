
.. _log_analyzer_module:

=================================
Log Analyzer module documentation
=================================

Overview
--------

.. automodule:: pylinac.log_analyzer
    :no-members:

.. _log_concepts:

Concepts
--------

Because the ``log_analyzer`` module functions without an end goal, the data has been formatted for easy exploration. However, there are a
few concepts that should be grasped before diving in.

* **Log Sections** - Upon log parsing, all data is placed into data structures. Varian has designated 4 sections for Trajectory logs:
  *Header*, *Axis Data*, *Subbeams*, and *CRC*. The *Subbeams* are only applicable for auto-sequenced beams and all v3.0 logs, and the *CRC* is specific to
  the Trajectory log. The *Header* and *Axis Data* however, are common to both Trajectory logs and Dynalogs.

   .. note::
    Dynalogs do not have explicit sections like the Trajectory logs,
    but pylinac formats them to have these two data structures for consistency.

* **Leaf Indexing & Positions** - Varian leaf identification is 1-index based, over against Python's 0-based indexing. Thus,
  indexing the first MLC leaf would be ``[1]``, not ``[0]``.

    .. warning:: When slicing or analyzing leaf data, keep the Varian 1-index base in mind.

  Leaf data is stored in a dictionary, with the leaf number as the key, from 1 up to the number of MLC leaves. E.g. if the machine has a
  Millennium 120 standard MLC model, leaf data will have 120 dictionary items from 1 to 120. Leaves from each bank have an offset of half the
  number of leaves. I.e. leaves A1 and B1 = 1 and 61. Thus, leaves 61-120 correspond
  to the B-bank, while leaves 1-60 correspond to the A-bank. This can be described by a function
  :math:`(A_{n}, B_{n}) = (n, n + N_{leaves}/2)`, where :math:`n` is the leaf number and :math:`N_{leaves}` is the number of leaves.

* **Units** - Units follow the Trajectory log specification: linear axes are in cm, rotational axes in degrees, and MU for dose.

  .. note::
    Dynalog files are inherently in mm for collimator and gantry axes, tenths of degrees for rotational axes, and
    MLC positions are not at isoplane. For consistency, Dynalog values are converted to Trajectory log specs, meaning
    linear axes, both collimator and MLCs are in cm at isoplane, and rotational axes are in degrees. Dynalog MU is always
    from 0 to 25000 no matter the delivered MU (i.e. it's relative), unless it was a VMAT delivery, in which case the
    gantry position is substituted in the dose fraction column.

  .. warning::
    Dynalog VMAT files replace the dose fraction column with the gantry position. Unfortunately, because of the variable dose rate of Varian linacs the gantry position
    is not a perfect surrogate for dose, but there is no other choice. Thus, fluence calculations will use the relative gantry movement as the dose in fluence calculations.


* **All data Axes are similar** - Log files capture machine data in "control cycles", aka "snapshots" or "heartbeats". Let's assume a
  log has captured 100 control cycles. Axis data that was captured will all be similar (e.g. gantry, collimator, jaws). They will all have an *actual* and sometimes an
  *expected* value for each cycle. Pylinac formats these as 1D numpy arrays along with a difference array if applicable. Each of these
  arrays can be quickly plotted for visual analysis. See :class:`~pylinac.log_analyzer.Axis` for more info.

Running the Demos
-----------------

As usual, the module comes with demo files and methods:

.. code-block:: python

    from pylinac import Dynalog
    Dynalog.run_demo()

Which will output the following::

    Results of file: C:\Users\James\Dropbox\Programming\Python\Projects\pylinac\pylinac\demo_files\AQA.dlg
    Average RMS of all leaves: 0.037 cm
    Max RMS error of all leaves: 0.076 cm
    95th percentile error: 0.088 cm
    Number of beam holdoffs: 20
    Gamma pass %: 18.65
    Gamma average: 0.468

.. plot::
    :include-source: false

    import pylinac
    pylinac.Dynalog.run_demo()

Your file location will be different, but the values should be the same.
The same can be done using the demo Trajectory log:

.. code-block:: python

    from pylinac import TrajectoryLog
    TrajectoryLog.run_demo()

Which will give::

    Results of file: C:\Users\James\Dropbox\Programming\Python\Projects\pylinac\pylinac\demo_files\Tlog.bin
    Average RMS of all leaves: 0.001 cm
    Max RMS error of all leaves: 0.002 cm
    95th percentile error: 0.002 cm
    Number of beam holdoffs: 19
    Gamma pass %: 100.00
    Gamma average: 0.002

.. plot::
    :include-source: false

    import pylinac
    pylinac.TrajectoryLog.run_demo()


Note that you can also save data in a PDF report:

.. code-block:: python

    tlog = ...
    tlog.publish_pdf('mytlog.pdf')

Loading Data
------------

Loading Single Logs
^^^^^^^^^^^^^^^^^^^

Logs can be loaded two ways.
The first way is through the main helper function :func:`~pylinac.log_analyzer.load_log`.

.. note::

    If you've used pylinac versions <1.6
    the helper function is new and can be a replacement for ``MachineLog`` and ``MachineLogs``, depending on the context
    as discussed below.

.. code-block:: python

    from pylinac import load_log

    log_path = "C:/path/to/tlog.bin"
    log = load_log(log_path)

In addition, a folder, ZIP archive, or URL can also be passed:

.. code-block:: python

    log2 = load_log('http://myserver.com/logs/2.dlg')

.. note:: If loading from a URL the object can be a file or ZIP archive.

Pylinac will automatically infer the log type and load it into the appropriate data structures for analysis.
The :func:`~pylinac.log_analyzer.load_log` function is a convenient wrapper around the classes within the log analysis module.
However, logs can be instantiated a second way: directly through the classes.

.. code-block:: python

    from pylinac import Dynalog, TrajectoryLog

    dlog_path = "C:/path/to/dlog.dlg"
    dlog = Dynalog(dlog_path)

    tlog_path = "C:/path/to/tlog.bin"
    tlog = TrajectoryLog(tlog_path)

Loading Multiple Logs
^^^^^^^^^^^^^^^^^^^^^

Loading multiple files is also possible using the :func:`~pylinac.log_analyzer.load_log` function as listed above.
The logs can also be directly instantiated by using :class:`~pylinac.log_analyzer.MachineLogs`. Acceptable inputs include a folder and zip archive.

.. code-block:: python

    from pylinac import load_log, MachineLogs

    path_to_folder = "C:/path/to/dir"

    # from folder; equivalent
    logs = MachineLogs(path_to_folder)
    logs = load_log(path_to_folder)

    # from ZIP archive
    logs = load_log('path/to/logs.zip')


Working with the Data
---------------------

Working with the log data is straightforward once the data structures and Axes are understood
(See :ref:`log_concepts` for more info). Pylinac follows the data structures specified by Varian for
trajectory logs, with a *Header* and *Axis Data* structure, and possibly a *Subbeams* structure if
the log is a Trajectory log and was autosequenced. For accessible attributes, see :class:`TrajectoryLog`. The
following sections explore each major section of log data and the data structures pylinac creates to assist
in data analysis.

.. note:: It may be helpful to also read the log specification format in parallel with this guide. It is easier to
 see that pylinac follows the log specifications and where the info comes from. Log specifications are on MyVarian.com.

Working with the Header
^^^^^^^^^^^^^^^^^^^^^^^

Header information is essentially anything that isn't axis measurement data; it's metadata about the file, format,
machine configuration, etc. Because of the different file formats, there are separate classes for Trajectory log and
Dynalog headers. The classes are:

    * :class:`~pylinac.log_analyzer.TrajectoryLogHeader`
    * :class:`~pylinac.log_analyzer.DynalogHeader`

Header attributes are listed in the class API docs by following the above links. For completeness they are also listed
here. For Trajectory logs:

    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.header`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.version`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.header_size`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.sampling_interval`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.num_axes`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.axis_enum`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.samples_per_axis`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.num_mlc_leaves`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.axis_scale`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.num_subbeams`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.is_truncated`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.num_snapshots`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogHeader.mlc_model`

For Dynalogs the following header information is available:

    * :attr:`~pylinac.log_analyzer.DynalogHeader.version`
    * :attr:`~pylinac.log_analyzer.DynalogHeader.patient_name`
    * :attr:`~pylinac.log_analyzer.DynalogHeader.plan_filename`
    * :attr:`~pylinac.log_analyzer.DynalogHeader.tolerance`
    * :attr:`~pylinac.log_analyzer.DynalogHeader.num_mlc_leaves`
    * :attr:`~pylinac.log_analyzer.DynalogHeader.clinac_scale`

.. rubric:: Example

Let's explore the header of the demo trajectory log:

.. code-block:: python

    >>> tlog = TrajectoryLog.from_demo()
    >>> tlog.header.header
    'VOSTL'
    >>> tlog.header.version
    2.1
    >>> tlog.header.num_subbeams
    2

Working with Axis Data
^^^^^^^^^^^^^^^^^^^^^^

Axis data is all the information relating to the measurements of the various machine axes and is accessible
under the ``axis_data`` attribute. This includes the gantry,
collimator, MLCs, etc. Trajectory logs capture more information than Dynalogs, and additionally hold the expected
positions not only for MLCs but also for all axes. Every measurement axis has :class:`~pylinac.log_analyzer.Axis` as
its base; they all have similar methods to access and plot the data (see :ref:`plotting`). However, not all attributes
are axes. Pylinac adds properties to the axis data structure for ease of use (e.g. the number of snapshots)
For Trajectory logs the following attributes are available,
based on the :class:`~pylinac.log_analyzer.TrajectoryLogAxisData` class:

    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.collimator`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.gantry`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.jaws`

      .. note::
            The ``jaws`` attribute is a data structure to hold all 4 jaw axes; see
            :class:`~pylinac.log_analyzer.JawStruct`

    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.couch`

      .. note::
            The ``couch`` attribute is a data structure to hold lateral, longitudinal, etc couch positions; see
            :class:`~pylinac.log_analyzer.CouchStruct`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.mu`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.beam_hold`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.control_point`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.carriage_A`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.carriage_B`
    * :attr:`~pylinac.log_analyzer.TrajectoryLogAxisData.mlc`

      .. note::
        The ``mlc`` attribute is a data structure to hold leaf information; see
        :class:`~pylinac.log_analyzer.MLC` for attributes and the :ref:`mlc` section for more info.

Dynalogs have similar attributes, derived from the :class:`~pylinac.log_analyzer.DynalogAxisData` class:

    * :attr:`~pylinac.log_analyzer.DynalogAxisData.collimator`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.gantry`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.jaws`

      .. note::
        The ``jaws`` attribute is a data structure to hold all 4 jaw axes; see
        :class:`~pylinac.log_analyzer.JawStruct`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.num_snapshots`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.mu`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.beam_hold`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.beam_on`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.previous_segment_num`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.previous_dose_index`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.next_dose_index`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.carriage_A`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.carriage_B`
    * :attr:`~pylinac.log_analyzer.DynalogAxisData.mlc`

      .. note::
        The ``mlc`` attribute is a data structure to hold leaf information; see
        :class:`~pylinac.log_analyzer.MLC` for attributes and the :ref:`mlc` section for more info.

.. rubric:: Example

Let's access a few axis data attributes:

.. code-block:: python

    >>> log = Dynalog.from_demo()
    >>> log.axis_data.mu.actual  # a numpy array
    array([  0, 100, ...
    >>> log.axis_data.num_snapshots
    99
    >>> log.axis_data.gantry.actual
    array([ 180, 180, 180, ...

.. _mlc:

Working with MLC Data
^^^^^^^^^^^^^^^^^^^^^

Although MLC data is acquired and included in Trajectory logs and Dynalogs, it is not always easy to parse. Additionally,
a physicist may be interested in the MLC metrics of a log (RMS, etc). Pylinac provides tools for accessing MLC raw data as
well as helper methods and properties via the :class:`~pylinac.log_analyzer.MLC` class. Note that this class is consistent
between Trajectory logs and Dynalogs. This class is reachable through the `axis_data`
attribute as ``mlc``.

Accessing Leaf data
###################

Leaf data for any leaf is available under the :attr:`~pylinac.log_analyzer.MLC.leaf_axes` attribute which is a dict. The leaves are keyed
by the leaf number and the value is an :class:`~pylinac.log_analyzer.Axis`. Example::

    >>> log = Dynalog.from_demo()
    >>> log.axis_data.mlc.leaf_axes[1].actual  # numpy array of the 'actual' values for leaf #1
    array([ 7.56374, ...
    >>> log.axis_data.mlc.leaf_axes[84].difference  # actual values minus the planned values for leaf 84
    array([-0.001966, ...

MLC helper methods/properties
#############################

Beyond direct MLC data, pylinac provides a number of helper methods and properties to make working with MLC data easier
and more helpful. All the methods are listed in the :class:`~pylinac.log_analyzer.MLC` class, but some examples of use
are given here::

    >>> log = Dynalog.from_demo()
    >>> log.axis_data.mlc.get_error_percentile(percentile=95)  # get an MLC error percentile value
    0.08847
    >>> log.axis_data.mlc.leaf_moved(12)  # did leaf 12 move during treatment?
    False
    >>> log.axis_data.mlc.get_RMS_avg()  # get the average RMS error
    0.03733
    >>> log.axis_data.mlc.get_RMS_avg('A')  # get the average RMS error for bank A
    0.03746
    >>> log.axis_data.mlc.num_leaves  # the number of MLC leaves
    120
    >>> log.axis_data.mlc.num_moving_leaves  # the number of leaves that moved during treatment
    60

.. _fluences:

Working with Fluences
^^^^^^^^^^^^^^^^^^^^^

Fluences created by the MLCs can also be accessed and viewed. Fluences are accessible under the `fluence`
attribute. There are three subclasses that handle the fluences:
The fluence actually delivered is in :class:`~pylinac.log_analyzer.ActualFluence`, the fluence planned is in
:class:`~pylinac.log_analyzer.ExpectedFluence`, and the gamma of the fluences is in :class:`~pylinac.log_analyzer.GammaFluence`.
Each fluence must be calculated, however pylinac makes reasonable defaults and has a few shortcuts. The actual and
expected fluences can be calculated to any resolution in the leaf-moving direction. Some examples::

    >>> log = Dynalog.from_demo()
    >>> log.fluence.actual.calc_map()  # calculate the actual fluence; returns a numpy array
    array([ 0, 0, ...
    >>> log.fluence.expected.calc_map(resolution=1)  # calculate at 1mm resolution
    array([ 0, 0, ...
    >>> log.fluence.gamma.calc_map(distTA=0.5, doseTA=1, resolution=0.1)  # change the gamma criteria
    array([ 0, 0, ...
    >>> log.fluence.gamma.pass_prcnt  # the gamma passing percentage
    99.82
    >>> log.fluence.gamma.avg_gamma  # the average gamma value
    0.0208

.. _plotting:

Plotting & Saving Axes/Fluences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each and every axis of the log can be accessed as a numpy array and/or plotted.
For each axis the "actual" array/plot is always available. Dynalogs only have expected values for the MLCs.
Trajectory logs have the actual and expected values for all axes. Additionally, if an axis has actual and
expected arrays, then the difference is also available.

Example of plotting the MU actual:

.. plot::
    :context:

    import pylinac

    log = pylinac.TrajectoryLog.from_demo()
    log.axis_data.mu.plot_actual()

Plot the Gantry difference:

.. plot::
    :context:

    log.axis_data.gantry.plot_difference()

Axis plots are just as easily saved:

.. code-block:: python

    log.axis_data.gantry.save_plot_difference(filename='gantry diff.png')

Now, lets plot the actual fluence:

.. plot::
    :context:

    log.fluence.actual.calc_map()
    log.fluence.actual.plot_map()

And the fluence gamma. But note we must calculate the gamma first, passing in any DoseTA or DistTA parameters:

.. plot::
    :context:

    log.fluence.gamma.calc_map()
    log.fluence.gamma.plot_map()

Additionally, you can calculate and view the fluences of subbeams if you're working with trajectory logs:

.. plot::
    :context:

    log.subbeams[0].fluence.gamma.calc_map()
    log.subbeams[0].fluence.actual.plot_map()

Converting Trajectory logs to CSV
---------------------------------

If you already have the log files, you obviously have a record of treatment. However, trajectory logs are in binary
format and are not easily readable without tools like pylinac. You can save trajectory logs in a more readable format
through the :meth:`~pylinac.log_analyzer.TrajectoryLog.to_csv()` method. This will write the log to a comma-separated
variable (CSV) file, which can be read with Excel and many other programs. You can do further or specialized analysis
with the CSV files if you wish, without having to use pylinac:

.. code-block:: python

    log = TrajectoryLog.from_demo()
    log.to_csv()

Anonymizing Logs
----------------

Machine logs can be anonymized two ways. The first is using the :meth:`~pylinac.log_analyzer.TrajectoryLog.anonymize` method, available to
both Trajectory logs and Dynalogs. Example script:

.. code-block:: python

    tlog = TrajectoryLog.from_demo()
    tlog.anonymize()
    dlog = Dynalog.from_demo()
    dlog.anonymize()

The other way is the use the module function :func:`~pylinac.log_analyzer.anonymize`. This function will anonymize a single
log file or a whole directory. If you plan on anonymizing a lot of logs, use this method as it is threaded and is much faster:

.. code-block:: python

    from pylinac.log_analyzer import anonymize

    log_file = 'path/to/tlog.bin'
    anonymize(log_file)
    log_dir = 'path/to/log/folder'
    anonymize(log_dir) # VERY fast

Batch Processing
----------------

Batch processing/loading of log files is helpful when dealing with one file at a time is too cumbersome. Pylinac allows you
to load logs of an entire directory via :class:`~pylinac.log_analyzer.MachineLogs`; individual log files can be accessed, and a handful of
batch methods are included.

.. rubric:: Example

Let's assume all of your logs for the past week are in a folder. You'd like to quickly see what the average gamma is of the files::

    >>> from pylinac import MachineLogs
    >>> log_dir = r"C:\path\to\log\directory"
    >>> logs = MachineLogs(log_dir)
    >>> logs.avg_gamma(resolution=0.2)
    0.03  # or whatever

You can also append to :class:`~pylinac.log_analyzer.MachineLogs` to have two or more different folders combined::

    >>> other_log_dir = r"C:\different\path"
    >>> logs.append(other_log_dir)

Trajectory logs in a MachineLogs instance can also be converted to CSV, just as for a single instance of TrajectoryLog::

    >>> logs.to_csv()  # only converts trajectory logs; dynalogs are already basically CSV files

.. note::
    Batch processing methods (like :meth:`~pylinac.log_analyzer.MachineLogs.avg_gamma` can take a while if numerous logs have been
    loaded, so be patient. You can also
    use the ``verbose=True`` argument in batch methods to see how the process is going.

API Documentation
-----------------

.. autofunction:: pylinac.log_analyzer.load_log

.. autofunction:: pylinac.log_analyzer.anonymize

.. autoclass:: pylinac.log_analyzer.Dynalog

.. autoclass:: pylinac.log_analyzer.TrajectoryLog

.. autoclass:: pylinac.log_analyzer.MachineLogs

.. autoclass:: pylinac.log_analyzer.Axis

.. autoclass:: pylinac.log_analyzer.MLC

.. autoclass:: pylinac.log_analyzer.DynalogHeader

.. autoclass:: pylinac.log_analyzer.DynalogAxisData

.. autoclass:: pylinac.log_analyzer.TrajectoryLogHeader

.. autoclass:: pylinac.log_analyzer.TrajectoryLogAxisData

.. autoclass:: pylinac.log_analyzer.SubbeamManager

.. autoclass:: pylinac.log_analyzer.Subbeam

.. autoclass:: pylinac.log_analyzer.FluenceStruct

.. autoclass:: pylinac.log_analyzer.FluenceBase

.. autoclass:: pylinac.log_analyzer.ActualFluence

.. autoclass:: pylinac.log_analyzer.ExpectedFluence

.. autoclass:: pylinac.log_analyzer.GammaFluence

.. autoclass:: pylinac.log_analyzer.JawStruct

.. autoclass:: pylinac.log_analyzer.CouchStruct

.. autoclass:: pylinac.log_analyzer.NotALogError

.. autoclass:: pylinac.log_analyzer.NotADynalogError

.. autoclass:: pylinac.log_analyzer.DynalogMatchError
