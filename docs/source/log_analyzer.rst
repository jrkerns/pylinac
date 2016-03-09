
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
  Millennium 120 standard MLC model, leaf data will have 120 dictionary items from 1 to 120. Leaf numbers have an offset of half the
  number of leaves. I.e. leaves 1 and 120 are a pair, as are 2 and 119, on up to leaves 60 and 61. In such a case, leaves 61-120 correspond
  to the B-bank, while leaves 1-60 correspond to the A-bank. This can be described by a function
  :math:`(A_{leaf}, B_{leaf}) = (n, N_{leaves} + 1 - n)`, where :math:`n` is the leaf number and :math:`N_{leaves}` is the number of leaves.

* **Units** - Units follow the Trajectory log specification: linear axes are in cm, rotational axes in degrees, and MU for dose.

  .. note::
    Dynalog files are inherently in mm for collimator and gantry axes, tenths of degrees for rotational axes, and
    MLC positions are not at isoplane. For consistency, Dynalog values are converted to Trajectory log specs, meaning
    linear axes, both collimator and MLCs are in cm at isoplane, and rotational axes are in degrees. Dynalog MU is always
    from 0 to 25000 no matter the delivered MU (i.e. it's relative), unless it was a VMAT delivery, in which case the MU is actually the gantry position.


* **All data Axes are similar** - Log files capture machine data in "control cycles", aka "snapshots" or "heartbeats". Let's assume a
  log has captured 100 control cycles. Axis data that was captured will all be similar (e.g. gantry, collimator, jaws). They will all have an *actual* and sometimes an
  *expected* value for each cycle. Pylinac formats these as 1D numpy arrays along with a difference array if applicable. Each of these
  arrays can be quickly plotted for visual analysis. See :class:`~pylinac.log_analyzer.Axis` for more info.

Running the Demo
----------------

As usual, the module comes with demo files and methods::

    >>> from pylinac import MachineLog
    >>> MachineLog.run_dlog_demo()

Which will output the following::

    MLC log type: Dynalog
    Average RMS of all leaves: 0.074 cm
    Max RMS error of all leaves: 0.076 cm
    95th percentile error: 0.088 cm
    Number of beam holdoffs: 20
    Gamma pass %: 99.83
    Gamma average: 0.021

.. image:: images/logs/dlog_results.png

The same can be done using the demo Trajectory log::

    >>> MachineLog.run_tlog_demo()
    MLC log type: Trajectory log
    Average RMS of all leaves: 0.001 cm
    Max RMS error of all leaves: 0.002 cm
    95th percentile error: 0.002 cm
    Number of beam holdoffs: 0
    Gamma pass %: 100.00
    Gamma average: 0.002

.. image:: images/logs/tlog_analyzed.png

Loading Data
------------

Logs are loaded upon class initialization::

    from pylinac import MachineLog

    log_path = "C:/path/to/tlog.bin"
    log = MachineLog(log_path)

Or load from a URL::

    log = MachineLog.from_url('http://myserver.com/logs/1.dlg')

When dealing with multiple logs, use :class:`~pylinac.log_analyzer.MachineLogs`::

    from pylinac import MachineLogs  # note the `s`

    log_folder = "path/to/folder"
    logs = MachineLogs(log_folder)

Or load from a ZIP archive::

    logs = MachineLogs.from_zip('path/to/logs.zip')

Pylinac will automatically infer the log type and load it into data structures for analysis.

Working with the Data
---------------------

Working with the log data is straightforward once the data structures and Axes are understood
(See :ref:`log_concepts` for more info). Pylinac follows the data structures specified by Varian for
trajectory logs, with a *Header* and *Axis Data* structure, and possibly a *Subbeams* structure if
the log is a Trajectory log and was autosequenced. For accessible attributes, see :class:`MachineLog`. The
following sections explore each major section of log data and the data structures pylinac creates to assist
in data analysis.

.. note:: It may be helpful to also read the log specification format in parallel with this guide. It is easier to
 see that pylinac follows the log specifications and where the info comes from. Log specifications are on MyVarian.com.

Working with the Header
^^^^^^^^^^^^^^^^^^^^^^^

Header information is essentially anything that isn't axis measurement data; it's metadata about the file, format,
machine configuration, etc. Because of the different file formats, there are separate classes for Trajectory log and
Dynalog headers. The classes are:

    * :class:`~pylinac.log_analyzer.TlogHeader`
    * :class:`~pylinac.log_analyzer.DlogHeader`

Header attributes are listed in the class API docs by following the above links. For completeness they are also listed
here. For Trajectory logs:

    * :attr:`~pylinac.log_analyzer.TlogHeader.header`
    * :attr:`~pylinac.log_analyzer.TlogHeader.version`
    * :attr:`~pylinac.log_analyzer.TlogHeader.header_size`
    * :attr:`~pylinac.log_analyzer.TlogHeader.sampling_interval`
    * :attr:`~pylinac.log_analyzer.TlogHeader.num_axes`
    * :attr:`~pylinac.log_analyzer.TlogHeader.axis_enum`
    * :attr:`~pylinac.log_analyzer.TlogHeader.samples_per_axis`
    * :attr:`~pylinac.log_analyzer.TlogHeader.num_mlc_leaves`
    * :attr:`~pylinac.log_analyzer.TlogHeader.axis_scale`
    * :attr:`~pylinac.log_analyzer.TlogHeader.num_subbeams`
    * :attr:`~pylinac.log_analyzer.TlogHeader.is_truncated`
    * :attr:`~pylinac.log_analyzer.TlogHeader.num_snapshots`
    * :attr:`~pylinac.log_analyzer.TlogHeader.mlc_model`

For Dynalogs the following header information is available:

    * :attr:`~pylinac.log_analyzer.DlogHeader.version`
    * :attr:`~pylinac.log_analyzer.DlogHeader.patient_name`
    * :attr:`~pylinac.log_analyzer.DlogHeader.plan_filename`
    * :attr:`~pylinac.log_analyzer.DlogHeader.tolerance`
    * :attr:`~pylinac.log_analyzer.DlogHeader.num_mlc_leaves`
    * :attr:`~pylinac.log_analyzer.DlogHeader.clinac_scale`

.. rubric:: Example

Let's explore the header of the demo trajectory log::

    >>> tlog = MachineLog.from_demo_trajectorylog()
    >>> tlog.header.header
    'VOSTL'
    >>> tlog.header.version
    2.1
    >>> tlog.header.num_subbeams
    2

Working with Axis Data
^^^^^^^^^^^^^^^^^^^^^^

Axis data is all the information relating to the measurements of the various machine axes and is accessible
under the :attr:`~pylinac.log_analyzer.MachineLog.axis_data` attribute. This includes the gantry,
collimator, MLCs, etc. Trajectory logs capture more information than Dynalogs, and additionally hold the expected
positions not only for MLCs but also for all axes. Every measurement axis has :class:`~pylinac.log_analyzer.Axis` as
its base; they all have similar methods to access and plot the data (see :ref:`plotting`). However, not all attributes
are axes. Pylinac adds properties to the axis data structure for ease of use (e.g. the number of snapshots)
For Trajectory logs the following attributes are available,
based on the :class:`~pylinac.log_analyzer.TlogAxisData` class:

    * :attr:`~pylinac.log_analyzer.TlogAxisData.collimator`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.gantry`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.jaws`

      .. note::
            The ``jaws`` attribute is a data structure to hold all 4 jaw axes; see
            :class:`~pylinac.log_analyzer.JawStruct`

    * :attr:`~pylinac.log_analyzer.TlogAxisData.couch`

      .. note::
            The ``couch`` attribute is a data structure to hold lateral, longitudinal, etc couch positions; see
            :class:`~pylinac.log_analyzer.CouchStruct`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.mu`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.beam_hold`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.control_point`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.carriage_A`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.carriage_B`
    * :attr:`~pylinac.log_analyzer.TlogAxisData.mlc`

      .. note::
        The ``mlc`` attribute is a data structure to hold leaf information; see
        :class:`~pylinac.log_analyzer.MLC` for attributes and the :ref:`mlc` section for more info.

Dynalogs have similar attributes, derived from the :class:`~pylinac.log_analyzer.DlogAxisData` class:

    * :attr:`~pylinac.log_analyzer.DlogAxisData.collimator`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.gantry`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.jaws`

      .. note::
        The ``jaws`` attribute is a data structure to hold all 4 jaw axes; see
        :class:`~pylinac.log_analyzer.JawStruct`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.num_snapshots`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.mu`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.beam_hold`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.beam_on`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.previous_segment_num`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.previous_dose_index`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.next_dose_index`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.carriage_A`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.carriage_B`
    * :attr:`~pylinac.log_analyzer.DlogAxisData.mlc`

      .. note::
        The ``mlc`` attribute is a data structure to hold leaf information; see
        :class:`~pylinac.log_analyzer.MLC` for attributes and the :ref:`mlc` section for more info.

.. rubric:: Example

Let's access a few axis data attributes::

    >>> log = MachineLog.from_demo_dynalog()
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
between Trajectory logs and Dynalogs. This class is reachable through the :attr:`~pylinac.log_analyzer.MachineLog.axis_data`
attribute as ``mlc``.

Accessing Leaf data
###################

Leaf data for any leaf is available under the :attr:`~pylinac.log_analyzer.MLC.leaf_axes` attribute which is a dict. The leaves are keyed
by the leaf number and the value is an :class:`~pylinac.log_analyzer.Axis`. Example::

    >>> log = MachineLog.from_demo_trajectorylog()
    >>> log.axis_data.mlc.leaf_axes[1].actual  # numpy array of the 'actual' values for leaf #1
    array([ 7.56374, ...
    >>> log.axis_data.mlc.leaf_axes[84].difference  # actual values minus the planned values for leaf 84
    array([-0.001966, ...

MLC helper methods/properties
#############################

Beyond direct MLC data, pylinac provides a number of helper methods and properties to make working with MLC data easier
and more helpful. All the methods are listed in the :class:`~pylinac.log_analyzer.MLC` class, but some examples of use
are given here::

    >>> log = MachineLog.from_demo_dynalog()
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

Fluences created by the MLCs can also be accessed and viewed. Fluences are accessible under the :class:`~pylinac.log_analyzer.MachineLog.fluence`
attribute. There are three subclasses that handle the fluences:
The fluence actually delivered is in :class:`~pylinac.log_analyzer.ActualFluence`, the fluence planned is in
:class:`~pylinac.log_analyzer.ExpectedFluence`, and the gamma of the fluences is in :class:`~pylinac.log_analyzer.GammaFluence`.
Each fluence must be calculated, however pylinac makes reasonable defaults and has a few shortcuts. The actual and
expected fluences can be calculated to any resolution in the leaf-moving direction. Some examples::

    >>> log = MachineLog.from_demo_dynalog()
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

Example of plotting the MU actual::

    log = MachineLog.from_demo_trajectorylog()
    log.axis_data.mu.plot_actual()

.. raw:: html
    :file: images/logs/tlog_mu_actual.html

Plot the Gantry difference::

    log.axis_data.gantry.plot_difference()

.. raw:: html
    :file: images/logs/gantry_difference.html

Axis plots are just as easily saved::

    log.axis_data.gantry.save_plot_difference(filename='gantry diff.png')

Now, lets plot the actual fluence::

    log.fluence.actual.plot_map()

.. plot::

    from pylinac import MachineLog
    log = MachineLog.from_demo_trajectorylog()
    log.fluence.actual.calc_map()
    log.fluence.actual.plot_map()

And the fluence gamma::

    log.fluence.gamma.plot_map()

.. image:: images/logs/log_gamma.png

Additionally, you can calculate and view the fluences of subbeams if you're working with trajectory logs::

    log = MachineLog.from_demo_trajectorylog()
    log.subbeams[0].fluence.actual.calc_map()
    log.subbeams[0].fluence.actual.plot_map()

.. plot::

    from pylinac import MachineLog
    log = MachineLog.from_demo_trajectorylog()
    log.subbeams[0].fluence.gamma.calc_map()
    log.subbeams[0].fluence.actual.plot_map()

Converting Trajectory logs to CSV
---------------------------------

If you already have the log files, you obviously have a record of treatment. However, trajectory logs are in binary
format and are not easily readable without tools like pylinac. You can save trajectory logs in a more readable format
through the :meth:`~pylinac.log_analyzer.MachineLog.to_csv()` method. This will write the log to a comma-separated
variable (CSV) file, which can be read with Excel and many other programs. You can do further or specialized analysis
with the CSV files if you wish, without having to use pylinac::

    log = MachineLog.from_demo_trajectorylog()
    log.to_csv()

Anonymizing Logs
----------------

Machine logs can be anonymized two ways. The first is using the :meth:`~pylinac.log_analyzer.MachineLog.anonymize` method, available to
both :class:`~pylinac.log_analyzer.MachineLog` and :class:`~pylinac.log_analyzer.MachineLogs`. Example script::

    log = MachineLog.from_demo_trajectorylog()
    log.anonymize()

The other way is the use the module function :func:`~pylinac.log_analyzer.anonymize`. This function will anonymize a single
log file or a whole directory. If you plan on anonymizing a lot of logs, use this method as it is threaded and is much faster::

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

Trajectory logs in a MachineLogs instance can also be converted to CSV, just as for a MachineLog::

    >>> logs.to_csv()  # only converts trajectory logs; dynalogs are already basically CSV files

.. note::
    Batch processing methods (like :meth:`~pylinac.log_analyzer.MachineLogs.avg_gamma` can take a while if numerous logs have been
    loaded, so be patient. You can also
    use the ``verbose=True`` argument in batch methods to see how the process is going.

API Documentation
-----------------

.. autoclass:: pylinac.log_analyzer.MachineLogs

.. autoclass:: pylinac.log_analyzer.MachineLog

.. autoclass:: pylinac.log_analyzer.Axis

.. autoclass:: pylinac.log_analyzer.MLC

.. autoclass:: pylinac.log_analyzer.DlogHeader

.. autoclass:: pylinac.log_analyzer.DlogAxisData

.. autoclass:: pylinac.log_analyzer.TlogHeader

.. autoclass:: pylinac.log_analyzer.TlogAxisData

.. autoclass:: pylinac.log_analyzer.SubbeamManager

.. autoclass:: pylinac.log_analyzer.Subbeam

.. autoclass:: pylinac.log_analyzer.FluenceStruct

.. autoclass:: pylinac.log_analyzer.Fluence

.. autoclass:: pylinac.log_analyzer.ActualFluence

.. autoclass:: pylinac.log_analyzer.ExpectedFluence

.. autoclass:: pylinac.log_analyzer.GammaFluence

.. autoclass:: pylinac.log_analyzer.JawStruct

.. autoclass:: pylinac.log_analyzer.CouchStruct

.. autofunction:: pylinac.log_analyzer.anonymize
