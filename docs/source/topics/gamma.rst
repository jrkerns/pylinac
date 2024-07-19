.. _gamma:

Gamma
=====

Overview
--------

Gamma the famous dose distribution calculation tool introduced by Low et al [1]_ [2]_.

Within pylinac, gamma can be calculated for 1D dose distributions reliably.
While there is a gamma function for 2d gamma, it is in the stages of being refactored for
easier use; it is thus not get documented fully here yet.

.. note::

  This discussion is for 1D gamma only.

As it relates to pylinac, the gamma is focused on comparing QA data measurements.

.. danger::

    The gamma functions in pylinac are not meant to be used for patient dose distributions.
    They are meant for machine QA data comparisons.

The gamma calculation in pylinac is also not fast. There are faster algorithms and libraries out there; e.g.
``pymedphys.gamma`` `here <https://docs.pymedphys.com/en/latest/users/howto/gamma/index.html>`__.

Assumptions
-----------

There is discussion in the literature about the assumptions of gamma and also assumptions
about the types of data to be evaluated. The following are a list of assumptions in pylinac:

#. The reference distribution is the ground truth distribution.
#. The reference distribution has more or equal resolution to the evaluated distribution.
#. We want to compute the gamma over each evaluation point within the context of a more
   dense reference distribution.
   E.g. comparing an IC Profiler to a finely-sampled water tank scan where the tank data is the commissioning data.
#. The evaluation distribution is within the physical bounds of the reference distribution.
   E.g. comparing a 10x10 field to a 20x20 reference field.

Usage
-----

We will use the geometric 1mm, DD example from `Agnew & McGarry <https://www.sciencedirect.com/science/article/abs/pii/S0167814015006660>`__.
Although their data is 2D, we can extract a 1D profile from their data.

.. plot::

  import matplotlib.pyplot as plt

  from pylinac.core.io import retrieve_demo_file
  from pylinac.core.gamma import gamma_1d
  from pylinac.core.image import DicomImage

  # get the files from the cloud
  ref_file = retrieve_demo_file(name='gamma/HN_VMAT_Reference_1mmPx.dcm')
  eval_file = retrieve_demo_file(name='gamma/HN_VMAT_Evaluated_1mmPx.dcm')
  # load them as DICOM images; we want to use the raw pixels
  ref_img = DicomImage(ref_file, raw_pixels=True)
  eval_img = DicomImage(eval_file, raw_pixels=True)
  # convert the images to float arrays; uints can result in overflow errors
  ref_array = ref_img.array.astype(float)
  eval_array = eval_img.array.astype(float)
  # take a random sample profile through the middle of the images
  ref_prof = ref_img[:, 90]
  eval_prof = eval_img[:, 90]

  # compute the gamma
  gamma_map, ref_vals, ref_x_vals = gamma_1d(ref_prof, eval_prof, resolution_factor=7, dose_threshold_percent=0)

  # plot the results
  fig, ax = plt.subplots(figsize=(10, 7))
  ax.plot(ref_prof, 'k+', label='Original Reference', )
  ax.plot(ref_x_vals, ref_vals, 'cx', label='Reference Interpolated',)
  ax.plot(eval_prof, 'bo', label='Evaluation')
  ax.set_ylabel('Pixel Value')
  ax.set_xlabel('Pixel Position')
  ax.legend(loc='upper left')
  g_ax = ax.twinx()
  g_ax.plot(gamma_map, 'rx', label='Gamma')
  g_ax.legend(loc='upper right')
  g_ax.set_ylabel('Gamma')
  fig.suptitle(f'1D Gamma Analysis; max \N{Greek Small Letter Gamma}:{gamma_map.max():.2f}; avg \N{Greek Small Letter Gamma}: {gamma_map.mean():.2f}; pass rate: {100*gamma_map[gamma_map <= 1].size/gamma_map.size:.2f}%')
  plt.show()

Explanation
------------

Gamma function
^^^^^^^^^^^^^^

The gamma function of Table 1 in Low et al [2]_ is used. The generalized gamma function is defined as:

.. math::

    \Gamma(\vec{r}_{e},\vec{r}_{r}) = \sqrt{ \frac{r^2(\vec{r}_{e},\vec{r}_{r})}{\Delta d^{2}} + \frac{ \delta(\vec{r}_{e},\vec{r}_{r}) }{\Delta D^2} }

computed for each evaluation point and reference point.

The final gamma value for a given evaluation point is the minimum gamma value of the above function
for all reference points within the search radius:

.. math::

    \gamma(\vec{r}_{e}) = \min \{ \Gamma(\vec{r}_{e},\vec{r}_{r}) \} \forall \{ \vec{r}_{r} \}

.. important::

    Per assumption #3, we will perform the gamma on each **evaluation** point, not the reference point
    based on the other assumptions which is the inverse of Low.

For any residual definitions, see Table I of Low et al [2]_.

Logic
^^^^^

.. note::

    The actual gamma calculation is relatively straightforward. However, the implementation
    when the reference and evaluation data are not the same resolution or at the same physical
    locations are where the details matter.

#. The threshold for computing gamma is set as :math:`\max \{ \text{reference} \} \cdot \text{threshold parameter} \cdot 100` . E.g. 10% is a common value under which the gamma is not calculated.
#. If using global dose: The dose-to-agreement :math:`\Delta D` is computed as: :math:`\max \{ \text{reference} \} \cdot \text{DTA parameter} \cdot 100`. E.g. 3% is a common value.
#. If using local dose: The dose-to-agreement :math:`\Delta D` is computed as: :math:`\vec{r}_{r} \cdot \text{DTA parameter} \cdot 100`.
#. The reference distribution is mapped to an linear interpolator. This is so that the reference data can be evaluated at any desired point, not just the discrete points passed.
#. For each evaluation point :math:`D_{e}(\vec{r}_{e})`:

   #. The reference points to evaluate over are extracted. This will be from :math:`-\Delta d` to :math:`+\Delta d` from the evaluation point.
      The number of reference points depend on the ``resolution_factor`` parameter.
      E.g. a factor of 3 will result in 7 reference points sitting about the evaluation point. E.g. (-3, -2, -1, 0, 1, 2, 3).
   #. For each reference point in the array above above :math:`D_{r}(\vec{r}_{r})`, the gamma function is computed :math:`\Gamma(\vec{r}_{e},\vec{r}_{r})`.
   #. The minimum gamma is taken as the gamma value :math:`\gamma(\vec{r}_{e})` for that evaluation point. If the minimum gamma
      is larger than ``gamma_cap_value``, the ``gamma_cap_value`` is used instead.

The resulting gamma array will be the same size as the evaluation distribution. I.e. gamma is calculated at each evaluation point.


.. [1] Low, Harms, Mutic, Purdy. A technique for the quantitative evaluation of dose distributions. Med Phys. 1998;25(5):656-61.

.. [2] Low, Dempsey. Evaluation of the gamma dose distribution comparison method. Med Phys. 2003;30(9):2455-64.

API
---

.. automodule:: pylinac.core.gamma
   :members:
