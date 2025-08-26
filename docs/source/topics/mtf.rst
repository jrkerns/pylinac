
.. _mtf_topic:

Modulation Transfer Function
============================

The **Modulation Transfer Function (MTF)** is a standard way of describing how well
an imaging system preserves detail and contrast. It measures the ratio of image contrast
to object contrast at different levels of detail, expressed in terms of spatial frequency.

In simpler terms, MTF shows how clearly an imaging system can reproduce fine patterns:
large features are usually transferred well, while very small details tend to blur and lose contrast.

.. _peak-valley-mtf:

Peak-Valley MTF
---------------

The Peak-Valley MTF is used in CBCT and planar imaging metrics to describe high-contrast characteristics of the imaging system.
An excellent introduction is `here <https://www.edmundoptics.com/knowledge-center/application-notes/optics/introduction-to-modulation-transfer-function/>`__.
In pylinac, Peak-Valley MTF is calculated using equation 3 of the above reference, which is also the :ref:`Michelson <michelson>` contrast definition.

.. math:: contrast = \frac{I_{max} - I_{min}}{I_{max} + I_{min}}

Then, all the contrasts are normalized to the largest one, resulting in a normalized MTF or rMTF (relative).
Pylinac only reports rMTF values. This is the first of two inputs. The other is the line pair spacing. The spacing
is usually provided by the phantom manufacturer. The rMTF is the plotted against the line pair/mm values. Also from
this data the MTF at a certain percentage (e.g. 50%) can be determined in units of lp/mm.

However, it's important to know what :math:`I_{max}` and :math:`I_{min}` means here. For a line pair set, each bar and space-between
is one contrast value. Thus, one contrast value is calculated for each bar/space combo. For phantoms with areas of the
same spacing (e.g. the Leeds), all bars and spaces are the same and thus we can use an area-based ROI for the input to
the contrast equation.

.. _ESF_MTF:

ESF-based MTF
-------------

Logic
~~~~~

The Modulation Transfer Function (MTF) can be determined from an **edge spread function** (ESF)
by first calculating the **line spread function** (LSF) from the ESF,
then taking the Fourier transform of the LSF. The magnitude of the resulting complex function is the MTF.

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.ndimage import gaussian_filter1d
    from scipy.fft import fft, fftfreq

    n = 100
    n_pad = 1024
    x = np.arange(0,n)
    y = np.zeros(n)
    y[n//2:] = 1
    esf = gaussian_filter1d(y, sigma=10)
    lsf = np.gradient(esf)
    lsf_padded = np.pad(lsf, n_pad, mode='constant', constant_values=0)
    mtf = np.abs(fft(lsf_padded))
    mtf /= mtf[0]  # Normalize
    freq = fftfreq(len(lsf_padded))
    half = len(freq) // 2
    freq = freq[:half]
    mtf = mtf[:half]

    fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(6, 2.5))

    ax = axs[0]
    ax.plot(x,esf)
    ax.set_title("ESF")
    ax.set_xlabel("Pixel #")
    ax.set_ylabel("Pixel count")
    ax.set_yticks([])

    ax = axs[1]
    ax.plot(x,lsf)
    ax.set_title("LSF")
    ax.set_xlabel("Pixel #")
    ax.set_ylabel('$\mathit{dESF/dx}$')
    ax.set_yticks([])

    ax = axs[2]
    ax.plot(freq,mtf)
    ax.set_title("MTF")
    ax.set_xlabel("Frequency [cycles/pixel]")
    ax.set_ylabel(r"$\mathcal{F}[LSF]$")
    ax.set_xlim(0,0.07)

    plt.show()

Implementation
~~~~~~~~~~~~~~

In **pylinac**, the MTF is computed from a **list of profiles (edge spread functions)**
for example, multiple profiles taken at different edge positions of a phantom.
The profiles can have different lengths (by default they are padded to a common size).
Since it is not guaranteed that all profiles are aligned, the MTF is calculated for each profile individually,
and then the results are averaged to obtain a single MTF.

It is assumed that the **sample spacing** (input parameter in mm) is the same for all profiles.

* If ``sample_spacing`` is None (default) the frequency axis is expressed in cycles per sample, up to 0.5 (Nyquist frequency).
* If ``sample_spacing`` is provided, the frequency axis is expressed in line pairs per mm,
  up to the Nyquist limit determined by the sample spacing.

.. code-block::

    # Default: frequency axis is expressd in cycles/sample
    EdgeSpreadFunctionMTF(...)
    # Same as above
    EdgeSpreadFunctionMTF(..., sample_spacing = None)
    # Sample_spacing in mm: frequency axis is expressd in line pairs/mm
    EdgeSpreadFunctionMTF(..., sample_spacing = 0.5)

Each profile MTF is computed using the following sequence of operations:

1. **Differentiate ESF → LSF**

   Compute the **line spread function (LSF)** as the discrete derivative of the ESF,
   implemented with a central-difference method (``lsf = np.gradient(esf)``).

2. **Apply windowing**

   Multiply the LSF by a windowing function to taper the edges and reduce spectral
   leakage introduced by the finite measurement length. The default window is
   ``scipy.signal.windows.hann``. This can be modified using the ``windowing`` parameter.
   Input arguments to the windowing function can also be passed in the constructor through ``kwargs``.

   (see more windowing functions here: https://docs.scipy.org/doc/scipy/reference/signal.windows.html)

   .. code-block::

    # Default window: Hann
    EdgeSpreadFunctionMTF(...)
    # No windowing
    EdgeSpreadFunctionMTF(..., windowing=None)
    # Custom window with default parameters
    EdgeSpreadFunctionMTF(..., windowing=scipy.signal.windows.tukey)
    # Custom window with custom parameters
    EdgeSpreadFunctionMTF(..., windowing=scipy.signal.windows.tukey, alpha=0.3)

3. **Zero-pad for frequency resolution**

   Pad the LSF to the target length (default: ``1024`` samples). Padding does not
   add information, but interpolates the frequency spectrum more finely, improving the
   resolution of the displayed MTF curve.
   This can be modified using the parameters ``padding_mode`` and ``num_samples``:

   .. code-block::

    # no padding
    EdgeSpreadFunctionMTF(..., padding_mode="none")
    # array zero-padded to 2048 elements (must be larger than the largest array)
    EdgeSpreadFunctionMTF(..., padding_mode="fixed", num_samples=2048)
    # array padded to the next power of 2 or num_samples
    (i.e. len(largest_array) == 200 => num_samples=1024, len(largest_array) == 1025 => num_samples=2048)
    EdgeSpreadFunctionMTF(..., padding_mode="auto", num_samples=1024)

4. **Fourier transform**

   Compute the fast Fourier transform (FFT) of the windowed and padded LSF.

5. **Magnitude spectrum → MTF**

   Take the modulus (magnitude) of the FFT result to obtain the (unnormalized) MTF.

6. **Normalization**

   Normalize the MTF by dividing by the first value so that ``MTF(0) = 1``.


Moments-based MTF
-----------------

The MTF can also be calculated using the moments of the line pair spread function (LPSF).
This algorithm is based on the work of Hander et al [1]_. Specifically, equation 8:

.. math::

   MTF = \frac{\sqrt{2 * (\sigma^{2} - \mu)}}{\mu}

where :math:`\mu` is the mean pixel value of the ROI and :math:`\sigma` is the standard deviation of the ROI pixel values.


.. [1] `Hander et al. <https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.597928>`__ "Rapid objective measurement of gamma camera resolution using statistical moments" (1998).


API Documentation
-----------------

.. automodule:: pylinac.core.mtf
