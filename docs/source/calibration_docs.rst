.. _calibration_module:

================================
Calibration module documentation
================================

Overview
--------

The calibration module actually consists of two submodules: ``tg51`` and ``trs398``, each addressing their respective protocol.
Both modules contain functions and classes for calculation the protocol dose. The modules have some overlap, especially with
basic functions as well as helper functions. The modules have tried to use the vocabulary of the respective protocol, but
occasionally there are differences when we felt that using the same name was clearer. See the vocabulary section for full definitions.

.. note::

    Besides the typical calculations one would expect, the modules also include helper functions, such as a PDD to TPR
    converter so that a TG-51 calculation can determine kQ from TPR and avoid the tedious PDDx. Additionally, pressure
    unit converters exist to go from the various units of pressure to kPa which is what pylinac uses.


Vocabulary that may be different than the protocol
--------------------------------------------------

* ``voltage_reference``: Used in both TG-51 and TRS-398 for the voltage used when taking a reference reading; commonly -300V.
* ``voltage_reduced``: Use in both TG-51 and TRS-398 for the lower voltage used to determine ks/Pion; commonly -150V.
* ``m_reference``: The ion chamber reading at the reference voltage.
* ``m_reduced``: The ion chamber reading at the reduced voltage.
* ``m_opposite``: The ion chamber reading at the opposite polarity of the reference voltage: commonly +300V.

Vocabulary not listed here should be the same as the respective protocol.

TG-51
-----


Equation Definitions
^^^^^^^^^^^^^^^^^^^^

Equation definitions are as follows:

* Ptp (Temp/Pressure correction) - TG-51 Eqn. 10:

  .. math:: \frac{273.2+T}{273.2+22} * \frac{101.33}{P}

  .. warning::

        Temperature is in Celsius and pressure is in kPa. Use the helper functions
        :func:`~pylinac.calibration.tg51.fahrenheit2celsius`,
        :func:`~pylinac.calibration.tg51.mmHg2kPa`, and
        :func:`~pylinac.calibration.tg51.mbar2kPa` as needed.

* Ppol (Polarity correction) - Rather than using TG-51 Eqn. 9, we opt instead for TRS-398 Eqn xx, which uses the absolute values
  of the positive and negative voltages. This is the same results as Eqn. 9 but without worrying about signs.:

  .. math:: \frac{|M^+_{raw}|+|M^-_{raw}|}{2*M_{raw}}

* Pion (Ion collection correction; only for pulsed beams) - TG-51 Eqn. 12:

  .. math:: \frac{1-\frac{V_{H}}{V_{L}}}{\frac{M^H_{raw}}{M^L_{raw}} - \frac{V_{H}}{V_{L}}}

* Dref (Reference electron depth; cm) - TG-51 Eqn. 18:

  .. math:: 0.6*R_{50} - 0.1

* R50 (Beam quality specifier; 50% dose depth; cm) - TG-51 Eqn. 16 & 17:

  .. math::
       \begin{cases}
          1.029*I_{50}-0.06 (cm) & 2\leq I_{50}\leq 10 \\
          1.059*I_{50}-0.37 (cm) & I_{50}\gt 10 \\
       \end{cases}

* k'R50 (k'R50 for cylindrical chambers) - TG-51 Eqn. 19:

  .. math::
        0.9905+0.0710e^\frac{-R_{50}}{3.67}

* PQ_gr (PQ gradient correction for cylindrical chambers) - TG-51 Eqn. 21:

  .. math::
        \frac{M_{raw}(d_{ref}+0.5*r_{cav})}{M_{raw}*d_{ref}}

* PDDx (PDD from photons only) - TG-51 Eqns. 13, 14 & 15:

  .. math::
      \begin{cases}
          PDD(10) & energy < 10 \\
          1.267*PDD(10)-20.0 & 75\leq PDD(10)\leq 89 \\
          PDD(10)_{Pb} & lead@50cm, PDD(10)_{Pb} < 73 \\
          (0.8905+0.00150*PDD(10)_{Pb})*PDD(10)_{Pb} & lead@50cm, PDD(10)_{Pb} \geq 73 \\
          PDD(10)_{Pb} & lead@30cm, PDD(10)_{Pb} < 71 \\
          (0.8116+0.00264*PDD(10)_{Pb})*PDD(10)_{Pb} & lead@30cm, PDD(10)_{Pb} \geq 71
      \end{cases}

* M-corrected (corrected chamber reading) - TG-51 Eqn. 8:

  .. math::
       P_{ion}*P_{TP}*P_{elec}*P_{pol}*M_{raw}

* kQ for Photons (cylindrical chamber-specific quality conversion factor) - TG-51 Addendum Eqn 1 & Table I:

  .. math::
        \begin{cases}
            A + B*10^-3*PDD(10)x+C*10^-5*(PDD(10)x)^2 & 63 < PDD(10)x < 86  \\
        \end{cases}

  Where A, B, and C are chamber-specific fitting values as given in Table I.
  Pylinac automatically retrieves values based on the chamber model passed to the function.

* kQ for Electrons (cylindrical chamber-specific quality conversion factor) - `Muir & Rogers 2014`_

  The study of Muir & Rogers was to find kecal values that could be determined soley from R50. Through
  Monte Carlo experiments, the optimal Pgradient was determined as well as fitting parameters for numerous
  common ion chambers. That study eliminates the need for Pgradient measurements. These kecal values will
  very likely be incorporated into the next TG-51 addendum (as has their kQ values for photons in the first
  addendum). From the paper, we can start with the known relationship given in Eqn. 9:

  .. math::
        k_{Q} = k_{Q,ecal} * k'_{Q}

  where Eqn. 11 states:

  .. math::
        k'_{Q} = a + b * R_{50}^{-c}

  Where a, b, and c are chamber-specific fitting values as given in Table VII
  and where :math:`k_{Q,ecal}` is given in Table VI.

* :math:`D^Q_{w}` photon (Dose to water at 10cm from a photon beam of quality Q - TG-51 Eqn. 3:

  .. math::
       M*k_{Q}*N^{60Co}_{D,w}  (Gy)

* :math:`D^Q_{w}` electron (Dose to water at 10cm from an electron beam of quality Q - TG-51 Eqn. 6:

  .. math::
       M*P^Q_{gr}*k'_{R_{50}}*k_{ecal}*N^{60Co}_{D,w}  (Gy)

.. _Muir & Rogers 2014: http://onlinelibrary.wiley.com/doi/10.1118/1.4893915/abstract


Function-based Use
^^^^^^^^^^^^^^^^^^

Using the TG-51 module can be complementary to your existing workflow, or completely replace it. For example,
you could use the kQ function to calculate kQ and then calculate the other corrections and values yourself.
If you want something a little more complete, you can use the :class:`~pylinac.tg51.TG51Photon`,
:class:`~pylinac.tg51.TG51ElectronLegacy` and :class:`~pylinac.tg51.TG51ElectronModern` classes which will calculate all necessary corrections and
values.

.. note::
    The Photon class uses kQ values from the TG-51 addendum.
    The Legacy Electron class will make the user specify a kecal value and measure Pgradient.
    The Modern Electron class will calculate kQ completely from R50 and the chamber from Muir & Rogers 2014 paper,
    no kecal or Pgradient needed.


.. literalinclude:: code_snippets/tg51_function.py


Class-based Use
^^^^^^^^^^^^^^^

.. literalinclude:: code_snippets/tg51_class.py

.. _trs398:

TRS-398
-------

.. warning:: Pylinac does not calculate electron dose in any other conditions than water; i.e. no solid water.

Equation Definitions
^^^^^^^^^^^^^^^^^^^^

* Ktp (Temp/Pressure correction):

  .. math:: \frac{273.2+T}{273.2+22} * \frac{101.33}{P}

  .. warning:: Temperature is in Celsius and pressure is in kPa. Use the helper functions fahrenheit2celsius, mmHg2kPa, and mbar2kPa as needed.

* Kpol (Polarity correction):

  .. math:: \frac{|M^+_{raw}|+|M^-_{raw}|}{2*M_{raw}}

* Ks (Ion collection correction; only for pulsed beams):

  .. math:: a_{0} + a_{1}*(\frac{M_{1}}{M_{2}}) + a_{2}*(\frac{M_{1}}{M_{2}})^2

* Zref (Reference electron depth; cm) - TRS-398 7.2:

  .. math:: 0.6*R_{50} - 0.1

* R50 (Beam quality specifier; 50% dose depth; cm) - TRS-398 7.1:

  .. math::
       \begin{cases}
          1.029*I_{50}-0.06 (cm) & 2\leq I_{50}\leq 10 \\
          1.059*I_{50}-0.37 (cm) & I_{50}\gt 10 \\
       \end{cases}

* :math:`D^Q_{w}` photon (Dose to water at Zref from a photon or electron beam of quality Q - TRS-398 7.3:

  .. math::
       D_{w,Q} = M_{Q}*N_{D,w,Qo}*k_{Q,Qo}  (Gy)

* M-corrected (corrected chamber reading):

  .. math::
       M_{Q} = k_{s}*k_{TP}*K_{elec}*K_{pol}*M_{1}

* kQ,Qo for Photons (cylindrical chamber-specific quality conversion factor): TRS-398 Table 6.III

* kQ for Electrons (cylindrical chamber-specific quality conversion factor; calibrated in Co-60): TRS-398 Table 7.III


Function-based Use
^^^^^^^^^^^^^^^^^^

.. literalinclude:: code_snippets/trs398_function.py


Class-based Use
^^^^^^^^^^^^^^^

.. literalinclude:: code_snippets/trs398_class.py


TG-51 API Documentation
-----------------------

.. autofunction:: pylinac.calibration.tg51.mmHg2kPa

.. autofunction:: pylinac.calibration.tg51.mbar2kPa

.. autofunction:: pylinac.calibration.tg51.fahrenheit2celsius

.. autofunction:: pylinac.calibration.tg51.tpr2010_from_pdd2010

.. autofunction:: pylinac.calibration.tg51.p_tp

.. autofunction:: pylinac.calibration.tg51.p_pol

.. autofunction:: pylinac.calibration.tg51.p_ion

.. autofunction:: pylinac.calibration.tg51.d_ref

.. autofunction:: pylinac.calibration.tg51.r_50

.. autofunction:: pylinac.calibration.tg51.kp_r50

.. autofunction:: pylinac.calibration.tg51.pq_gr

.. autofunction:: pylinac.calibration.tg51.m_corrected

.. autofunction:: pylinac.calibration.tg51.pddx

.. autofunction:: pylinac.calibration.tg51.kq_photon_pddx

.. autofunction:: pylinac.calibration.tg51.kq_photon_tpr

.. autofunction:: pylinac.calibration.tg51.kq_electron

.. autoclass:: pylinac.calibration.tg51.TG51Photon

.. autoclass:: pylinac.calibration.tg51.TG51ElectronLegacy

.. autoclass:: pylinac.calibration.tg51.TG51ElectronModern

TRS-398 API Documentation
-------------------------

.. autofunction:: pylinac.calibration.trs398.k_s

.. autofunction:: pylinac.calibration.trs398.kq_photon

.. autofunction:: pylinac.calibration.trs398.kq_electron

.. autofunction:: pylinac.calibration.trs398.m_corrected

.. autoclass:: pylinac.calibration.trs398.TRS398Photon

.. autoclass:: pylinac.calibration.trs398.TRS398Electron

