
==========================
TG-51 module documentation
==========================

Overview
--------

.. automodule:: pylinac.tg51
    :no-members:

Equation Definitions
--------------------

Equation definitions are as follows:

* Ptp (Temp/Pressure correction) - TG-51 Eqn. 10:

  .. math:: \frac{273.2+T}{273.2+22} * \frac{101.33}{P}

* Ppol (Polarity correction) - TG-51 Eqn. 9:

  .. math:: |\frac{M^+_{raw}-M^-_{raw}}{2*M_{raw}}|

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

* kQ for Electrons (cylindrical chamber-specific quality conversion factor) - `Muir & Rodgers 2014`_

  The study of Muir & Rodgers was to find kecal values that could be determined soley from R50. Through
  Monte Carlo experiments, the optimal Pgradient was determined as well as fitting parameters for numerous
  common ion chambers. That study eliminates the need for Pgradient measurements. These kecal values will
  very likely be incorporated into the next TG-51 addendum (as has their kQ values for photons in the first
  addendum). From the paper, we can start with the known relationship given in Eqn. 9:

  .. math::
        k_{Q} = k_{Q,ecal} * k'_{Q}

  where Eqn. 11 states:

  .. math::
        k'_{Q} = a + b * R_{50}^{-c}

  Where a, b, and c are chamber-specific fitting values as given in Table VII.
  and where :math:`k_{Q,ecal}` is given in Table VI.

* :math:`D^Q_{w}` photon (Dose to water at 10cm from a photon beam of quality Q - TG-51 Eqn. 3:

  .. math::
       M*k_{Q}*N^{60Co}_{D,w}  (Gy)

* :math:`D^Q_{w}` electron (Dose to water at 10cm from an electron beam of quality Q - TG-51 Eqn. 6:

  .. math::
       M*P^Q_{gr}*k'_{R_{50}}*k_{ecal}*N^{60Co}_{D,w}  (Gy)

.. _Muir & Rodgers 2014: http://onlinelibrary.wiley.com/doi/10.1118/1.4893915/abstract


Typical Use
-----------

Using the TG-51 module can be complementary to your existing workflow, or completely replace it. For example,
you could use the kQ function to calculate kQ and then calculate the other corrections and values yourself.
If you want something a little more complete, you can use the :class:`~pylinac.tg51.TG51Photon` and
:class:`~pylinac.tg51.TG51Electron` classes which will completely calculate all necessary corrections and
values.

.. note::
    The Photon class uses kQ values from the TG-51 addendum and the Electron class uses kQ values from Muir & Rodgers
    2014.

Function-based Use
^^^^^^^^^^^^^^^^^^

.. literalinclude:: code_snippets/tg51_function.py


Class-based Use
^^^^^^^^^^^^^^^

.. literalinclude:: code_snippets/tg51_class.py


API Documentation
-----------------

.. autofunction:: pylinac.tg51.p_tp

.. autofunction:: pylinac.tg51.p_pol

.. autofunction:: pylinac.tg51.p_ion

.. autofunction:: pylinac.tg51.d_ref

.. autofunction:: pylinac.tg51.r_50

.. autofunction:: pylinac.tg51.kp_r50

.. autofunction:: pylinac.tg51.pq_gr

.. autofunction:: pylinac.tg51.m_corrected

.. autofunction:: pylinac.tg51.pddx

.. autofunction:: pylinac.tg51.kq

.. autoclass:: pylinac.tg51.TG51Photon

.. autoclass:: pylinac.tg51.TG51Electron
