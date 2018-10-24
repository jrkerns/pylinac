"""A script to calculate TRS-398 dose using pylinac functions and following the TRS-398 photon form"""
from pylinac.calibration import trs398


TEMP = 22.1
PRESS = trs398.mmHg2kPa(755.0)
CHAMBER = '30013'  # PTW
K_ELEC = 1.000
ND_w = 5.443  # Gy/nC
MU = 200


# Section 3 (dosimeter corrections)
k_tp = trs398.k_tp(temp=TEMP, press=PRESS)
k_pol = trs398.k_pol(m_reference=(25.66, 25.67, 25.66), m_opposite=(25.65, 25.66, 25.66))
k_s = trs398.k_s(voltage_reference=300, voltage_reduced=150,
                 m_reference=(25.66, 25.67, 25.66), m_reduced=(25.63, 25.65, 25.64))
m_corrected = trs398.m_corrected(m_reference=(25.66, 25.67, 25.66),
                                 k_tp=k_tp, k_elec=K_ELEC, k_pol=k_pol, k_s=k_s) \
              / MU

# Section 4 (kQ + dose at zref)
kq = trs398.kq_photon(chamber=CHAMBER, tpr=(39.2/68.1))
dose_mu_zref = m_corrected * ND_w * kq

# Section 5 (Dose at zmax)
# SSD setup
CLINICAL_PDD = 66.5
dose_mu_zmax = dose_mu_zref * 100 / CLINICAL_PDD

# SAD setup
CLINICAL_TMR = 0.666
dose_mu_zmax = dose_mu_zref / CLINICAL_TMR

# Done!
print(dose_mu_zmax)
