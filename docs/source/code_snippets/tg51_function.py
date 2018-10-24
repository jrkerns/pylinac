"""A script to calculate TG-51 dose using pylinac functions and following the TG-51 photon form"""
from pylinac.calibration import tg51


ENERGY = 6
TEMP = 22.1
PRESS = tg51.mmHg2kPa(755.0)
CHAMBER = '30013'  # PTW
P_ELEC = 1.000
ND_w = 5.443  # Gy/nC
MU = 200
CLINICAL_PDD = 66.5

# Section 4 (beam quality)
# since energy is 6MV, PDDx == PDD, but we'll run it through anyway just for show
pdd10x = tg51.pddx(pdd=66.4, energy=ENERGY)

# Section 5 (kQ)
kq = tg51.kq_photon_pddx(chamber=CHAMBER, pddx=pdd10x)
# Alternatively, get kQ from TPR (way quicker to measure, without needing to measure TPR!)
tpr = tg51.tpr2010_from_pdd2010(pdd2010=(38.0/66.4))
kq = tg51.kq_photon_tpr(chamber=CHAMBER, tpr=tpr)

# Section 6 (Temp/Press)
p_tp = tg51.p_tp(temp=TEMP, press=PRESS)

# Section 7 (polarity)
m_reference = (25.66, 25.67, 25.66)
m_opposite = (25.67, 25.67, 25.68)
p_pol = tg51.p_pol(m_reference=m_reference, m_opposite=m_opposite)

# Section 8 (ionization)
m_reduced = (25.61, 25.62)
p_ion = tg51.p_ion(voltage_reference=300, voltage_reduced=150, m_reference=m_reference, m_reduced=m_reduced)

# Section 9 (M corrected)
m_corr = tg51.m_corrected(p_ion=p_ion, p_tp=p_tp, p_elec=P_ELEC, p_pol=p_pol, m_reference=m_reference)

# Section 10 (dose to water @ 10cm)
dose_10 = m_corr*kq*ND_w
dose_10_per_mu = dose_10 / MU

# Section 11 (dose/MU to water @ dmax)
dose_ddmax = dose_10_per_mu / CLINICAL_PDD

# Done!
print(dose_ddmax)
