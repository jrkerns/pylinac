"""A script to calculate TG-51 dose using pylinac classes and following the TG-51 photon form"""
from pylinac import tg51


ENERGY = 6
TEMP = 22.1
PRESS = 755.0
CHAMBER = '30013'  # PTW
P_ELEC = 1.000
ND_w = 5.443  # Gy/nC
MU = 200
CLINICAL_PDD = 66.5

tg51_6x = tg51.TG51Photon(temp=TEMP, press=PRESS, model=CHAMBER,
                          n_dw=ND_w, p_elec=P_ELEC,
                          measured_pdd=66.4, lead_foil=None,
                          clinical_pdd=66.5, energy=ENERGY,
                          volt_high=-300, volt_low=-150,
                          m_raw=(25.65, 25.66, 25.65),
                          m_opp=(25.64, 25.65, 25.65),
                          m_low=(25.64, 25.63, 25.63),
                          mu=MU, tissue_correction=1.0)

# Done!
print(tg51_6x.dose_mu_dmax)

# examine other parameters
tg51_6x.pddx
tg51_6x.kq
tg51_6x.p_ion

# change readings if you adjust output
tg51_6x.m_raw = (25.44, 25.44, 25.43)
# print new dose value
print(tg51_6x.dose_mu_dmax)

