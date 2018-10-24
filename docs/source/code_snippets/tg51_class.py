"""A script to calculate TG-51 dose using pylinac classes and following the TG-51 photon form"""
from pylinac.calibration import tg51


ENERGY = 6
TEMP = 22.1
PRESS = tg51.mmHg2kPa(755.0)
CHAMBER = '30013'  # PTW
P_ELEC = 1.000
ND_w = 5.443  # Gy/nC
MU = 200
CLINICAL_PDD = 66.5

tg51_6x = tg51.TG51Photon(
    unit='TrueBeam1',
    chamber=CHAMBER,
    temp=TEMP, press=PRESS,
    n_dw=ND_w, p_elec=P_ELEC,
    measured_pdd10=66.4, lead_foil=None,
    clinical_pdd10=66.5, energy=ENERGY,
    voltage_reference=-300, voltage_reduced=-150,
    m_reference=(25.65, 25.66, 25.65),
    m_opposite=(25.64, 25.65, 25.65),
    m_reduced=(25.64, 25.63, 25.63),
    mu=MU, tissue_correction=1.0
)

# Done!
print(tg51_6x.dose_mu_dmax)

# examine other parameters
print(tg51_6x.pddx)
print(tg51_6x.kq)
print(tg51_6x.p_ion)

# change readings if you adjust output
tg51_6x.m_reference_adjusted = (25.44, 25.44, 25.43)
# print new dose value
print(tg51_6x.dose_mu_dmax_adjusted)

# generate a PDF for record-keeping
tg51_6x.publish_pdf('TB1 6MV TG-51.pdf', notes=['My notes', 'I used Pylinac to do this; so easy!'], open_file=False)
