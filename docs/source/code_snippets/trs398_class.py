"""A script to calculate TRS-398 dose using pylinac classes and following the TRS-398 photon form"""
from pylinac.calibration import trs398


ENERGY = 6
TEMP = 22.1
PRESS = trs398.mmHg2kPa(755.0)
CHAMBER = '30013'  # PTW
K_ELEC = 1.000
ND_w = 5.443  # Gy/nC
MU = 200
CLINICAL_PDD = 66.5

trs398_6x = trs398.TRS398Photon(
    unit='TrueBeam1',
    setup='SSD',
    chamber=CHAMBER,
    temp=TEMP, press=PRESS,
    n_dw=ND_w,
    clinical_pdd_zref=CLINICAL_PDD,
    tpr2010=(38.2/66.6),
    energy=ENERGY,
    fff=False,
    k_elec=K_ELEC,
    voltage_reference=-300, voltage_reduced=-150,
    m_reference=(25.65, 25.66, 25.65),
    m_opposite=(25.64, 25.65, 25.65),
    m_reduced=(25.64, 25.63, 25.63),
    mu=MU, tissue_correction=1.0
)

# Done!
print(trs398_6x.dose_mu_zmax)

# examine other parameters
print(trs398_6x.kq)
print(trs398_6x.k_s)
print(trs398_6x.k_tp)

# change readings if you adjust output
trs398_6x.m_reference_adjusted = (25.44, 25.44, 25.43)
# print new dose value
print(trs398_6x.dose_mu_zmax_adjusted)

# generate a PDF for record-keeping
trs398_6x.publish_pdf('TB1 6MV TRS-398.pdf', notes=['My notes', 'I used Pylinac to do this; so easy!'], open_file=False)
