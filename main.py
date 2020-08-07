from pylinac import WinstonLutz

wl = WinstonLutz.from_zip('N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 05/2020_08_04/QA Winston-Lutz U5/RI.QA Winston-Lutz U5.MV_90_0a.zip')
#Print 3D Target Offset
#print(wl.bb_shift_instructions())
#my_directory = 'N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 05/2020_08_04/QA Winston-Lutz U5'

wl.publish_pdf('reportBuilder3.pdf')
