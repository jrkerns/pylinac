from pylinac import WinstonLutz


my_directory = 'N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 05/2020_06_30/QA Winston-Lutz U5/'

wl = WinstonLutz(my_directory)

#wl.publish_pdf('reportBuilder.pdf', my_directory)
for k, v in (wl.plot_deltas()).items():
    print("{}: {}mm".format(k, v))
