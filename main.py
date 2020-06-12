from pylinac import WinstonLutz
from pylinac import settings
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2019/Unit 3 SN 3484/WINSTON LUTZ/U3 WL/6X_Final/Z20130715"
#my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2019/Unit 6 SN 1245/WINSTON LUTZ/U6 WL/Full WL/Z20190312"
#my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2019/Unit 5 SN 2899/WINSTON LUTZ/6X_Final/Z20130715"

#my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 01/2020_01_27/U1/WL/Z20190814"

wl = WinstonLutz(my_directory)
print(wl.bb_shift_instructions())

wl.plot_images()
#wl.save_images('test.png')
#plt.savefig('test.png')#

wl.publish_pdf('unit1.pdf')
