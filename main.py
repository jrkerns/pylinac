from pylinac import WinstonLutz
from pylinac import settings
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

#my_directory = "C:/Users/adnanhafeez/Documents/Python Scripts/TBCCpylinac/WL_Data/Annual/Unit 06/2018-01-01/4Mar2018 WL - corrected/Z20130715/"
my_directory = "C:/Users/adnanhafeez/Documents/Python Scripts/TBCCpylinac/WL_Data/Unit 05/2020_01_30/QA Winston-Lutz U5 2019"
#my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2019/Unit 5 SN 2899/WINSTON LUTZ/6X_Final/Z20130715"

#my_directory = "N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 01/2020_01_27/U1/WL/Z20190814"

wl = WinstonLutz(my_directory)
wl.plot_deltas()
