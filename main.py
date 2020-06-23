from pylinac import WinstonLutz
from pylinac import settings
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

my_directory ="N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2019/Unit 3 SN 3484/WINSTON LUTZ/U3 WL/6X_Final/Z20130715"

wl = WinstonLutz(my_directory)
wl.plot_images()
wl.plot_deltas()
wl.results()
#wl.bb_shift_instructions()
