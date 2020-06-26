from pylinac import WinstonLutz


#my_directory ="N:/MedicalPhysics/Quality Control Program/Treatment/Annual/QC/2018/Unit 8 SN 1495/Winston-Lutz/Z20130715"
#my_directory ="N:/MedicalPhysics/Quality Control Program/Treatment/Weekly/Winston-Lutz/Unit 03/2020_06_23 redo/QA Winston-Lutz U3"

my_directory ="H:/WL Testing on U2/Forth Run/Z20190312/"

wl = WinstonLutz(my_directory)
wl.plot_images()
print(wl.plot_deltas())
#print(wl.results()[3][-7:-2])
#print(wl.bb_shift_instructions())
print(wl.results())