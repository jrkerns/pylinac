from pylinac import WinstonLutz

wl = WinstonLutz.from_zip('C:/Users/adnanhafeez/Documents/Python Scripts/TBCCpylinac/WL_Data/Unit 01/2019_08_12/Z20150323/Z20150323.zip')

for k, v in (wl.plot_deltas()).items():
    print("{}: {}mm".format(k, v))



