import json
import tempfile

import numpy as np
import matplotlib.pyplot as plt
import os.path as osp
# import requests
from blackfire import probe
import plotly.express as px


import pylinac
import pydicom
import base64
import os
from pylinac.core import image_generator
from pylinac.core.image_generator import layers, simulators, GaussianFilterLayer
from pylinac.core.geometry import cos
from pylinac.core import image
from blackfire import probe

import pylinac
import pydicom
import base64
import os
from pylinac.core import image_generator
from pylinac.core.image_generator import layers, simulators, GaussianFilterLayer, generate_winstonlutz
from pylinac.core.geometry import cos
from pylinac.core import image

sim = image_generator.simulators.AS1000Image()
field_layer = image_generator.layers.FilteredFieldLayer  # could also do FilterFreeLayer
# create a set of WL images
# this will create 4 images (via image_axes len) with an offset of 3mm to the left
# the function is smart enough to correct for the offset w/r/t gantry angle.
# generate_winstonlutz(simulator=sim, field_layer=field_layer,
#                      final_layers=[GaussianFilterLayer()], gantry_tilt=0,
#                      dir_out='./wl_dir', offset_mm_left=0,
#                      image_axes=[[0, 0, 0], [180, 0, 0], [90, 0, 0], [270, 0, 0]])
#
# wl = pylinac.WinstonLutz('./wl_dir')
# print(wl.results())
# wl.publish_pdf("Perfect.pdf")
# config = dict({'scrollZoom': True})
# leeds = pylinac.LeedsTOR.from_demo_image()
# leeds.analyze()
# wl_mod = pylinac.winston_lutz
pylinac.winston_lutz.MAX_BB_SIZE = 3
wl = pylinac.WinstonLutz.from_demo_images()
wl.plot_summary()


# plt.plot(range(len(leeds.mtf.norm_mtfs)), list(leeds.mtf.norm_mtfs.values()))
# plt.show()



path = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac manuscript\VMAT\2 DRGS\Trajectory Logs\RA QA_T2_DR_GS_HDTB_T2_DR_GS_20160213173812.bin"
# im = image.load(path)
# im.plot()
# log = pylinac.TrajectoryLog(path)
l = pylinac.PicketFence.from_demo_image()
print(l.image.dpmm)
l.image.plot()
data = json.dumps(l.image.array.tolist())
ttt = 1
with open('pf_data.json', 'w') as f:
    f.write(data)
# plt.plot(im[400, 500])

# probe.initialize()
# with probe.run():
#     for _ in range(5):
#         # ds = image.DicomImage(path)
#         stack = image.DicomImageStack.from_zip(path)
        # ct = pylinac.CatPhan600.from_demo_images()
# plt.plot(im[400, 500])
import gc
total = 3
images = 100
probe.initialize()
with probe.run():
    for i in range(total):
        # items = list()
        #
        # for r in range(images):
        #     items.append(pydicom.dcmread(r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac\pylinac\demo_files\Doselab_kV.dcm"))
        # del items
        # ds = image.DicomImage(path)
        def f():
            stack = image.DicomImageStack.from_zip(path)
            return stack[0]
        f()
        # del stack
        # ct = pylinac.CatPhan600.from_demo_images()
        # del ct
        gc.collect()
        print(f"loop {i+1} of {total}")
print("Done")
    # ct.analyze()
    # ct = pylinac.CatPhan600.from_demo_images()
    # ct.analyze()
    # ct = pylinac.CatPhan600.from_demo_images()
    # ct.analyze()
    # ct = pylinac.CatPhan600.from_demo_images()
    # ct.analyze()
# ct = pylinac.CatPhan600.from_demo_images()
ct.analyze()
# ct.plot_analyzed_subimage('lc')
ct.plot_analyzed_image()
# fs = pylinac.FlatSym.from_demo_image()
# fs.analyze('varian', 'varian', invert=True)
# fs.plot_analyzed_image()

# ds = pydicom.dcmread(r'E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac\playground_dir\WL G=0, C=0, P=0; BB=5mm @ left=3, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm')
# plt.imshow(ds.pixel_array)

#%%
# from pylinac.core import image_generator
# import pydicom
# import matplotlib.pyplot as plt
sim = image_generator.simulators.AS1000Image()
field_layer = image_generator.layers.FilteredFieldLayer
cone_layer = image_generator.layers.FilterFreeConeLayer
# bb_layer = image_generator.layers.PerfectBBLayer
fs = image_generator.generate_winstonlutz_cone(simulator=sim, cone_layer=cone_layer, final_layers=[GaussianFilterLayer(sigma_mm=1)], gantry_tilt=0, dir_out='./playground_dir',
                                               cone_size_mm=15,
                                               offset_mm_left=0,
                                               offset_mm_in=0,
                                               offset_mm_up=2.5,
                                               image_axes=[[0, 0, 0], [90, 0, 0], [180, 0, 0], [270, 0, 0]])
wl = pylinac.WinstonLutz(directory='./playground_dir')
print(wl.results())
wl.plot_images()
# wl.plot_summary()
for file in os.scandir('./playground_dir'):
    os.remove(file.path)
print("Files deleted")
# for f in fs:
#     ds = pydicom.dcmread(f)
#     plt.imshow(ds.pixel_array, origin='lower')
#     plt.show()

#%%
# import pylinac
# tlog = pylinac.TrajectoryLog(r"E:\Radformation Repos\RadLib\RadLib\Rad.Dynalog.Tests\data\Tlog2-1.bin")
# f = tlog.fluence.actual.calc_map()
# tlog.fluence.actual.plot_map()

#%%
import pylinac


#%%
# import pylinac
# ct = pylinac.CatPhan504.from_demo_images()
# ct.analyze()
# ct.publish_pdf('stufpid.pdf')

#%%

import pylinac
# from pylinac.vmat import Segment
# Segment._nominal_width_mm = 15
# Segment._nominal_height_mm = 150
vm = pylinac.DRGS.from_demo_images()
vm.analyze(segment_size_mm=(10, 150))
print(vm.results())
vm.plot_analyzed_image()


#%%
# from pylinac.core import image_generator
# import pydicom
# import matplotlib.pyplot as plt
# sim = image_generator.simulators.AS1000Image()
# field_layer = image_generator.layers.FilteredFieldLayer
# image_generator.generate_picketfence(simulator=sim, field_layer=field_layer, file_out='stufff.dcm')
# ds = pydicom.dcmread('stufff.dcm')
# plt.imshow(ds.pixel_array)

#%%
# import io
# import pylinac
# path = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac test files\Picket Fences\AS1000-HD-small.dcm"
# with open(path, 'rb') as f:
#     s = io.BytesIO(f.read())
# pf = pylinac.PicketFence(s, mlc='HD')
# pf.analyze()
# print(pf.results())
# pf.plot_analyzed_image()

#%%
# import pylinac
# leeds = pylinac.LasVegas.from_demo_image()
# print(leeds.image.dpmm)
# leeds.image.plot()

#%%

import pylinac
pylinac.CatPhan504.run_demo()
# ct = pylinac.CatPhan504.from_demo_images()
# ct.analyze()
# ct.publish_pdf('asdf.pdf')
# print(ct.results())
# ct.plot_analyzed_image()

#%%

# for (name, angle) in zip(('g0.dcm', 'g90.dcm', 'g180.dcm', 'g270.dcm'), (0, 90, 180, 270)):
#     as500 = simulators.AS500Image()
#     as500.add_layer(layers.FilteredFieldLayer(field_size_mm=(30, 30)))
#     as500.add_layer(layers.PerfectBBLayer(cax_offset_mm=(0, 0*cos(angle))))
#     as500.add_layer(layers.GaussianFilterLayer())
#     as500.add_layer(layers.RandomNoiseLayer())
#     as500.generate_dicom(gantry_angle=angle, file_out_name='./playground_dir/' + name)
    # plt.title(name)
    # plt.imshow(as500.image)
    # plt.show()


# import pylinac
# path = r"S:\Downloads\WL ix10\20210127"
# # path = './playground_dir'
# wl = pylinac.WinstonLutz(path, use_filenames=True)
# print(wl.results())
# wl.plot_summary()
# wl.plot_images()
# wl.publish_pdf("dead on.pdf")

    # plt.title(name)
    # plt.imshow(as500.image)
    # plt.show()

# fs = pylinac.FlatSym.from_demo_image()
# fs.analyze('varian', 'varian', invert=True)
# fs.plot_analyzed_image()
# from pylinac import py_gui
# py_gui.gui()
# ct = pylinac.PicketFence.from_demo_images()
# ct.analyze()
# ct.plot_analyzed_image()
# path = r"S:\Downloads\RP_all_dyn_MLC.dcm"
# ds = pydicom.dcmread(path)
ttt =1
# pylinac.Starshot.run_demo()
# %%
# tlog = pylinac.TrajectoryLog.from_demo()
# tt = 1
# from pylinac.picketfence import MLCArrangement
# pf = pylinac.PicketFence.from_demo_image()
# folder = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac test files\21-01-18_15-36-15"
# for pdir, sdir, files in os.walk(folder):
#     for file in files:
#         full_path = osp.join(pdir, file)
#         ds = pydicom.dcmread(full_path)
#         plt.imshow(ds.pixel_array)
#         plt.title(ds[0x5000, 0x2500])
#         plt.xlabel(file)
#         plt.show()
# path = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac test files\Picket Fences\AS500#4.dcm"
# with open(path, 'rb') as p, open("dummyfile", 'wb') as w:
#     e = base64.b64encode(p.read())
# log_path = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac\tests_basic\test_files\Picket Fence\PF_log.bin"
# path = r"E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac\tests_basic\test_files\Picket Fence\PF, Elekta, pickets near edges.dcm"
# pf = pylinac.PicketFence.from_demo_image()
# plt.plot(pf.image[200, :])
# # plt.show()
# pf.image.ground()
# pf.image.normalize()
# plt.plot(pf.image[200, :])
# plt.show()
# pf.analyze()
# print(pf.results())
# # pf.publish_pdf('stupid.pdf')
# # print("done w/ pdf")
# pf.plot_analyzed_image(leaf_error_subplot=False)
# print("done w/ plot")
# m = MLC([(10, 10), (40, 5), (10, 10)])
# ttt = 1


# from matplotlib import pyplot as plt
# from pylinac.core.image_generator import AS1200Image
# from pylinac.core.image_generator.layers import FilteredFieldLayer, GaussianFilterLayer
#
# as1200 = AS1200Image()
# height = 350
# width = 4
# offsets = range(-100, 100, 20)
# for offset in offsets:
#     as1200.add_layer(FilteredFieldLayer((height, width), cax_offset_mm=(0, offset)))
# as1200.add_layer(GaussianFilterLayer())
# plt.imshow(as1200.image)
# plt.show()



TEST_DIR = osp.join(osp.dirname(__file__), 'tests_basic', 'test_files', 'Picket Fence')
# path1 = osp.join(TEST_DIR, 'combo-jaw.dcm')
# path2 = osp.join(TEST_DIR, 'combo-mlc.dcm')
# pf = pylinac.PicketFence.from_multiple_images([path1, path2], stretch_each=True)
# pf = pylinac.PicketFence(osp.join(TEST_DIR, 'PF.dcm'), log=osp.join(TEST_DIR, 'PF_log.bin'))
# pf.analyze()
# pf.plot_analyzed_image()
# pf.publish_pdf('deleteme.pdf')
# pf.plot_analyzed_image()
# fs = pylinac.FlatSym.run_demo()
# fs.analyze('varian', 'varian')
# print(fs.results())
# fs.plot_analyzed_image()
# wl = pylinac.WinstonLutz.from_demo_images()
# print(wl.results())
# wl.publish_pdf('deleteme.pdf')
# wl.publish_pdf('deleteme2.pdf')
# from scipy import signal
(TEST_DIR)
#
# img = pylinac.image.load(osp.join(TEST_DIR, 'Starshot#1.tif'))
# xdata = np.linspace(0, 1.7 * np.pi, num=200)
# ydata = signal.sawtooth(xdata, width=0.5)
# prof = pylinac.profile.SingleProfile(ydata)
# print(prof.field_edges())
# prof.filter(size=0.01, kind='gaussian')
# prof.find_valleys(min_distance=0.08)
# prof.plot()
# plt.show()
# fs.image.filter(size=5)
# fs.analyze('varian', 'varian')
# print(fs.results())
# fs.plot_analyzed_image()
ttt =1