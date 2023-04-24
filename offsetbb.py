import pylinac
from pylinac.core.image_generator import simulators, layers, generate_winstonlutz_multi_bb_multi_field

wl_dir = 'wl_dir'
generate_winstonlutz_multi_bb_multi_field(
        simulator=simulators.AS1200Image(sid=1000),
        field_layer=layers.PerfectFieldLayer,
        final_layers=[layers.GaussianFilterLayer(sigma_mm=1),],
        dir_out=wl_dir,
        field_offsets=((0, 0, 0), (20, -20, 60)),
        field_size_mm=(20, 20),
        bb_offsets=[[1, 0, 0], [19, -20, 60]],  # here's the offset
)
arrange = (
    {'name': 'Iso', 'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 5, 'rad_size_mm': 20},
    {'name': 'Left,Down,In', 'offset_left_mm': 20, 'offset_up_mm': -20, 'offset_in_mm': 60, 'bb_size_mm': 5, 'rad_size_mm': 20},)

wl = pylinac.WinstonLutzMultiTargetMultiField(wl_dir)
wl.analyze(bb_arrangement=arrange)
print(wl.results())