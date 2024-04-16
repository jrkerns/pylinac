import copy
from itertools import count, cycle
from math import ceil, floor
from typing import Iterable

import numpy as np
from pydicom import Dataset

from .dicom import BeamType, DicomLinacPlan, FluenceMode, GantryDirection
from .mlc import MLCShaper, inject_sacrifices


def create_dlg_2_side_test(
    nominal_dlg: float = 1.5,
    mlc_type: str = "Millennium",
    gantry_angle: int = 0,
    beam_mu: int = 50,
    dose_rate: int = 600,
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lat: int = 0,
    couch_lng: int = 100,
    couch_rot: int = 0,
    n_gaps: int = 6,
    step_size: float = 0.2,
    position: float = 0,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "DLG",
    beam_name: str = "DLG",
) -> Dataset:
    """Create a single-field DLG evaluation test. This test will create two halves of a field, with the MLCs
    almost perfectly aligned along the CAX inline axis. Some MLCs will be slightly retracted or slightly protruding
    depending on the gap being tested. When the two halves are delivered, there will be some "missing" dose if the gaps
    overlap too much and "extra" dose if the gap has an underlap. By evaluating the dose at the middle of the gap for
    all the DLG test gaps, a linear regression line can be fitted to the gap values. Where this fitted line crosses the
    0-point is the ideal DLG.

    Note that the DLG gaps do not have to be exactly aligned to the ideal gap to determine the ideal DLG. The point is to
    measure several gaps and interpolate. In fact, it's better from an image analysis perspective to **not** have any of the DLG
    gaps be dead on.

    Parameters
    ----------
    nominal_dlg
        This is the DLG that the gaps will be centered about. This should be roughly equal to the DLG of the machine in ARIA.
    n_gaps
        The number of DLG gaps to assess. These will be roughly centered about the nominal DLG parameter.
    step_size
        The size in between each DLG measurement. E.g. 0.2 may result in gaps of 1.2, 1.4, 1.6...
    position
        The X-position of the DLG gaps. Typically this should be 0 (along the CAX).
    """
    # TODO: check ngaps must be even
    # TODO: check total range isn't an obscene value like 500
    # TODO: check n_steps isn't more than the # of leaves, preferably each step is 2 leaves
    # TODO: incorporate MLC config max MLC positions
    mlc = MLCShaper(mlc_type)
    ref_mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = [0, 0.5, 0.5, 1.0]
    field_size = 100

    # create MLC positions
    mlc_positions = [
        mlc.create_dlg_pattern(
            nominal_dlg=nominal_dlg,
            n_gaps=n_gaps,
            step_size=step_size,
            field_size=field_size,
            position=position,
            side="left",
        ),
        mlc.create_dlg_pattern(
            nominal_dlg=nominal_dlg,
            n_gaps=n_gaps,
            step_size=step_size,
            field_size=field_size,
            position=position,
            side="left",
        ),
        mlc.create_dlg_pattern(
            nominal_dlg=nominal_dlg,
            n_gaps=n_gaps,
            step_size=step_size,
            field_size=field_size,
            position=position,
            side="right",
        ),
        mlc.create_dlg_pattern(
            nominal_dlg=nominal_dlg,
            n_gaps=n_gaps,
            step_size=step_size,
            field_size=field_size,
            position=position,
            side="right",
        ),
    ]
    ref_mlc_positions = [
        ref_mlc.add_rectangle(
            left_position=-field_size / 2,
            right_position=field_size / 2,
            top_position=mlc.mlc.max_y,
            bottom_position=-mlc.mlc.max_y,
        ),
    ] * 4

    # create jaw positions
    x1 = -field_size / 2 - jaw_padding
    x2 = field_size / 2 + jaw_padding
    y1 = -field_size
    y2 = field_size

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    # create DMLC field
    v.add_beam(
        beam_name=beam_name,
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=mlc_positions,
        meter_sets=meter_sets,
        beam_mu=beam_mu,
        fluence_mode=fluence_mode,
    )
    # create reference field
    v.add_beam(
        beam_name=f"{beam_name} Ref",
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=ref_mlc_positions,
        meter_sets=meter_sets,
        beam_mu=beam_mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


def create_dose_rate_test(
    dose_rates: tuple = (100, 200, 300, 400, 500, 600),
    default_dose_rate: int = 600,
    mlc_type: str = "Millennium",
    mlc_speed_max: int = 25,
    gantry_angle: int = 0,
    desired_mu: int = 50,
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lat: int = 0,
    couch_lng: int = 100,
    couch_rot: int = 0,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "DoseRateROI",
) -> Dataset:
    """Create a single-field test to deliver several different dose rates of the machine. This creates several ROIs and
    all of the ROIs should ideally give the same signal response (after accounting for the natural beam shape via reference).
    The dose rate slow-down is accomplished via "sacrificial" moves of the outside MLC pair, which hinder the machine
    from delivering at its maximum.

    Unlike the Varian beams, this **only** tests the dose rate and the gantry is static. Additionally, the desired dose
    rates can be configured unlike Varian.

    The goal of this test is the same as the "Dose Rate Linearity Test". The difference between the two is that this is
    one field + reference, thus quicker. The latter delivers several open fields. The latter is easier to process when using primitive
    image analysis techniques and does not need a reference image and also corresponds to what the mental default of a
    dose rate analysis test is for a physicist where an ion chamber is normally used. E.g. a physicist can directly compare
    the results of the latter test with the EPID vs an ion chamber easily. This current test is proposed as a quicker
    alternative.
    """
    mlc = MLCShaper(mlc_type)
    ref_mlc = MLCShaper(mlc_type)
    # create MU weighting
    start_meter_sets = np.linspace(0, 1, len(dose_rates) + 1)
    end_meter_sets = np.linspace(0, 1, len(dose_rates) + 1)
    meter_sets = np.concatenate((start_meter_sets, end_meter_sets))[
        :-1
    ]  # drop last one as it's a 1.00 repeat
    meter_sets.sort()

    # TODO: add validation of passed speed ROIs fit in MLC window
    roi_size = 20  # mm wide

    # calculate MU
    mlc_transition_time = roi_size / mlc_speed_max
    min_mu = mlc_transition_time * max(dose_rates) * len(dose_rates) / 60
    mu = max(desired_mu, ceil(min_mu))

    # create MLC positions
    start_pos = -len(dose_rates) / 2 * roi_size  # center the ROIs about the CAX
    center_pos = [
        idx * roi_size + start_pos + roi_size / 2 for idx in range(len(dose_rates))
    ]
    times_to_transition = [
        mu * 60 / (dose_rate * len(dose_rates)) for dose_rate in dose_rates
    ]
    # start off at the left-most edge
    mlc_positions = [mlc.create_vmat_start_position(start_pos)]
    # TODO: use MLCconfig max for below
    ref_mlc_positions = [
        ref_mlc.create_vmat_position(start_pos, y_upper=200, y_lower=-200)
    ]
    sacrificial_movements = [tt * mlc_speed_max for tt in times_to_transition]
    for center, sacrifice, _ in zip(center_pos, sacrificial_movements, meter_sets):
        # two identical positions; one is the start and the other is the end. MU is delivered in between.
        mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center, strip_width=roi_size, sacrifice=0
            )
        )  # cut off top and bottom pair as they are sacrificial
        mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center, strip_width=roi_size, sacrifice=sacrifice
            )
        )  # cut off top and bottom pair as they are sacrificial
        # TODO: use MLCconfig max/min
        ref_mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center,
                strip_width=roi_size,
                y_upper=200,
                y_lower=-200,
            )
        )
        ref_mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center,
                strip_width=roi_size,
                y_upper=200,
                y_lower=-200,
            )
        )

    # break up sacrifices
    new_meter_sets, new_mlc_positions = inject_sacrifices(
        meter_sets,
        mlc,
        mlc_positions,
        sacrificial_movements,
        sacrifice_start_pos=start_pos,
        sacrifice_max=50,
    )

    # create jaw positions
    x1 = min(center_pos) - roi_size / 2 - jaw_padding  # give mm extra room
    x2 = max(center_pos) + roi_size / 2 + jaw_padding
    y1 = -mlc.mlc.max_y
    y2 = mlc.mlc.max_y

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    # dmlc beam
    v.add_beam(
        beam_name=f"DR{min(dose_rates)},{max(dose_rates)}",
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max(dose_rates),
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=new_mlc_positions,
        meter_sets=new_meter_sets,
        beam_mu=mu,
        fluence_mode=fluence_mode,
    )
    # ref beam
    v.add_beam(
        beam_name="DR Ref",
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=default_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=ref_mlc_positions,
        meter_sets=meter_sets,
        beam_mu=mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


def create_gantry_speed_test(
    speeds: tuple = (2, 3, 4, 4.8),
    mlc_type: str = "Millennium",
    mlc_speed_max: int = 25,
    max_dose_rate: int = 400,
    start_gantry_angle: int = 179,
    max_gantry_speed: float = 4.8,
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lat: int = 0,
    couch_lng: int = 100,
    couch_rot: int = 0,
    beam_name: str = "GS",
    gantry_rot_dir: GantryDirection = GantryDirection.CLOCKWISE,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "GSpeed",
) -> Dataset:
    """Create a single-field gantry speed test plan. This plan will slow down the gantry to the desired speeds and deliver
    an ROI, similar to the classical VMAT tests, but this allows the user to specify the speeds. Inevitably, the dose
    rate will also vary, but it is a dependent variable in the dicom plan generation. The gantry speeds are achieved
    via sacrificial MLC movements of the outside pairs."""
    gantry_rot_dir = (
        GantryDirection(gantry_rot_dir)
        if isinstance(gantry_rot_dir, str)
        else gantry_rot_dir
    )
    mlc = MLCShaper(mlc_type)
    ref_mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = np.linspace(0, 1, len(speeds) * 6 + 1)
    ref_meter_sets = copy.copy(meter_sets)

    # TODO: add validation of passed speed ROIs fit in MLC window
    # TODO: validate sacrificial movements (esp slow speeds) don't go over max overreach of MLCs (~14cm?)

    # create MLC positions
    roi_size = 30  # mm wide
    step_size = 10  # control point step size. I.e. 10mm to 20mm to 30mm, etc
    start_pos = -len(speeds) / 2 * roi_size  # center the ROIs
    strip_widths = cycle((10, 20, roi_size, 20, 10, 1))  # TODO: make this parametric
    center_pos = count(start_pos + step_size / 2, step_size / 2)
    time_to_finish = [max_gantry_speed / speed for speed in np.repeat(speeds, 6)]
    sacrificial_movements = [ttf * mlc_speed_max for ttf in time_to_finish]
    # start off at the left-most edge
    mlc_positions = [mlc.create_vmat_start_position(start_pos, strip_width=1)]
    ref_mlc_positions = [
        ref_mlc.create_vmat_position(
            start_pos, y_upper=mlc.mlc.max_y, y_lower=-mlc.mlc.max_y, strip_width=1
        )
    ]
    for width, center, sacrifice, _ in zip(
        strip_widths, center_pos, sacrificial_movements, meter_sets
    ):
        mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center, strip_width=width, sacrifice=sacrifice
            )
        )  # cut off top and bottom pair as they are sacrificial
        ref_mlc_positions.append(
            ref_mlc.create_vmat_position(
                x_infield_position=center,
                strip_width=width,
                y_upper=mlc.mlc.max_y,
                y_lower=-mlc.mlc.max_y,
            )
        )

    # break up sacrifices
    # meter_sets, mlc_positions = inject_sacrifices(meter_sets, mlc, mlc_positions, sacrificial_movements,
    #                                               sacrifice_start_pos=start_pos, sacrifice_max=50)

    # calc max allowable MU
    min_ttf = min(time_to_finish)
    meter_diff = 1 / (len(meter_sets) - 1)
    max_mu = floor((min_ttf * (max_dose_rate / 60)) / meter_diff)

    # create jaw positions
    x1 = min(mlc_positions[0]) - roi_size / 2 - jaw_padding
    x2 = np.max(mlc_positions[-1]) + roi_size / 2 + jaw_padding
    y1 = -mlc.mlc.max_y
    y2 = mlc.mlc.max_y

    # create gantry angles
    gantry_angles = [
        start_gantry_angle - idx * max_gantry_speed for idx in range(len(meter_sets))
    ]

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    # add the DMLC beam
    v.add_beam(
        beam_name=beam_name,
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angles,
        gantry_direction=gantry_rot_dir,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=mlc_positions,
        meter_sets=meter_sets,
        beam_mu=max_mu,
        fluence_mode=fluence_mode,
    )
    # add the reference beam (same movements, but static gantry)
    v.add_beam(
        beam_name="Reference" + beam_name,
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angles[-1],
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=ref_mlc_positions,
        meter_sets=ref_meter_sets,
        beam_mu=max_mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


def create_mlc_speed_test(
    speeds: tuple = (17.14, 10, 15, 20),
    mlc_type: str = "Millennium",
    mlc_speed_max: int = 25,
    max_dose_rate: int = 400,
    gantry_angle: int = 0,
    energy: int = 6,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lat: int = 0,
    couch_lng: int = 100,
    couch_rot: int = 0,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "MLC Speed",
    beam_name: str = "MLC Speed",
) -> Dataset:
    """Create an MLC speed test. This is similar to the standard Varian DRMLC test where multiple ROIs are created,
    but the gantry isn't moving.
    The user can also configure the speeds desired. The speeds are achieved via sacrificial moves of the outside MLC pairs.
    """
    mlc = MLCShaper(mlc_type)
    ref_mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = np.linspace(0, 1, len(speeds) * 6 + 1)

    # TODO: add validation of passed speed ROIs fit in MLC window
    # TODO: validate sacrificial movements (esp slow speeds) don't go over max overreach of MLCs (~14cm?)

    # create MLC positions
    roi_size = 30  # mm wide
    step_size = 10  # control point step size. I.e. 10mm to 20mm to 30mm, etc
    strip_width = 1
    start_pos = -len(speeds) / 2 * roi_size  # center the ROIs
    strip_widths = cycle((10, 20, 30, 20, 10, strip_width))
    center_pos = count(
        start_pos + step_size / 2, step_size / 2
    )  # TODO: change to add spacing between ROIs to minimize spikes on profile
    time_to_finish = [step_size / speed for speed in np.repeat(speeds, 6)]
    sacrificial_movements = [ttf * mlc_speed_max for ttf in time_to_finish]
    # start off at the left-most edge
    mlc_positions = [
        mlc.create_vmat_start_position(start_pos, strip_width=strip_width, to_left=True)
    ]
    ref_mlc_position = [
        ref_mlc.create_vmat_start_position(
            start_pos, strip_width=strip_width, to_left=True
        )
    ]
    for width, center, sacrifice in zip(
        strip_widths, center_pos, sacrificial_movements
    ):
        mlc_positions.append(
            mlc.create_vmat_position(
                x_infield_position=center, strip_width=width, sacrifice=sacrifice
            )
        )  # cut off top and bottom pair as they are sacrificial
        ref_mlc_position.append(
            ref_mlc.create_vmat_position(
                x_infield_position=center,
                strip_width=width,
                y_upper=mlc.mlc.max_y,
                y_lower=-mlc.mlc.max_y,
            )
        )

    # break up sacrifices
    meter_sets, mlc_positions = inject_sacrifices(
        meter_sets,
        mlc,
        mlc_positions,
        sacrificial_movements,
        sacrifice_start_pos=start_pos,
        sacrifice_max=50,
    )

    # calc max allowable MU
    min_ttf = min(time_to_finish)
    meter_diff = 1 / (len(meter_sets) - 1)
    max_mu = floor((min_ttf * (max_dose_rate / 60)) / meter_diff)

    # create jaw positions
    x1 = min(mlc_positions[0]) - roi_size / 2 - jaw_padding
    x2 = np.max(mlc_positions[-1]) + roi_size / 2 + jaw_padding
    y1 = -mlc.mlc.max_y
    y2 = mlc.mlc.max_y

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    # add DMLC beam
    v.add_beam(
        beam_name=beam_name,
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=mlc_positions,
        meter_sets=meter_sets,
        beam_mu=max_mu,
        fluence_mode=fluence_mode,
    )
    # add reference beam
    v.add_beam(
        beam_name="MLC Sp Ref",
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=ref_mlc_position,
        meter_sets=meter_sets,
        beam_mu=max_mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


# TODO: add VMAT-type option
def create_picketfence_test(
    strip_width: float = 3,
    strip_positions: tuple = (-50, -30, -10, 10, 30, 50),
    mlc_type: str = "Millennium",
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    max_dose_rate: int = 600,
    gantry_angle: int = 0,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lng: int = 100,
    couch_lat: int = 0,
    couch_rot: int = 0,
    beam_mu: int = 200,
    transition_dose: float = 0.01,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "PicketF",
    beam_name: str = "PF",
) -> Dataset:
    """Create a typical picket fence test. This will create N pickets at the given positions with a specified width.
    The gantry is static for this delivery."""
    mlc = MLCShaper(mlc_type)
    # create MU weighting
    start_meter_sets = np.linspace(0, 1, len(strip_positions) + 1)
    end_meter_sets = np.linspace(0, 1, len(strip_positions) + 1)
    meter_sets = np.concatenate((start_meter_sets, end_meter_sets))[:-1]
    meter_sets.sort()
    meter_sets[1:-1:2] = [
        s + transition_dose for idx, s in enumerate(meter_sets[1:-1:2])
    ]

    # create MLC positions
    # TODO: validate this doesn't go outside MLC window
    mlc_positions = [
        mlc.add_strip(
            position=strip_positions[0] - 10,
            y_bound_upper=mlc.mlc.max_y,
            y_bound_lower=-mlc.mlc.max_y,
            strip_width=strip_width,
        ),
    ]  # align just to the left of first picket
    for pos in strip_positions:
        mlc_positions.append(
            mlc.add_strip(
                position=pos,
                y_bound_lower=-mlc.mlc.max_y,
                y_bound_upper=mlc.mlc.max_y,
                strip_width=strip_width,
            )
        )
        mlc_positions.append(
            mlc.add_strip(
                position=pos,
                y_bound_lower=-mlc.mlc.max_y,
                y_bound_upper=mlc.mlc.max_y,
                strip_width=strip_width,
            )
        )

    # create jaw positions
    x1 = min(strip_positions) - strip_width / 2 - jaw_padding
    x2 = max(strip_positions) + strip_width / 2 + jaw_padding
    y1 = -mlc.mlc.max_y
    y2 = mlc.mlc.max_y

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    v.add_beam(
        beam_name=beam_name,
        beam_type=BeamType.DYNAMIC,
        energy=energy,
        dose_rate=max_dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=mlc_positions,
        meter_sets=meter_sets,
        beam_mu=beam_mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


def create_open_field(
    x1: float,
    x2: float,
    y1: float,
    y2: float,
    defined_by_mlcs: bool = True,
    mlc_type: str = "Millennium",
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    dose_rate: int = 600,
    gantry_angle: int = 0,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lng: int = 100,
    couch_lat: int = 0,
    couch_rot: int = 0,
    beam_mu: int = 200,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "Open Field",
    beam_name: str = "Open",
    machine_sn: str = "H191111",
    tolerance_table: str = "T1",
    patient_name: str = "RadMachine",
    patient_id: str = "rad",
) -> Dataset:
    """Create a simple, static open field defined either by the jaws or MLCs."""
    mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = [0, 1]

    # create MLC positions
    if defined_by_mlcs:
        # TODO: validate this doesn't go outside MLC window
        # create MLC square
        mlc_positions = mlc.add_rectangle(
            left_position=x1,
            right_position=x2,
            top_position=y2,
            bottom_position=y1,
            outer_strip_width=3,
            x_outfield_position=x1 - jaw_padding * 2,
        )
        # move jaws just outside MLC window
        x1 -= jaw_padding
        x2 += jaw_padding
        y1 -= jaw_padding
        y2 += jaw_padding
    else:
        # park the MLCs
        mlc_positions = mlc.add_rectangle(
            left_position=-mlc.mlc.max_x, right_position=mlc.mlc.max_x
        )

    # create plan
    v = DicomLinacPlan(
        machine_name=machine_name,
        plan_label=plan_label,
        patient_name=patient_name,
        patient_id=patient_id,
        machine_sn=machine_sn,
        tolerance_table=tolerance_table,
    )
    v.add_beam(
        beam_name=beam_name,
        beam_type=BeamType.STATIC,
        energy=energy,
        dose_rate=dose_rate,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        gantry_angles=gantry_angle,
        gantry_direction=GantryDirection.NONE,
        coll_angle=coll_angle,
        couch_vrt=couch_vrt,
        couch_lat=couch_lat,
        couch_lng=couch_lng,
        couch_rot=couch_rot,
        mlc_positions=[
            mlc_positions,
        ],
        meter_sets=meter_sets,
        beam_mu=beam_mu,
        fluence_mode=fluence_mode,
    )
    return v.ds


def create_MU_linearity_set(
    x1: float = -100,
    x2: float = 100,
    y1: float = -100,
    y2: float = 100,
    mu_set: tuple = (5, 50, 100, 200, 500, 1000),
    mlc_type: str = "Millennium",
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    dose_rate: int = 600,
    gantry_angle: int = 0,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lng: int = 100,
    couch_lat: int = 0,
    couch_rot: int = 0,
    machine_name: str = "TrueBeam",
    plan_label: str = "MU Linearity",
) -> Dataset:
    """Create a set of open beams that deliver varying amounts of MUs."""
    mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = [0, 1]

    # park the MLCs
    mlc_positions = mlc.add_rectangle(
        left_position=-mlc.mlc.max_x, right_position=mlc.mlc.max_x
    )

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    for mu in mu_set:
        v.add_beam(
            beam_name=f"{mu} MU",
            beam_type=BeamType.STATIC,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=[
                mlc_positions,
            ],
            meter_sets=meter_sets,
            beam_mu=mu,
            fluence_mode=fluence_mode,
        )
    return v.ds


def create_dose_rate_linearity_set(
    x1: float = -100,
    x2: float = 100,
    y1: float = -100,
    y2: float = 100,
    dr_set: tuple = (5, 50, 100, 200, 400, 600),
    mlc_type: str = "Millennium",
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    mu: int = 100,
    gantry_angle: int = 0,
    coll_angle: int = 0,
    couch_vrt: int = 0,
    couch_lng: int = 100,
    couch_lat: int = 0,
    couch_rot: int = 0,
    machine_name: str = "TrueBeam",
    plan_label: str = "DR Linearity",
) -> Dataset:
    """Create a set of open beams that deliver the same dose but with varying dose rates. See the other dose rate test for more info."""
    mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = [0, 1]

    # park the MLCs
    mlc_positions = mlc.add_rectangle(
        left_position=-mlc.mlc.max_x, right_position=mlc.mlc.max_x
    )
    # mlc_positions = mlc.create_rectangle(left_position=-200, right_position=50, top_position=50, bottom_position=-50)

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    for dr in dr_set:
        v.add_beam(
            beam_name=f"{dr} MU/min",
            beam_type=BeamType.STATIC,
            energy=energy,
            dose_rate=dr,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=[
                mlc_positions,
            ],
            meter_sets=meter_sets,
            beam_mu=mu,
            fluence_mode=fluence_mode,
        )
    return v.ds


def create_winston_lutz_set(
    x1: float = -10,
    x2: float = 10,
    y1: float = -10,
    y2: float = 10,
    defined_by_mlcs: bool = True,
    mlc_type: str = "Millennium",
    energy: int = 6,
    fluence_mode: FluenceMode = FluenceMode.STANDARD,
    dose_rate: int = 600,
    axes_positions: Iterable[dict] = ({"gantry": 0, "collimator": 0, "couch": 0},),
    couch_vrt: int = 0,
    couch_lng: int = 100,
    couch_lat: int = 0,
    beam_mu: int = 10,
    jaw_padding: int = 5,
    machine_name: str = "TrueBeam",
    plan_label: str = "Winston Lutz",
) -> Dataset:
    """Create a set of simple, static open fields, defined either by the jaws or MLCs with G/C/P settings as done in WL."""
    mlc = MLCShaper(mlc_type)
    # create MU weighting
    meter_sets = [0, 1]

    # create MLC positions
    if defined_by_mlcs:
        # TODO: validate this doesn't go outside MLC window
        # create MLC square
        mlc_positions = mlc.add_rectangle(
            left_position=x1,
            right_position=x2,
            top_position=y2,
            bottom_position=y1,
            outer_strip_width=3,
            x_outfield_position=x1 - jaw_padding * 2,
        )
        # move jaws just outside MLC window
        x1 -= jaw_padding
        x2 += jaw_padding
        y1 -= jaw_padding
        y2 += jaw_padding
    else:
        # park the MLCs
        mlc_positions = mlc.add_rectangle(
            left_position=-mlc.mlc.max_x, right_position=mlc.mlc.max_x
        )

    # create plan
    v = DicomLinacPlan(machine_name=machine_name, plan_label=plan_label)
    # create beams
    for axes in axes_positions:
        v.add_beam(
            beam_name=f"G{axes['gantry']:g}C{axes['collimator']:g}P{axes['couch']:g}",
            beam_type=BeamType.STATIC,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=axes["gantry"],
            gantry_direction=GantryDirection.NONE,
            coll_angle=axes["collimator"],
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=axes["couch"],
            mlc_positions=[
                mlc_positions,
            ],
            meter_sets=meter_sets,
            beam_mu=beam_mu,
            fluence_mode=fluence_mode,
        )
    return v.ds
