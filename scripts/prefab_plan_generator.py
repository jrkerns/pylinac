from pathlib import Path

from pylinac.plan_generator.dicom import FluenceMode, PlanGenerator

rt_plan_dir = (
    Path(__file__).parent.parent.absolute()
    / "Files"
    / "RapidArc QA Test Procedures and Files for TrueBeam"
    / "Plans"
)

REVISION = 2
MLCS = (
    ("HD120", rt_plan_dir / "HDMLC" / "T0.2_PicketFenceStatic_HD120_TB_Rev02.dcm"),
    ("M120", rt_plan_dir / "Millennium" / "T0.2_PicketFenceStatic_M120_TB_Rev02.dcm"),
)
CONTEXTS = [
    # energy, mode, max dose rate
    (6, FluenceMode.STANDARD, 600),
    (6, FluenceMode.FFF, 1400),
    (10, FluenceMode.STANDARD, 600),
    (10, FluenceMode.FFF, 2400),
    (15, FluenceMode.STANDARD, 600),
    (18, FluenceMode.STANDARD, 600),
]

for mlc, rt_plan_file in MLCS:
    for context in CONTEXTS:
        energy, fluence_mode, max_dose_rate = context
        plan_name = f"{energy}MV{'' + fluence_mode.value if fluence_mode != FluenceMode.STANDARD else ''}"
        generator = PlanGenerator.from_rt_plan_file(
            rt_plan_file, plan_name=plan_name, plan_label=plan_name
        )
        # picket fence
        for gantry_angle in (0, 90, 180, 270):
            generator.add_picketfence_beam(
                strip_width_mm=3,
                strip_positions_mm=(-60, -30, 0, 30, 60),  # 5 pickets
                mu=100,
                beam_name=f"PF 3mm G{gantry_angle}",
                gantry_angle=gantry_angle,
                fluence_mode=fluence_mode,
                energy=energy,
                dose_rate=max_dose_rate,
            )
        # winston lutz
        generator.add_winston_lutz_beams(
            axes_positions=(
                # basic gantry
                {"gantry": 0, "collimator": 0, "couch": 0},
                {"gantry": 90, "collimator": 0, "couch": 0},
                {"gantry": 180, "collimator": 0, "couch": 0},
                {"gantry": 270, "collimator": 0, "couch": 0},
                # collimator-rotation
                {"gantry": 0, "collimator": 270, "couch": 0},
                {"gantry": 0, "collimator": 90, "couch": 0},
                {"gantry": 0, "collimator": 315, "couch": 0},
                {"gantry": 0, "collimator": 45, "couch": 0},
                # couch-rotation
                {"gantry": 0, "collimator": 0, "couch": 45},
                {"gantry": 0, "collimator": 0, "couch": 90},
                {"gantry": 0, "collimator": 0, "couch": 315},
                {"gantry": 0, "collimator": 0, "couch": 270},
                # combo
                {"gantry": 45, "collimator": 15, "couch": 15},
                {"gantry": 5, "collimator": 330, "couch": 60},
                {"gantry": 330, "collimator": 350, "couch": 350},
            ),
            x1=-10,
            x2=10,
            y1=-10,
            y2=10,
            defined_by_mlcs=True,
            mu=5,
            energy=energy,
            dose_rate=max_dose_rate,
            fluence_mode=fluence_mode,
        )

        # dose rate
        generator.add_dose_rate_beams(
            dose_rates=(100, 200, 400, 600),
            y1=-50,
            y2=50,
            default_dose_rate=max_dose_rate,
            desired_mu=100,
            fluence_mode=fluence_mode,
            energy=energy,
        )

        # MLC speed
        generator.add_mlc_speed_beams(
            speeds=(5, 10, 15, 20),
            roi_size_mm=20,
            y1=-50,
            y2=50,
            mu=100,
            default_dose_rate=max_dose_rate,
            fluence_mode=fluence_mode,
            energy=energy,
        )

        # gantry speed
        generator.add_gantry_speed_beams(
            speeds=(1, 2, 3, 4),
            max_dose_rate=max_dose_rate,
            start_gantry_angle=179,
            roi_size_mm=20,
            y1=-50,
            y2=50,
            mu=100,
            fluence_mode=fluence_mode,
            energy=energy,
        )

        # save file
        file_name = f"R{REVISION}_{plan_name}_{mlc}_prefab.dcm".replace(" ", "_")
        generator.to_file(filename=file_name)
        print(f"Generated {file_name} plan")
