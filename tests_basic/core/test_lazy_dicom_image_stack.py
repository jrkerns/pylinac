from __future__ import annotations

import tempfile
from pathlib import Path
from unittest import TestCase

import numpy as np
import pydicom

from pylinac.core.array_utils import create_dicom_files_from_3d_array
from pylinac.core.image import LazyDicomImageStack


class TestLazyDicomImageStackUidFiltering(TestCase):
    def test_check_uid_preserves_path_metadata_pairing_for_interleaved_series(self):
        """Filtering to the dominant UID should keep each path paired with its metadata."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary_dir = create_dicom_files_from_3d_array(
                np.stack(
                    [
                        np.full((2, 2), fill_value=value, dtype=np.uint16)
                        for value in (10, 20, 30)
                    ],
                    axis=-1,
                ),
                out_dir=root / "primary",
            )
            secondary_dir = create_dicom_files_from_3d_array(
                np.stack(
                    [
                        np.full((2, 2), fill_value=value, dtype=np.uint16)
                        for value in (100, 200)
                    ],
                    axis=-1,
                ),
                out_dir=root / "secondary",
            )

            interleaved_paths = [
                primary_dir / "0.dcm",
                secondary_dir / "0.dcm",
                primary_dir / "1.dcm",
                secondary_dir / "1.dcm",
                primary_dir / "2.dcm",
            ]

            stack = LazyDicomImageStack(interleaved_paths, min_number=3, check_uid=True)

            expected_paths = [
                primary_dir / "0.dcm",
                primary_dir / "1.dcm",
                primary_dir / "2.dcm",
            ]
            self.assertEqual(
                [Path(path) for path in stack._image_path_keys], expected_paths
            )
            self.assertEqual(len(stack.metadatas), len(expected_paths))

            for path, metadata in zip(stack._image_path_keys, stack.metadatas):
                ds = pydicom.dcmread(path, force=True, stop_before_pixels=True)
                self.assertEqual(ds.SeriesInstanceUID, metadata.SeriesInstanceUID)
                self.assertEqual(ds.SOPInstanceUID, metadata.SOPInstanceUID)
