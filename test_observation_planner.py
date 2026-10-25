import tempfile
from copy import deepcopy
import unittest
from datetime import datetime
from pathlib import Path

from observation_planner import (
    PlannerStore,
    estimate_sequence_duration,
    export_plan_text,
    field_outline,
    format_duration,
    generate_exposure_script,
    make_mosaic,
    safe_horizon_grid_steps,
    suggested_archive_name,
)
from coverage import planner_field_identity


class ObservationPlannerTests(unittest.TestCase):
    def setUp(self):
        self.profile = {
            "name": "Test camera",
            "width_deg": 10.0,
            "height_deg": 7.0,
            "rotation_deg": 27.0,
            "dither_deg": 0.08,
            "readout_seconds": 2.0,
            "slew_seconds": 30.0,
        }

    def test_rotated_mosaic_crosses_ra_zero_safely(self):
        fields = make_mosaic(359.0, 72.0, self.profile, rows=2, columns=3, overlap_percent=15)
        self.assertEqual(len(fields), 6)
        for field in fields:
            self.assertGreaterEqual(field["ra_deg"], 0)
            self.assertLess(field["ra_deg"], 360)
            self.assertGreaterEqual(field["dec_deg"], -90)
            self.assertLessEqual(field["dec_deg"], 90)
            self.assertEqual(len(field_outline(field)), 33)

    def test_planner_field_identity_deduplicates_loaded_archive_geometry(self):
        field = make_mosaic(15.0, -10.0, self.profile, 1, 1, 10)[0]
        archived_copy = deepcopy(field)
        archived_copy["mosaic_orientation"] = "horizon_eq"
        self.assertEqual(
            planner_field_identity(field), planner_field_identity(archived_copy)
        )

    def test_script_uses_each_fields_dither(self):
        fields = make_mosaic(15.0, -10.0, self.profile, 1, 2, 10)
        fields[1]["profile"]["dither_deg"] = 0.2
        script = generate_exposure_script(fields, "R", 45, 12)
        self.assertIn("loopdisable filter=R for 12", script)
        self.assertEqual(script.count("T0.WOFFS="), 2)
        self.assertIn(" E 45 ", script)

    def test_horizon_grid_spacing_accounts_for_tilted_eq_sensor(self):
        horizontal, vertical = safe_horizon_grid_steps(10, 5, 0, 10)
        self.assertAlmostEqual(horizontal, 9.0)
        self.assertAlmostEqual(vertical, 4.5)
        horizontal, vertical = safe_horizon_grid_steps(10, 5, 45, 10)
        self.assertAlmostEqual(horizontal, 4.5 / (2 ** -0.5))
        self.assertAlmostEqual(vertical, 4.5 / (2 ** -0.5))

    def test_sequence_duration_includes_readout_and_inter_field_slews(self):
        fields = make_mosaic(15.0, -10.0, self.profile, 1, 2, 10)
        duration = estimate_sequence_duration(fields, exposure_seconds=10, loops=3)
        self.assertEqual(duration["frame_count"], 6)
        self.assertEqual(duration["exposure_seconds"], 60)
        self.assertEqual(duration["readout_seconds"], 12)
        # A-B, B-A, A-B, B-A, A-B: five actual field changes.
        self.assertEqual(duration["slew_seconds"], 150)
        self.assertEqual(duration["total_seconds"], 222)
        self.assertEqual(format_duration(duration["total_seconds"]), "3 min 42 s")

    def test_sequence_duration_does_not_slew_repeated_single_field(self):
        fields = make_mosaic(15.0, -10.0, self.profile, 1, 1, 10)
        duration = estimate_sequence_duration(fields, exposure_seconds=10, loops=3)
        self.assertEqual(duration["slew_seconds"], 0)
        self.assertEqual(duration["total_seconds"], 36)

    def test_suggested_archive_name_uses_evening_night_date(self):
        self.assertEqual(
            suggested_archive_name(datetime(2026, 10, 23, 5, 0), ["CTA-S1"], True),
            "2026-10-22_CTA-S1_Morning",
        )
        self.assertEqual(
            suggested_archive_name(datetime(2026, 10, 22, 20, 0), ["CTA S1"], False),
            "2026-10-22_CTA-S1_Evening",
        )

    def test_each_plan_exports_to_its_own_text_file(self):
        fields = make_mosaic(15.0, -10.0, self.profile, 1, 2, 10)
        archive = {
            "name": "2026-10-22_CTA-S1_Morning",
            "filter": "N",
            "exposure_seconds": 30,
            "loops": 15,
            "fields": fields,
        }
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            path = export_plan_text(directory, archive, "2026-10-23 05:00", "Paranal")
            self.assertEqual(Path(path).name, "2026-10-22_CTA-S1_Morning.txt")
            contents = Path(path).read_text(encoding="utf-8")
            self.assertIn("Počet polí: 2", contents)
            self.assertIn("Odhad délky sekvence:", contents)
            self.assertIn("vyčítání:", contents)
            self.assertIn("přejezdy mezi poli:", contents)
            self.assertIn("Pozorovací skript:", contents)
            with self.assertRaises(FileExistsError):
                export_plan_text(directory, archive)

    def test_store_round_trip(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            path = Path(directory) / "planner.json"
            store = PlannerStore(str(path))
            store.data["profiles"][0]["rotation_deg"] = 42.0
            store.save()
            loaded = PlannerStore(str(path))
            self.assertEqual(loaded.data["profiles"][0]["rotation_deg"], 42.0)


if __name__ == "__main__":
    unittest.main()
