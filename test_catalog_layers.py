from datetime import datetime
import json

import numpy as np
from astropy.time import Time

from coverage import (
    MainWindow, PLANET_STYLES, PointingDataReader, SupplementCatalogReader,
    format_moon_readout, format_sun_altitude, moonlight_background_active,
)


def test_catalog_coordinate_and_comet_element_helpers():
    assert np.isclose(SupplementCatalogReader.sexagesimal("12:00:00", hours=True), 180.0)
    assert np.isclose(SupplementCatalogReader.sexagesimal("-24:30:00"), -24.5)

    perihelion = Time("2026-10-22", scale="tt").jd
    vector, radius = MainWindow.comet_heliocentric_ecliptic({
        "tp_jd": perihelion,
        "q": 1.25,
        "e": 1.0,
        "arg_peri": 0.0,
        "node": 0.0,
        "inclination": 0.0,
    }, perihelion)
    assert np.allclose(vector, [1.25, 0.0, 0.0])
    assert np.isclose(radius, 1.25)


def test_skycoverage_window_is_anchored_to_displayed_night():
    records = [
        (datetime(2026, 9, 30), "old", "500"),
        (datetime(2026, 10, 1), "newest", "500"),
    ]
    assert PointingDataReader.files_in_date_window(
        records, datetime(2026, 10, 22), 10
    ) == []


def test_sun_altitude_readout_describes_horizon_side():
    assert format_sun_altitude(12.34) == "Slunce: +12.3° nad obzorem"
    assert format_sun_altitude(-13.06) == "Slunce: -13.1° pod obzorem"
    assert format_sun_altitude(0.01) == "Slunce: 0.0° na obzoru"


def test_moon_readout_contains_altitude_and_illuminated_percent():
    assert format_moon_readout(18.26, 0.734) == "Měsíc: +18.3° | fáze: 73 %"


def test_moonlight_background_requires_moon_above_horizon_away_from_new():
    assert moonlight_background_active(12.0, 3.1)
    assert moonlight_background_active(12.0, 26.0)
    assert not moonlight_background_active(-0.1, 14.0)
    assert not moonlight_background_active(12.0, 2.9)
    assert not moonlight_background_active(12.0, 26.6)


def test_mpc_packed_date_and_bright_asteroid_parser(tmp_path):
    line = list(" " * 195)

    def put(start, end, value):
        line[start:end] = list(str(value).rjust(end - start))

    put(0, 7, "00004")
    put(8, 13, "3.25")
    put(14, 19, "0.15")
    line[20:25] = list("K25BL")
    put(26, 35, "26.80968")
    put(37, 46, "151.53712")
    put(48, 57, "103.70232")
    put(59, 68, "7.14406")
    put(70, 79, "0.0901676")
    put(92, 103, "2.3615413")
    line[179:188] = list("(4) Vesta")
    path = tmp_path / "bright.txt"
    path.write_text("".join(line) + "\n", encoding="ascii")

    objects = SupplementCatalogReader.load_bright_asteroids(path)
    assert len(objects) == 1
    assert objects[0]["number"] == 4
    assert np.isclose(objects[0]["a"], 2.3615413)
    assert Time(objects[0]["epoch_jd"], format="jd").datetime.date().isoformat() == "2025-11-21"


def test_asteroid_propagation_at_epoch_for_circular_orbit():
    vector, radius = MainWindow.asteroid_heliocentric_ecliptic({
        "a": 2.0,
        "e": 0.0,
        "epoch_jd": 2460000.5,
        "mean_anomaly": 0.0,
        "arg_peri": 0.0,
        "node": 0.0,
        "inclination": 0.0,
    }, 2460000.5)
    assert np.allclose(vector, [2.0, 0.0, 0.0])
    assert np.isclose(radius, 2.0)


def test_milky_way_geojson_loader_keeps_contour_levels(tmp_path):
    path = tmp_path / "mw.json"
    path.write_text(json.dumps({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {},
            "geometry": {
                "type": "MultiPolygon",
                "coordinates": [[[
                    [-10.0, -2.0], [0.0, 3.0], [10.0, -2.0], [-10.0, -2.0]
                ]]],
            },
        }],
    }), encoding="utf-8")
    levels = SupplementCatalogReader.load_milky_way_contours(path)
    assert len(levels) == 1
    assert levels[0][0][0] == (350.0, -2.0)


def test_planet_layer_uses_short_labels_and_omits_pluto():
    assert list(PLANET_STYLES) == [
        "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune",
    ]
    assert [style[0] for style in PLANET_STYLES.values()] == [
        "Me", "Ve", "Ma", "Ju", "Sa", "Ur", "Ne",
    ]
