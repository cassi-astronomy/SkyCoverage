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


def _mpcorb_line(packed, h, e, a, readable, mean_anomaly="26.80968"):
    line = list(" " * 195)

    def put(start, end, value):
        line[start:end] = list(str(value).rjust(end - start))

    put(0, 7, packed)
    put(8, 13, h)
    put(14, 19, "0.15")
    line[20:25] = list("K25BL")
    put(26, 35, mean_anomaly)
    put(37, 46, "151.53712")
    put(48, 57, "103.70232")
    put(59, 68, "7.14406")
    put(70, 79, e)
    put(92, 103, a)
    line[166:166 + len(readable)] = list(readable)
    return "".join(line).rstrip() + "\n"


def test_asteroid_parser_keeps_names_and_unnumbered_objects(tmp_path):
    path = tmp_path / "mpcorb.txt"
    path.write_text(
        _mpcorb_line("00004", "3.25", "0.0901676", "2.3615413", "(4) Vesta")
        + _mpcorb_line("K24A01B", "12.10", "0.1500000", "2.5000000", "2024 AB1"),
        encoding="ascii",
    )
    vesta, unnumbered = SupplementCatalogReader.load_bright_asteroids(path)
    assert (vesta["number"], vesta["name"]) == (4, "Vesta")
    assert (unnumbered["number"], unnumbered["name"]) == (None, "2024 AB1")
    label = SupplementCatalogReader.asteroid_label
    assert label(vesta) == "4"
    assert label(vesta, detailed=True, magnitude=6.04) == "4 Vesta  6.0m"
    assert label(unnumbered, detailed=True, magnitude=13.96) == "2024 AB1  14.0m"


def test_mpcorb_filter_keeps_only_orbits_that_can_reach_limit(tmp_path):
    import gzip as gzip_module

    bright_inner = _mpcorb_line("00004", "3.25", "0.0901676", "2.3615413", "(4) Vesta")
    faint_inner = _mpcorb_line("99999", "16.00", "0.1000000", "2.3000000", "(99999) Faint")
    faint_outer = _mpcorb_line("88888", "11.50", "0.0500000", "3.2000000", "(88888) Outer")
    gz_path = tmp_path / "MPCORB.DAT.gz"
    with gzip_module.open(gz_path, "wt", encoding="ascii") as output:
        output.write("MPCORB header\n-----------\n")
        output.write(bright_inner + faint_inner + faint_outer)
    output_path = tmp_path / "bright.txt"
    count = SupplementCatalogReader.write_bright_mpcorb(
        str(gz_path), str(output_path), 14.5
    )
    assert count == 1
    objects = SupplementCatalogReader.load_bright_asteroids(output_path)
    assert [item["name"] for item in objects] == ["Vesta"]


def test_vectorised_asteroid_propagation_matches_scalar_version():
    from orbital_mechanics import (
        asteroid_heliocentric_ecliptic, asteroid_heliocentric_ecliptic_array,
        build_asteroid_arrays, calendar_jd,
    )

    assert calendar_jd(2000, 1, 1.5) == 2451545.0
    assert calendar_jd(2025, 11, 21) == Time("2025-11-21").jd
    elements = [
        {"a": 2.36, "e": 0.09, "epoch_jd": 2461000.5, "mean_anomaly": 26.8,
         "arg_peri": 151.5, "node": 103.7, "inclination": 7.1, "h": 3.2, "g": 0.15},
        {"a": 1.46, "e": 0.22, "epoch_jd": 2461000.5, "mean_anomaly": 310.0,
         "arg_peri": 178.8, "node": 304.3, "inclination": 10.8, "h": 10.4, "g": 0.15},
    ]
    vectors, radii = asteroid_heliocentric_ecliptic_array(
        build_asteroid_arrays(elements), 2461150.25
    )
    for index, element in enumerate(elements):
        vector, radius = asteroid_heliocentric_ecliptic(element, 2461150.25)
        assert np.allclose(vectors[index], vector, atol=1e-10)
        assert np.isclose(radii[index], radius)
