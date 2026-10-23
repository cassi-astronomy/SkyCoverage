from datetime import datetime

import numpy as np
from astropy.time import Time

from coverage import (
    MainWindow, PLANET_STYLES, PointingDataReader, SupplementCatalogReader,
    format_sun_altitude,
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


def test_planet_layer_uses_short_labels_and_omits_pluto():
    assert list(PLANET_STYLES) == [
        "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune",
    ]
    assert [style[0] for style in PLANET_STYLES.values()] == [
        "Me", "Ve", "Ma", "Ju", "Sa", "Ur", "Ne",
    ]
