"""Tests for the fast map helpers that replace per-point astropy calls."""

from datetime import datetime, timezone

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers
from PyQt6.QtCore import QDateTime, QTimeZone

from coverage import (
    MainWindow, mpc_file_year_is_plausible, qdatetime_to_time,
    qdatetime_to_utc,
)
from observation_planner import (
    field_outline, field_outline_arrays, offset_coordinate, offset_coordinates,
)
from sky_geometry import (
    angular_separation, append_projected_segment, append_sky_segment,
    fit_rotation, polyline_segments, projected_segments, rotate_angles,
    sky_segments, unit_vectors,
)

iers.conf.auto_download = False


def circular_difference(first, second):
    return (np.asarray(first) - np.asarray(second) + 180.0) % 360.0 - 180.0


def test_vectorised_offset_matches_astropy_directional_offset():
    rng = np.random.default_rng(7)
    ra = rng.uniform(0.0, 360.0, 50)
    dec = rng.uniform(-88.0, 88.0, 50)
    east = rng.uniform(-12.0, 12.0, 50)
    north = rng.uniform(-12.0, 12.0, 50)
    result_ra, result_dec = offset_coordinates(ra, dec, east, north)
    expected = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).directional_offset_by(
        np.arctan2(east, north) * u.rad, np.hypot(east, north) * u.deg
    )
    assert np.allclose(circular_difference(result_ra, expected.ra.deg), 0.0, atol=1e-9)
    assert np.allclose(result_dec, expected.dec.deg, atol=1e-9)
    assert offset_coordinate(10.0, 20.0, 0.0, 0.0) == (10.0, 20.0)


def test_field_outline_arrays_match_list_outline():
    field = {
        "ra_deg": 359.6, "dec_deg": 71.0,
        "profile": {"width_deg": 10.0, "height_deg": 7.0, "rotation_deg": 27.0},
    }
    ra, dec = field_outline_arrays(field)
    vertices = field_outline(field)
    assert len(vertices) == len(ra) == 33
    assert np.allclose(ra, [point[0] for point in vertices])
    assert np.allclose(dec, [point[1] for point in vertices])
    assert np.isclose(ra[0], ra[-1]) and np.isclose(dec[0], dec[-1])


def test_vectorised_segments_match_scalar_helpers():
    rng = np.random.default_rng(3)
    x = rng.uniform(0.0, 360.0, 200)
    y = rng.uniform(-90.0, 90.0, 200)
    expected_x, expected_y = [], []
    for index in range(len(x) - 1):
        append_sky_segment(expected_x, expected_y, x[index], y[index], x[index + 1], y[index + 1])
    result_x, result_y = sky_segments(x[:-1], y[:-1], x[1:], y[1:])
    assert np.allclose(result_x, expected_x, equal_nan=True)
    assert np.allclose(result_y, expected_y, equal_nan=True)

    small_x, small_y = x / 200.0, y / 200.0
    expected_x, expected_y = [], []
    for index in range(len(x) - 1):
        append_projected_segment(
            expected_x, expected_y, small_x[index], small_y[index],
            small_x[index + 1], small_y[index + 1], 0.35,
        )
    result_x, result_y = projected_segments(
        small_x[:-1], small_y[:-1], small_x[1:], small_y[1:], 0.35
    )
    assert np.allclose(result_x, expected_x, equal_nan=True)
    assert np.allclose(result_y, expected_y, equal_nan=True)


def test_polyline_segments_never_join_separate_rows():
    x = np.array([[10.0, 20.0], [300.0, 310.0]])
    y = np.array([[0.0, 1.0], [5.0, 6.0]])
    result_x, _ = polyline_segments(x, y, projected=False)
    finite = result_x[np.isfinite(result_x)]
    assert np.allclose(finite, [10.0, 20.0, 300.0, 310.0])
    assert np.isnan(result_x[2])


def test_fitted_horizontal_rotation_matches_astropy_altaz():
    location = EarthLocation.from_geodetic(-70.4042, -24.6272, 2635)
    frame = AltAz(obstime=Time("2026-10-29T23:30:00", scale="utc"), location=location)
    samples_ra = np.array([0, 60, 120, 180, 240, 300, 30, 150, 270, 90, 210, 330, 0, 0], float)
    samples_dec = np.array([0, 0, 0, 0, 0, 0, 55, 55, 55, -55, -55, -55, 89, -89], float)
    horizontal = SkyCoord(ra=samples_ra * u.deg, dec=samples_dec * u.deg).transform_to(frame)
    matrix = fit_rotation(
        unit_vectors(samples_ra, samples_dec),
        unit_vectors(horizontal.az.deg, horizontal.alt.deg),
    )

    rng = np.random.default_rng(11)
    ra = rng.uniform(0.0, 360.0, 300)
    dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, 300)))
    azimuth, altitude = rotate_angles(matrix, ra, dec)
    expected = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).transform_to(frame)
    error = angular_separation(azimuth, altitude, expected.az.deg, expected.alt.deg)
    assert np.max(error) < 0.01

    restored_ra, restored_dec = rotate_angles(matrix.T, azimuth, altitude)
    assert np.allclose(circular_difference(restored_ra, ra), 0.0, atol=1e-8)
    assert np.allclose(restored_dec, dec, atol=1e-8)


def test_angular_separation_matches_astropy():
    first = SkyCoord(ra=12.0 * u.deg, dec=-40.0 * u.deg)
    second = SkyCoord(ra=250.0 * u.deg, dec=33.0 * u.deg)
    assert np.isclose(
        angular_separation(12.0, -40.0, 250.0, 33.0), first.separation(second).deg
    )


def test_map_time_uses_real_utc_instant_in_any_display_zone():
    epoch = int(datetime(2026, 10, 30, 0, 0, tzinfo=timezone.utc).timestamp())
    for offset in (0, 3600, 7200):
        value = QDateTime.fromSecsSinceEpoch(epoch, QTimeZone.fromSecondsAheadOfUtc(offset))
        assert qdatetime_to_utc(value) == datetime(2026, 10, 30, 0, 0, tzinfo=timezone.utc)
        assert qdatetime_to_time(value).isot == "2026-10-30T00:00:00.000"


def test_mpc_file_year_limit_follows_current_year():
    now = datetime(2035, 6, 1, tzinfo=timezone.utc)
    assert mpc_file_year_is_plausible(2035, now)
    assert mpc_file_year_is_plausible(2036, now)
    assert not mpc_file_year_is_plausible(2037, now)
    assert not mpc_file_year_is_plausible(1989, now)


def test_catalog_geometry_is_prepared_once_as_arrays():
    segments = MainWindow.build_constellation_segments([[[370.0, 1.0], [20.0, 2.0], [30.0, 3.0]]])
    assert segments.shape == (2, 4)
    assert np.isclose(segments[0, 0], 10.0)
    coordinates, lengths = MainWindow.build_boundary_arrays([
        [[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0]], [[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]],
    ])
    assert lengths.tolist() == [2, 3]
    assert coordinates.shape == (5, 2)
    polygons = MainWindow.build_polygon_arrays([
        ([1.0, 2.0, 3.0, 4.0, 1.0], [5.0, 6.0, 7.0, 8.0, 5.0], 19.5),
    ])
    assert polygons["ra"].shape == (1, 5)
    assert polygons["mag"].tolist() == [19.5]
