import numpy as np

from sky_geometry import (
    append_sky_segment,
    inverse_global,
    inverse_local,
    milky_way_outer_envelope,
    project_global,
    project_local,
)


def circular_difference(first, second):
    return (np.asarray(first) - np.asarray(second) + 180.0) % 360.0 - 180.0


def test_global_projection_round_trips_supported_views():
    ra = np.array([25.0, 145.0, 275.0])
    dec = np.array([-35.0, 5.0, 42.0])
    for projection in range(5):
        x, y = project_global(projection, ra, dec, 160.0, 12.0)
        restored_ra, restored_dec = inverse_global(
            projection, x, y, 160.0, 12.0
        )
        assert np.allclose(circular_difference(restored_ra, ra), 0.0, atol=1e-8)
        assert np.allclose(restored_dec, dec, atol=1e-8)


def test_local_projection_round_trips_supported_views():
    azimuth = np.array([35.0, 155.0, 285.0])
    altitude = np.array([-20.0, 15.0, 55.0])
    for projection in range(5):
        x, y = project_local(projection, azimuth, altitude, 190.0, 25.0)
        restored_azimuth, restored_altitude = inverse_local(
            projection, x, y, 190.0, 25.0
        )
        expected_azimuth = (
            (azimuth - 180.0) % 360.0 if projection == 1 else azimuth
        )
        assert np.allclose(
            circular_difference(restored_azimuth, expected_azimuth),
            0.0,
            atol=1e-8,
        )
        assert np.allclose(restored_altitude, altitude, atol=1e-8)


def test_wrapped_segment_is_split_at_both_map_edges():
    x_values, y_values = [], []
    append_sky_segment(x_values, y_values, 359.0, 2.0, 1.0, 4.0)
    assert np.allclose(
        np.asarray(x_values)[np.isfinite(x_values)], [359.0, 360.0, 0.0, 1.0]
    )
    assert np.allclose(
        np.asarray(y_values)[np.isfinite(y_values)], [2.0, 3.0, 3.0, 4.0]
    )


def test_milky_way_envelope_keeps_only_outer_edges():
    right_ascension = np.linspace(0.0, 360.0, 25)
    lower = [(value, -12.0) for value in right_ascension]
    upper = [(value, 18.0) for value in right_ascension]
    envelope = milky_way_outer_envelope([[lower, upper]], sample_count=73)
    assert envelope is not None
    sample_ra, lower_dec, upper_dec = envelope
    assert len(sample_ra) == 73
    assert np.allclose(lower_dec, -12.0)
    assert np.allclose(upper_dec, 18.0)
