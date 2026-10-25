"""Pure map-projection and sky-overlay geometry helpers.

This module deliberately contains no Qt objects.  Keeping the numerical work
outside the main window makes it independently testable and avoids coupling
catalog geometry to widget state.
"""

from __future__ import annotations

import numpy as np


LAMBERT_Y_SCALE = 180.0 / np.pi
STEREO_VIEW_RADIUS = 7.5


def append_sky_segment(x_values, y_values, x1, y1, x2, y2):
    """Append a line segment without drawing across the azimuth/RA wrap."""
    delta = x2 - x1
    if abs(delta) <= 180.0:
        x_values.extend([x1, x2, np.nan])
        y_values.extend([y1, y2, np.nan])
        return

    original_x2 = x2
    if delta > 0:
        x2 -= 360.0
        boundary = 0.0
        opposite_boundary = 360.0
    else:
        x2 += 360.0
        boundary = 360.0
        opposite_boundary = 0.0

    fraction = (boundary - x1) / (x2 - x1)
    boundary_y = y1 + fraction * (y2 - y1)
    x_values.extend([x1, boundary, np.nan, opposite_boundary, original_x2, np.nan])
    y_values.extend([y1, boundary_y, np.nan, boundary_y, y2, np.nan])


def append_projected_segment(
        x_values, y_values, x1, y1, x2, y2, max_step=0.35):
    """Append a short segment in a projected X/Y plane."""
    if abs(x2 - x1) <= max_step and abs(y2 - y1) <= max_step:
        x_values.extend([x1, x2, np.nan])
        y_values.extend([y1, y2, np.nan])


def project_oblique_mercator(longitude_deg, latitude_deg,
                              center_longitude_deg, center_latitude_deg):
    longitude = np.radians(np.asarray(longitude_deg, dtype=float))
    latitude = np.radians(np.asarray(latitude_deg, dtype=float))
    center_longitude = np.radians(float(center_longitude_deg))
    center_latitude = np.radians(float(center_latitude_deg))
    delta_longitude = (
        longitude - center_longitude + np.pi
    ) % (2 * np.pi) - np.pi
    rotated_x = np.cos(latitude) * np.sin(delta_longitude)
    rotated_y = (
        np.cos(center_latitude) * np.sin(latitude)
        - np.sin(center_latitude) * np.cos(latitude) * np.cos(delta_longitude)
    )
    rotated_z = (
        np.sin(center_latitude) * np.sin(latitude)
        + np.cos(center_latitude) * np.cos(latitude) * np.cos(delta_longitude)
    )
    rotated_longitude = np.arctan2(rotated_x, rotated_z)
    rotated_latitude = np.arcsin(np.clip(rotated_y, -1.0, 1.0))
    safe_latitude = np.clip(
        rotated_latitude, np.radians(-89.5), np.radians(89.5)
    )
    return (
        np.degrees(rotated_longitude),
        np.degrees(np.arcsinh(np.tan(safe_latitude))),
    )


def inverse_oblique_mercator(x, y, center_longitude_deg,
                              center_latitude_deg):
    rotated_longitude = np.radians(np.asarray(x, dtype=float))
    rotated_latitude = np.arctan(
        np.sinh(np.radians(np.asarray(y, dtype=float)))
    )
    center_longitude = np.radians(float(center_longitude_deg))
    center_latitude = np.radians(float(center_latitude_deg))
    east = np.cos(rotated_latitude) * np.sin(rotated_longitude)
    north = np.sin(rotated_latitude)
    forward = np.cos(rotated_latitude) * np.cos(rotated_longitude)
    latitude = np.arcsin(np.clip(
        np.cos(center_latitude) * north + np.sin(center_latitude) * forward,
        -1.0,
        1.0,
    ))
    longitude = center_longitude + np.arctan2(
        east,
        np.cos(center_latitude) * forward - np.sin(center_latitude) * north,
    )
    return np.degrees(longitude) % 360.0, np.degrees(latitude)


def project_global(projection, ra, dec, center_ra=180.0, center_dec=0.0):
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    if projection == 0:
        return ra % 360.0, dec
    if projection == 2:
        return ra % 360.0, LAMBERT_Y_SCALE * np.sin(np.radians(dec))
    if projection == 4:
        return project_oblique_mercator(ra, dec, center_ra, center_dec)

    longitude = np.radians((ra - 180.0 + 180.0) % 360.0 - 180.0)
    latitude = np.radians(dec)
    if projection == 3:
        denominator = 1.0 + np.cos(latitude) * np.cos(longitude)
        safe = denominator > 1e-6
        safe_denominator = np.where(safe, denominator, 1.0)
        x = np.where(
            safe,
            2.0 * np.cos(latitude) * np.sin(longitude) / safe_denominator,
            np.nan,
        )
        y = np.where(
            safe, 2.0 * np.sin(latitude) / safe_denominator, np.nan
        )
        visible = np.hypot(x, y) <= STEREO_VIEW_RADIUS * 1.02
        return np.where(visible, x, np.nan), np.where(visible, y, np.nan)

    denominator = np.sqrt(
        1.0 + np.cos(latitude) * np.cos(longitude / 2.0)
    )
    return (
        2.0 * np.sqrt(2.0) * np.cos(latitude)
        * np.sin(longitude / 2.0) / denominator,
        np.sqrt(2.0) * np.sin(latitude) / denominator,
    )


def inverse_global(projection, x, y, center_ra=180.0, center_dec=0.0):
    if projection == 0:
        return np.asarray(x), np.asarray(y)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if projection == 2:
        valid = np.abs(y) <= LAMBERT_Y_SCALE
        dec = np.where(
            valid,
            np.degrees(np.arcsin(np.clip(
                y / LAMBERT_Y_SCALE, -1.0, 1.0
            ))),
            np.nan,
        )
        return x % 360.0, dec
    if projection == 4:
        return inverse_oblique_mercator(x, y, center_ra, center_dec)
    if projection == 3:
        radius = np.hypot(x, y)
        angular_distance = 2.0 * np.arctan2(radius, 2.0)
        radius_safe = np.where(radius > 1e-12, radius, 1.0)
        dec = np.arcsin(np.clip(
            y * np.sin(angular_distance) / radius_safe, -1.0, 1.0
        ))
        longitude = np.arctan2(
            x * np.sin(angular_distance),
            radius_safe * np.cos(angular_distance),
        )
        dec = np.where(radius > 1e-12, dec, 0.0)
        longitude = np.where(radius > 1e-12, longitude, 0.0)
        valid = radius <= STEREO_VIEW_RADIUS * 1.02
        return (
            np.where(
                valid, (np.degrees(longitude) + 180.0) % 360.0, np.nan
            ),
            np.where(valid, np.degrees(dec), np.nan),
        )
    z = np.sqrt(np.maximum(
        0.0, 1.0 - (x ** 2 / 16.0) - (y ** 2 / 4.0)
    ))
    longitude = 2.0 * np.arctan2(
        z * x, 2.0 * (2.0 * z ** 2 - 1.0)
    )
    latitude = np.arcsin(np.clip(z * y, -1.0, 1.0))
    return (np.degrees(longitude) + 180.0) % 360.0, np.degrees(latitude)


def project_local(projection, azimuth, altitude, center_azimuth=180.0,
                  center_altitude=0.0):
    azimuth = np.asarray(azimuth, dtype=float)
    altitude = np.asarray(altitude, dtype=float)
    if projection == 0:
        return azimuth, altitude
    if projection == 2:
        return azimuth, LAMBERT_Y_SCALE * np.sin(np.radians(altitude))
    if projection == 4:
        return project_oblique_mercator(
            azimuth, altitude, center_azimuth, center_altitude
        )
    if projection == 3:
        angle = np.radians((azimuth - 180.0 + 180.0) % 360.0 - 180.0)
        radius = 2.0 * np.tan(np.radians((90.0 - altitude) / 2.0))
        return radius * np.sin(angle), -radius * np.cos(angle)
    longitude = np.radians(
        (azimuth - 180.0 + 180.0) % 360.0 - 180.0
    )
    latitude = np.radians(altitude)
    denominator = np.sqrt(
        1.0 + np.cos(latitude) * np.cos(longitude / 2.0)
    )
    return (
        2.0 * np.sqrt(2.0) * np.cos(latitude)
        * np.sin(longitude / 2.0) / denominator,
        np.sqrt(2.0) * np.sin(latitude) / denominator,
    )


def inverse_local(projection, x, y, center_azimuth=180.0,
                  center_altitude=0.0):
    if projection == 0:
        return np.asarray(x), np.asarray(y)
    if projection == 2:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        valid = np.abs(y) <= LAMBERT_Y_SCALE
        altitude = np.where(
            valid,
            np.degrees(np.arcsin(np.clip(
                y / LAMBERT_Y_SCALE, -1.0, 1.0
            ))),
            np.nan,
        )
        return x, altitude
    if projection == 4:
        return inverse_oblique_mercator(
            x, y, center_azimuth, center_altitude
        )

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if projection == 3:
        radius = np.hypot(x, y)
        altitude = 90.0 - np.degrees(2.0 * np.arctan2(radius, 2.0))
        angle = np.degrees(np.arctan2(x, -y))
        azimuth = (angle + 180.0) % 360.0
        azimuth = np.where(radius > 1e-12, azimuth, 180.0)
        return azimuth, altitude

    z = np.sqrt(np.maximum(
        0.0, 1.0 - (x ** 2 / 16.0) - (y ** 2 / 4.0)
    ))
    longitude = 2.0 * np.arctan2(
        z * x, 2.0 * (2.0 * z ** 2 - 1.0)
    )
    latitude = np.arcsin(np.clip(z * y, -1.0, 1.0))
    return np.degrees(longitude), np.degrees(latitude)


def milky_way_outer_envelope(contours, sample_count=721):
    """Return the outer RA/Dec envelope of the faintest Milky Way level.

    The source's two main rings can fold back in right ascension.  Sampling
    every segment and retaining the extrema produces only the visible outside
    edges and intentionally discards internal structure.
    """
    if not contours or not contours[0] or len(contours[0]) < 2:
        return None
    rings = contours[0][:2]
    sample_ra = np.linspace(0.0, 360.0, sample_count)
    lower_dec = np.full(sample_ra.shape, np.inf)
    upper_dec = np.full(sample_ra.shape, -np.inf)
    for ring in rings:
        polygon = np.asarray(ring, dtype=float)
        for index in range(len(polygon) - 1):
            first_ra = polygon[index, 0] % 360.0
            delta_ra = (
                polygon[index + 1, 0] - first_ra + 180.0
            ) % 360.0 - 180.0
            if abs(delta_ra) < 1e-9:
                continue
            aligned_ra = (
                (sample_ra - first_ra + 180.0) % 360.0 - 180.0 + first_ra
            )
            second_ra = first_ra + delta_ra
            crosses = (
                (aligned_ra - first_ra) * (aligned_ra - second_ra) <= 1e-9
            )
            if not np.any(crosses):
                continue
            fraction = (aligned_ra[crosses] - first_ra) / delta_ra
            crossing_dec = (
                polygon[index, 1]
                + fraction * (polygon[index + 1, 1] - polygon[index, 1])
            )
            lower_dec[crosses] = np.minimum(lower_dec[crosses], crossing_dec)
            upper_dec[crosses] = np.maximum(upper_dec[crosses], crossing_dec)
    if not np.all(np.isfinite(lower_dec) & np.isfinite(upper_dec)):
        return None
    return sample_ra, lower_dec, upper_dec
