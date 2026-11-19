"""Pure two-body orbit propagation used by the sky-map overlays."""

from __future__ import annotations

import math

import numpy as np


GAUSSIAN_GRAVITATIONAL_CONSTANT = 0.01720209895


def comet_heliocentric_ecliptic(element, jd):
    """Approximate a comet's heliocentric J2000 ecliptic vector in AU."""
    q = element["q"]
    eccentricity = element["e"]
    elapsed_days = jd - element["tp_jd"]
    if q <= 0.0 or eccentricity < 0.0:
        return None

    if abs(eccentricity - 1.0) < 1e-4:
        barker = (
            GAUSSIAN_GRAVITATIONAL_CONSTANT * elapsed_days
            / np.sqrt(2.0 * q ** 3)
        )
        anomaly = 2.0 * np.sinh(np.arcsinh(1.5 * barker) / 3.0)
        true_anomaly = 2.0 * np.arctan(anomaly)
        radius = q * (1.0 + anomaly ** 2)
    elif eccentricity < 1.0:
        semi_major = q / (1.0 - eccentricity)
        mean_anomaly = (
            GAUSSIAN_GRAVITATIONAL_CONSTANT * elapsed_days
            / semi_major ** 1.5
        )
        mean_anomaly = (mean_anomaly + np.pi) % (2.0 * np.pi) - np.pi
        eccentric_anomaly = mean_anomaly
        if eccentricity > 0.8:
            eccentric_anomaly = np.pi if mean_anomaly >= 0.0 else -np.pi
        for _ in range(30):
            correction = (
                eccentric_anomaly
                - eccentricity * np.sin(eccentric_anomaly)
                - mean_anomaly
            ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
            eccentric_anomaly -= correction
            if abs(correction) < 1e-12:
                break
        radius = semi_major * (
            1.0 - eccentricity * np.cos(eccentric_anomaly)
        )
        true_anomaly = 2.0 * np.arctan2(
            np.sqrt(1.0 + eccentricity)
            * np.sin(eccentric_anomaly / 2.0),
            np.sqrt(1.0 - eccentricity)
            * np.cos(eccentric_anomaly / 2.0),
        )
    else:
        semi_major = q / (eccentricity - 1.0)
        mean_anomaly = (
            GAUSSIAN_GRAVITATIONAL_CONSTANT * elapsed_days
            / semi_major ** 1.5
        )
        hyperbolic_anomaly = np.arcsinh(mean_anomaly / eccentricity)
        for _ in range(30):
            correction = (
                eccentricity * np.sinh(hyperbolic_anomaly)
                - hyperbolic_anomaly - mean_anomaly
            ) / (eccentricity * np.cosh(hyperbolic_anomaly) - 1.0)
            hyperbolic_anomaly -= correction
            if abs(correction) < 1e-12:
                break
        radius = semi_major * (
            eccentricity * np.cosh(hyperbolic_anomaly) - 1.0
        )
        true_anomaly = 2.0 * np.arctan2(
            np.sqrt(eccentricity + 1.0)
            * np.sinh(hyperbolic_anomaly / 2.0),
            np.sqrt(eccentricity - 1.0)
            * np.cosh(hyperbolic_anomaly / 2.0),
        )

    argument = np.radians(element["arg_peri"]) + true_anomaly
    node = np.radians(element["node"])
    inclination = np.radians(element["inclination"])
    x = radius * (
        np.cos(node) * np.cos(argument)
        - np.sin(node) * np.sin(argument) * np.cos(inclination)
    )
    y = radius * (
        np.sin(node) * np.cos(argument)
        + np.cos(node) * np.sin(argument) * np.cos(inclination)
    )
    z = radius * np.sin(argument) * np.sin(inclination)
    return np.array([x, y, z], dtype=float), float(radius)


def asteroid_heliocentric_ecliptic(element, jd):
    """Propagate an elliptic MPC orbit and return its J2000 vector."""
    semi_major = float(element["a"])
    eccentricity = float(element["e"])
    if semi_major <= 0.0 or not 0.0 <= eccentricity < 1.0:
        return None
    mean_anomaly = math.radians(element["mean_anomaly"])
    mean_anomaly += (
        GAUSSIAN_GRAVITATIONAL_CONSTANT / semi_major ** 1.5
        * (float(jd) - float(element["epoch_jd"]))
    )
    mean_anomaly = (mean_anomaly + math.pi) % (2.0 * math.pi) - math.pi
    eccentric_anomaly = mean_anomaly
    for _ in range(20):
        correction = (
            eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)
            - mean_anomaly
        ) / (1.0 - eccentricity * math.cos(eccentric_anomaly))
        eccentric_anomaly -= correction
        if abs(correction) < 1e-12:
            break
    x_orbit = semi_major * (math.cos(eccentric_anomaly) - eccentricity)
    y_orbit = (
        semi_major * math.sqrt(1.0 - eccentricity ** 2)
        * math.sin(eccentric_anomaly)
    )
    radius = semi_major * (
        1.0 - eccentricity * math.cos(eccentric_anomaly)
    )
    node, inclination, argument = np.radians([
        element["node"], element["inclination"], element["arg_peri"]
    ])
    cos_node, sin_node = math.cos(node), math.sin(node)
    cos_inc, sin_inc = math.cos(inclination), math.sin(inclination)
    cos_arg, sin_arg = math.cos(argument), math.sin(argument)
    x = (
        (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit
        + (-cos_node * sin_arg - sin_node * cos_arg * cos_inc) * y_orbit
    )
    y = (
        (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit
        + (-sin_node * sin_arg + cos_node * cos_arg * cos_inc) * y_orbit
    )
    z = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit
    return np.array([x, y, z], dtype=float), float(radius)


EARTH_PERIHELION_AU = 0.9833
EARTH_APHELION_AU = 1.0167


def calendar_jd(year, month, day):
    """Julian date at 0h of a Gregorian calendar date (fractional day allowed)."""
    year = int(year)
    month = int(month)
    if month <= 2:
        year -= 1
        month += 12
    century = year // 100
    correction = 2 - century + century // 4
    return (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + float(day) + correction - 1524.5
    )


def asteroid_brightest_possible_magnitude(h, q, closest_earth_au=0.2):
    """Optimistic opposition magnitude an orbit with perihelion ``q`` can reach.

    Uses the smallest heliocentric and geocentric distance the orbit allows and
    no phase darkening, so the result is a lower bound of V. Near-Earth orbits
    are limited to ``closest_earth_au`` because two-body propagation is not
    reliable for very close approaches anyway.
    """
    h = np.asarray(h, dtype=float)
    q = np.asarray(q, dtype=float)
    radius = np.maximum(q, EARTH_PERIHELION_AU)
    distance = np.maximum(q - EARTH_APHELION_AU, closest_earth_au)
    return h + 5.0 * np.log10(radius * distance)


def build_asteroid_arrays(elements):
    """Convert parsed MPC asteroid elements into numpy arrays once per load."""
    def column(key, default=0.0):
        return np.asarray(
            [float(item.get(key, default)) for item in elements], dtype=float
        )

    if not elements:
        return None
    node, inclination, argument = (
        np.radians(column("node")),
        np.radians(column("inclination")),
        np.radians(column("arg_peri")),
    )
    cos_node, sin_node = np.cos(node), np.sin(node)
    cos_inc, sin_inc = np.cos(inclination), np.sin(inclination)
    cos_arg, sin_arg = np.cos(argument), np.sin(argument)
    semi_major = column("a")
    return {
        "a": semi_major,
        "e": column("e"),
        "h": column("h"),
        "g": column("g", 0.15),
        "epoch_jd": column("epoch_jd"),
        "mean_anomaly": np.radians(column("mean_anomaly")),
        "mean_motion": GAUSSIAN_GRAVITATIONAL_CONSTANT / semi_major ** 1.5,
        # Rotation from the orbital plane to J2000 ecliptic, columns P and Q.
        "p": np.stack([
            cos_node * cos_arg - sin_node * sin_arg * cos_inc,
            sin_node * cos_arg + cos_node * sin_arg * cos_inc,
            sin_arg * sin_inc,
        ], axis=1),
        "q": np.stack([
            -cos_node * sin_arg - sin_node * cos_arg * cos_inc,
            -sin_node * sin_arg + cos_node * cos_arg * cos_inc,
            cos_arg * sin_inc,
        ], axis=1),
    }


def asteroid_heliocentric_ecliptic_array(arrays, jd):
    """Vectorised elliptic propagation; returns (N×3 vectors, radii) in AU."""
    eccentricity = arrays["e"]
    mean_anomaly = arrays["mean_anomaly"] + arrays["mean_motion"] * (
        float(jd) - arrays["epoch_jd"]
    )
    mean_anomaly = (mean_anomaly + np.pi) % (2.0 * np.pi) - np.pi
    eccentric_anomaly = np.where(
        eccentricity > 0.8, np.pi * np.sign(mean_anomaly), mean_anomaly
    )
    for _ in range(30):
        correction = (
            eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)
            - mean_anomaly
        ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
        eccentric_anomaly = eccentric_anomaly - correction
        if np.max(np.abs(correction), initial=0.0) < 1e-12:
            break
    semi_major = arrays["a"]
    x_orbit = semi_major * (np.cos(eccentric_anomaly) - eccentricity)
    y_orbit = (
        semi_major * np.sqrt(1.0 - eccentricity ** 2)
        * np.sin(eccentric_anomaly)
    )
    radius = semi_major * (1.0 - eccentricity * np.cos(eccentric_anomaly))
    vectors = arrays["p"] * x_orbit[:, None] + arrays["q"] * y_orbit[:, None]
    return vectors, radius


def asteroid_apparent_positions(arrays, jd, earth_ecliptic, obliquity_deg=23.4392911):
    """Return geocentric RA, Dec (deg) and V magnitude (H, G system) arrays."""
    if arrays is None or len(arrays["a"]) == 0:
        empty = np.zeros(0, dtype=float)
        return empty, empty, empty
    vectors, radius = asteroid_heliocentric_ecliptic_array(arrays, jd)
    geocentric = vectors - np.asarray(earth_ecliptic, dtype=float)[None, :]
    distance = np.linalg.norm(geocentric, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_phase = np.sum(-vectors * -geocentric, axis=1) / (radius * distance)
        phase_angle = np.arccos(np.clip(cos_phase, -1.0, 1.0))
        tangent = np.tan(phase_angle / 2.0)
        phi1 = np.exp(-3.33 * tangent ** 0.63)
        phi2 = np.exp(-1.87 * tangent ** 1.22)
        g = arrays["g"]
        phase_function = np.maximum((1.0 - g) * phi1 + g * phi2, 1e-12)
        magnitude = (
            arrays["h"] + 5.0 * np.log10(radius * distance)
            - 2.5 * np.log10(phase_function)
        )
    obliquity = math.radians(obliquity_deg)
    x = geocentric[:, 0]
    y = geocentric[:, 1] * math.cos(obliquity) - geocentric[:, 2] * math.sin(obliquity)
    z = geocentric[:, 1] * math.sin(obliquity) + geocentric[:, 2] * math.cos(obliquity)
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    with np.errstate(divide="ignore", invalid="ignore"):
        dec = np.degrees(np.arcsin(np.clip(z / distance, -1.0, 1.0)))
    magnitude = np.where((distance > 0.0) & (radius > 0.0), magnitude, np.nan)
    return ra, dec, magnitude
