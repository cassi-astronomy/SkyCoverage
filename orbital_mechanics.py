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
