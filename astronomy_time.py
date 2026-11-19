"""Time-selection helpers for observing-night shortcuts."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.coordinates import AltAz, get_sun
from astropy.time import Time
from astropy.utils import iers


iers.conf.auto_download = False
iers.conf.auto_max_age = None
iers.conf.iers_degraded_accuracy = "ignore"


def next_twilight_crossing(location, target_altitude_deg, evening,
                           now_utc=None):
    """Return the nearest future evening or morning solar-altitude crossing.

    Searching forward from the current instant avoids the ambiguous civil-date
    boundary after midnight: at 02:55 the upcoming dawn still belongs to the
    current morning, rather than tomorrow morning.
    """
    current = now_utc or datetime.now(timezone.utc)
    if current.tzinfo is None:
        current = current.replace(tzinfo=timezone.utc)
    else:
        current = current.astimezone(timezone.utc)
    search_start = current.replace(second=0, microsecond=0)
    times = [
        search_start + timedelta(minutes=minute)
        for minute in range(36 * 60 + 1)
    ]
    astropy_times = Time(times)
    suns = get_sun(astropy_times).transform_to(
        AltAz(obstime=astropy_times, location=location)
    )
    altitudes = suns.alt.deg
    if evening:
        crossings = np.where(
            (altitudes[:-1] > target_altitude_deg)
            & (altitudes[1:] <= target_altitude_deg)
        )[0]
    else:
        crossings = np.where(
            (altitudes[:-1] < target_altitude_deg)
            & (altitudes[1:] >= target_altitude_deg)
        )[0]

    if len(crossings):
        index = int(crossings[0])
        if (abs(altitudes[index + 1] - target_altitude_deg)
                < abs(altitudes[index] - target_altitude_deg)):
            index += 1
    else:
        # Polar regions may not reach the requested twilight altitude. Keep
        # the former graceful behaviour and choose the closest future sample.
        index = int(np.argmin(np.abs(altitudes - target_altitude_deg)))
    return times[index], float(suns.az.deg[index])
