from datetime import datetime, timezone

from astropy.coordinates import EarthLocation

from astronomy_time import next_twilight_crossing


def test_morning_twilight_after_midnight_stays_on_current_utc_date():
    location = EarthLocation.from_geodetic(14.2858, 48.8639, 1070)
    now = datetime(2026, 10, 25, 0, 55, tzinfo=timezone.utc)
    twilight, _ = next_twilight_crossing(
        location, -18.0, evening=False, now_utc=now
    )
    assert twilight > now
    assert twilight.date() == now.date()


def test_morning_twilight_after_dawn_selects_the_next_day():
    location = EarthLocation.from_geodetic(14.2858, 48.8639, 1070)
    now = datetime(2026, 10, 25, 8, 0, tzinfo=timezone.utc)
    twilight, _ = next_twilight_crossing(
        location, -18.0, evening=False, now_utc=now
    )
    assert twilight > now
    assert twilight.date().isoformat() == "2026-10-26"
