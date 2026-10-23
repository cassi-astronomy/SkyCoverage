"""Persistent equipment profiles, observing plans and script generation."""

from __future__ import annotations

import json
import math
import os
import re
from copy import deepcopy
from datetime import datetime, timedelta, timezone

import astropy.units as u
from astropy.coordinates import SkyCoord


DEFAULT_DATA = {
    "profiles": [
        {
            "name": "Výchozí sestava",
            "width_deg": 10.0,
            "height_deg": 7.0,
            "rotation_deg": 0.0,
            "dither_deg": 0.05,
            "readout_seconds": 0.0,
            "slew_seconds": 0.0,
        }
    ],
    "archives": [],
}


class PlannerStore:
    def __init__(self, path):
        self.path = path
        self.data = deepcopy(DEFAULT_DATA)
        self.load()

    def load(self):
        try:
            with open(self.path, "r", encoding="utf-8") as source:
                loaded = json.load(source)
            if isinstance(loaded.get("profiles"), list) and loaded["profiles"]:
                self.data["profiles"] = loaded["profiles"]
                for profile in self.data["profiles"]:
                    profile.setdefault("readout_seconds", 0.0)
                    profile.setdefault("slew_seconds", 0.0)
            if isinstance(loaded.get("archives"), list):
                self.data["archives"] = loaded["archives"]
        except (OSError, ValueError, TypeError, AttributeError):
            pass

    def save(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        temporary_path = self.path + ".tmp"
        with open(temporary_path, "w", encoding="utf-8") as output:
            json.dump(self.data, output, ensure_ascii=False, indent=2)
        os.replace(temporary_path, self.path)


def offset_coordinate(ra_deg, dec_deg, east_deg, north_deg):
    """Offset an ICRS coordinate in its local tangent plane."""
    origin = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    distance = math.hypot(east_deg, north_deg)
    if distance == 0:
        return float(origin.ra.deg), float(origin.dec.deg)
    position_angle = math.atan2(east_deg, north_deg) * u.rad
    target = origin.directional_offset_by(position_angle, distance * u.deg)
    return float(target.ra.deg % 360.0), float(target.dec.deg)


def rotated_offset(x_deg, y_deg, rotation_deg):
    angle = math.radians(rotation_deg)
    return (
        x_deg * math.cos(angle) - y_deg * math.sin(angle),
        x_deg * math.sin(angle) + y_deg * math.cos(angle),
    )


def safe_horizon_grid_steps(width_deg, height_deg, camera_x_angle_deg,
                            overlap_percent=10.0):
    """Return safe Alt/Az grid spacing for a tilted EQ-mounted sensor."""
    width = max(1e-9, float(width_deg))
    height = max(1e-9, float(height_deg))
    factor = max(0.01, 1.0 - float(overlap_percent) / 100.0)
    angle = math.radians(float(camera_x_angle_deg))

    def limited_step(first_size, first_component, second_size, second_component):
        limits = []
        if abs(first_component) > 1e-9:
            limits.append(first_size * factor / abs(first_component))
        if abs(second_component) > 1e-9:
            limits.append(second_size * factor / abs(second_component))
        return min(limits)

    horizontal = limited_step(width, math.cos(angle), height, math.sin(angle))
    vertical = limited_step(width, math.sin(angle), height, math.cos(angle))
    return horizontal, vertical


def suggested_archive_name(observation_datetime, profile_names, is_morning):
    """Build a stable night/profile/half-of-night archive name."""
    night_date = observation_datetime.date()
    if is_morning:
        night_date -= timedelta(days=1)
    clean_names = []
    for raw_name in profile_names:
        cleaned = re.sub(r"\s+", "-", str(raw_name).strip())
        cleaned = cleaned.replace("/", "-").replace("\\", "-") or "Sestava"
        if cleaned not in clean_names:
            clean_names.append(cleaned)
    equipment = "+".join(clean_names) if clean_names else "Sestava"
    period = "Morning" if is_morning else "Evening"
    return f"{night_date:%Y-%m-%d}_{equipment}_{period}"


def safe_plan_filename(name):
    """Return a Windows-safe filename stem while retaining the archive name."""
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "-", str(name).strip())
    cleaned = cleaned.rstrip(" .")
    return cleaned or "Pozorovani"


def estimate_sequence_duration(fields, exposure_seconds, loops):
    """Return exposure, readout, slew and total seconds for a field sequence.

    The initial slew is unknown and therefore excluded.  With multiple fields,
    every later field change is counted, including last-to-first between loops.
    """
    field_count = len(fields)
    loop_count = max(0, int(loops))
    exposure = max(0.0, float(exposure_seconds))
    exposure_total = field_count * loop_count * exposure
    readout_total = loop_count * sum(
        max(0.0, float(field.get("profile", {}).get("readout_seconds", 0.0)))
        for field in fields
    )

    slew_total = 0.0
    if field_count > 1 and loop_count > 0:
        # Within each loop the destination fields are 2..N.
        slew_total += loop_count * sum(
            max(0.0, float(field.get("profile", {}).get("slew_seconds", 0.0)))
            for field in fields[1:]
        )
        # A repeated loop moves from the final field back to the first one.
        first_slew = max(
            0.0,
            float(fields[0].get("profile", {}).get("slew_seconds", 0.0)),
        )
        slew_total += max(0, loop_count - 1) * first_slew

    return {
        "exposure_seconds": exposure_total,
        "readout_seconds": readout_total,
        "slew_seconds": slew_total,
        "total_seconds": exposure_total + readout_total + slew_total,
        "frame_count": field_count * loop_count,
    }


def format_duration(seconds):
    """Format an estimated duration without hiding short overheads."""
    total = max(0, int(round(float(seconds))))
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    parts = []
    if hours:
        parts.append(f"{hours} h")
    if minutes or hours:
        parts.append(f"{minutes} min")
    parts.append(f"{seconds} s")
    return " ".join(parts)


def export_plan_text(directory, archive, observation_time_utc="", observatory=""):
    """Atomically export one archived observing plan to one UTF-8 text file."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, safe_plan_filename(archive["name"]) + ".txt")
    if os.path.exists(path):
        raise FileExistsError(path)

    fields = archive.get("fields", [])
    profile_names = []
    for field in fields:
        name = str(field.get("profile", {}).get("name", "Sestava"))
        if name not in profile_names:
            profile_names.append(name)
    duration = estimate_sequence_duration(
        fields,
        archive.get("exposure_seconds", 0),
        archive.get("loops", 0),
    )
    lines = [
        f"Název pozorování: {archive['name']}",
        f"Čas mapy (UTC): {observation_time_utc}",
        f"Observatoř: {observatory}",
        f"Sestava: {', '.join(profile_names)}",
        f"Filtr: {archive.get('filter', '')}",
        f"Expozice: {float(archive.get('exposure_seconds', 0)):g} s",
        f"Cykly: {int(archive.get('loops', 0))}",
        f"Počet polí: {len(fields)}",
        f"Počet snímků: {duration['frame_count']}",
        f"Odhad délky sekvence: {format_duration(duration['total_seconds'])}",
        f"  expozice: {format_duration(duration['exposure_seconds'])}",
        f"  vyčítání: {format_duration(duration['readout_seconds'])}",
        f"  přejezdy mezi poli: {format_duration(duration['slew_seconds'])}",
        "",
        "Pole:",
        "#\tRA (J2000)\tDec (J2000)\tSestava",
    ]
    for index, field in enumerate(fields, 1):
        ra_text, dec_text = format_ra_dec(field["ra_deg"], field["dec_deg"])
        profile_name = field.get("profile", {}).get("name", "Sestava")
        lines.append(f"{index}\t{ra_text}\t{dec_text}\t{profile_name}")
    lines.extend([
        "",
        "Pozorovací skript:",
        generate_exposure_script(
            fields, archive.get("filter", ""),
            archive.get("exposure_seconds", 0), archive.get("loops", 0),
        ),
        "",
    ])

    temporary_path = path + ".tmp"
    try:
        with open(temporary_path, "w", encoding="utf-8", newline="\n") as output:
            output.write("\n".join(lines))
        os.replace(temporary_path, path)
    finally:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)
    return path


def make_mosaic(ra_deg, dec_deg, profile, rows=1, columns=1, overlap_percent=10.0):
    """Return field centres for a centred, camera-rotated mosaic."""
    width = float(profile["width_deg"])
    height = float(profile["height_deg"])
    rotation = float(profile["rotation_deg"])
    factor = max(0.01, 1.0 - float(overlap_percent) / 100.0)
    result = []
    for row in range(rows):
        for column in range(columns):
            x = (column - (columns - 1) / 2.0) * width * factor
            y = ((rows - 1) / 2.0 - row) * height * factor
            east, north = rotated_offset(x, y, rotation)
            field_ra, field_dec = offset_coordinate(ra_deg, dec_deg, east, north)
            result.append({
                "ra_deg": field_ra,
                "dec_deg": field_dec,
                "profile": deepcopy(profile),
            })
    return result


def field_outline(field, samples_per_edge=8):
    """Return sampled ICRS vertices for a rotated rectangular field."""
    profile = field["profile"]
    half_width = float(profile["width_deg"]) / 2.0
    half_height = float(profile["height_deg"]) / 2.0
    rotation = float(profile["rotation_deg"])
    corners = [
        (-half_width, -half_height),
        (-half_width, half_height),
        (half_width, half_height),
        (half_width, -half_height),
        (-half_width, -half_height),
    ]
    vertices = []
    for start, end in zip(corners[:-1], corners[1:]):
        for index in range(samples_per_edge):
            fraction = index / float(samples_per_edge)
            x = start[0] + fraction * (end[0] - start[0])
            y = start[1] + fraction * (end[1] - start[1])
            east, north = rotated_offset(x, y, rotation)
            vertices.append(offset_coordinate(
                field["ra_deg"], field["dec_deg"], east, north
            ))
    east, north = rotated_offset(corners[-1][0], corners[-1][1], rotation)
    vertices.append(offset_coordinate(field["ra_deg"], field["dec_deg"], east, north))
    return vertices


def format_ra_dec(ra_deg, dec_deg):
    coordinate = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    ra = coordinate.ra.hms
    dec = coordinate.dec.signed_dms
    ra_text = f"{int(ra.h):02d} {int(ra.m):02d} {ra.s:05.2f}"
    sign = "+" if dec.sign >= 0 else "-"
    dec_text = f"{sign}{int(abs(dec.d)):02d} {int(abs(dec.m)):02d} {abs(dec.s):04.1f}"
    return ra_text, dec_text


def generate_exposure_script(fields, filter_name, exposure_seconds, loops):
    """Generate the telescope command used by the original off_CO.py script."""
    if not fields:
        return "Kliknutím do mapy přidejte alespoň jedno pozorovací pole."
    first = fields[0]
    ra_text, dec_text = format_ra_dec(first["ra_deg"], first["dec_deg"])
    commands = []
    for field in fields:
        raw_delta_ra = (float(field["ra_deg"]) - float(first["ra_deg"]) + 180.0) % 360.0 - 180.0
        delta_dec = float(field["dec_deg"]) - float(first["dec_deg"])
        dither = float(field["profile"].get("dither_deg", 0.05))
        if field is first:
            ra_low, ra_high = 0.0, dither
            dec_low, dec_high = 0.0, dither
        else:
            half = dither / 2.0
            ra_low, ra_high = raw_delta_ra - half, raw_delta_ra + half
            dec_low, dec_high = delta_dec - half, delta_dec + half
        commands.append(
            "for 1 { T0.WOFFS=(rand(%.3f,%.3f),rand(%.3f,%.3f)) E %g }"
            % (ra_low, ra_high, dec_low, dec_high, float(exposure_seconds))
        )
    filter_value = str(filter_name).strip() or "N"
    command = f"loopdisable filter={filter_value} for {int(loops)} {{ " + " ".join(commands) + " }"
    return f"# Výchozí souřadnice (RA Dec)\n{ra_text} {dec_text}\n\n{command}"


def new_archive(name, fields, filter_name, exposure_seconds, loops):
    return {
        "name": name,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "visible": True,
        "filter": filter_name,
        "exposure_seconds": float(exposure_seconds),
        "loops": int(loops),
        "fields": deepcopy(fields),
    }
