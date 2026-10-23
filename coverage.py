import sys
import os
import json
import csv
import gzip
import shutil
import tarfile
import tempfile
import urllib.request
import math
from datetime import datetime, timedelta, timezone
import numpy as np

from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QSpinBox, QDoubleSpinBox, 
                             QPushButton, QComboBox, QDateTimeEdit, QCheckBox,
                             QFrame, QTabWidget, QLineEdit, QPlainTextEdit,
                             QTableWidget, QTableWidgetItem, QHeaderView,
                             QInputDialog, QMessageBox, QAbstractItemView,
                             QGridLayout, QGraphicsItem, QGraphicsPathItem,
                             QSplitter)
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QDateTime, QTimer, QSize, QRectF
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPainterPath, QPen
import pyqtgraph as pg

import astropy.units as u
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import (get_sun, get_body, EarthLocation, SkyCoord,
                                 AltAz, Galactic, GeocentricTrueEcliptic,
                                 HeliocentricMeanEcliptic, get_body_barycentric)

from observation_planner import (
    PlannerStore, estimate_sequence_duration, export_plan_text, field_outline,
    format_duration, generate_exposure_script, make_mosaic, new_archive,
    offset_coordinate, rotated_offset, safe_horizon_grid_steps,
    suggested_archive_name,
)


class MinuteStepDateTimeEdit(QDateTimeEdit):
    """Date/time editor whose arrow buttons always move by one minute."""

    def stepBy(self, steps):
        self.setDateTime(self.dateTime().addSecs(int(steps) * 60))
        self.setCurrentSection(QDateTimeEdit.Section.MinuteSection)


class AntialiasedPathItem(QGraphicsPathItem):
    """A filled vector path whose curved edge stays smooth in the sky map."""

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        super().paint(painter, option, widget)

# Lokální Alt/Az musí fungovat i bez sítě; přibalená IERS tabulka je pro
# plánovací mapu dostatečná i po skončení jejího predikčního intervalu.
iers.conf.auto_download = False
iers.conf.auto_max_age = None
iers.conf.iers_degraded_accuracy = "ignore"


OBSERVATORIES = {
    '048 - Kleť (CZ)': EarthLocation.from_geodetic(14.2858, 48.8639, 1070),
    '950 - La Palma (ES)': EarthLocation.from_geodetic(-17.8816, 28.7603, 2396),
    'I05 - Paranal (CL)': EarthLocation.from_geodetic(-70.4042, -24.6272, 2635),
    'I47 - Pierre Auger (AR)': EarthLocation.from_geodetic(-69.3158, -35.2060, 1400)
}

SKYCOV_RETENTION_DAYS = 92
SKYCOV_URL = "https://www.minorplanetcenter.net/iau/data/skycov.tgz"
COMET_ELEMENTS_URL = "https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt"
OPENNGC_URL = "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv"
LAMBERT_Y_SCALE = 180.0 / np.pi
STEREO_VIEW_RADIUS = 7.5
MER_DEFAULT_HORIZONTAL_FOV_DEG = 116.0

BAYER_GREEK = {
    "Alp": "α", "Bet": "β", "Gam": "γ", "Del": "δ", "Eps": "ε",
    "Zet": "ζ", "Eta": "η", "The": "θ", "Iot": "ι", "Kap": "κ",
    "Lam": "λ", "Mu": "μ", "Nu": "ν", "Xi": "ξ", "Omi": "ο",
    "Pi": "π", "Rho": "ρ", "Sig": "σ", "Tau": "τ", "Ups": "υ",
    "Phi": "φ", "Chi": "χ", "Psi": "ψ", "Ome": "ω",
}
SUPERSCRIPT_DIGITS = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

PLANET_STYLES = {
    "mercury": ("Me", (190, 195, 205, 255), 8),
    "venus": ("Ve", (255, 235, 175, 255), 11),
    "mars": ("Ma", (255, 115, 75, 255), 9),
    "jupiter": ("Ju", (245, 205, 145, 255), 13),
    "saturn": ("Sa", (235, 195, 105, 255), 12),
    "uranus": ("Ur", (130, 225, 235, 255), 9),
    "neptune": ("Ne", (100, 145, 255, 255), 9),
}


def greek_bayer_label(value):
    """Convert the abbreviated HYG Bayer field to a Greek label."""
    if not value:
        return ""
    base, separator, suffix = value.partition("-")
    greek = BAYER_GREEK.get(base, "")
    if not greek:
        return ""
    return greek + (suffix.translate(SUPERSCRIPT_DIGITS) if separator else "")


def star_color_from_index(color_index):
    """Return a restrained display colour from the HYG B-V colour index."""
    if not np.isfinite(color_index):
        return (225, 232, 255, 235)
    if color_index < -0.05:
        return (170, 205, 255, 240)
    if color_index < 0.35:
        return (205, 222, 255, 240)
    if color_index < 0.80:
        return (255, 250, 225, 240)
    if color_index < 1.30:
        return (255, 220, 170, 240)
    return (255, 180, 125, 240)


def format_sun_altitude(altitude):
    """Format the Sun's apparent altitude for the coordinate readout."""
    if altitude is None or not np.isfinite(altitude):
        return "Slunce: --"
    altitude = float(altitude)
    if abs(altitude) < 0.05:
        return "Slunce: 0.0° na obzoru"
    position = "nad obzorem" if altitude > 0.0 else "pod obzorem"
    return f"Slunce: {altitude:+.1f}° {position}"


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


def append_projected_segment(x_values, y_values, x1, y1, x2, y2, max_step=0.35):
    """Append a short segment in a projected X/Y plane."""
    if abs(x2 - x1) <= max_step and abs(y2 - y1) <= max_step:
        x_values.extend([x1, x2, np.nan])
        y_values.extend([y1, y2, np.nan])


class PointingDataReader(QThread):
    finished = pyqtSignal(object)
    progress = pyqtSignal(str)

    def __init__(self, dir_path, days_back, target_date):
        super().__init__()
        self.dir_path = dir_path
        self.days_back = days_back
        self.target_date = datetime(
            target_date.year, target_date.month, target_date.day
        )

    @staticmethod
    def files_in_date_window(files, target_date, days_back):
        """Return only MPC files valid for the displayed night's look-back window."""
        window_start = target_date - timedelta(days=days_back)
        return [record for record in files if window_start <= record[0] <= target_date]

    def parse_mpc_date(self, filename):
        base = os.path.basename(filename).split('.')[0].strip()
        if len(base) < 7:
            return None
        try:
            year = int(base[:4])
            if not (1990 <= year <= 2030):
                return None
            doy = int(base[4:])
            if 1 <= doy <= 366:
                return datetime(year, 1, 1) + timedelta(days=doy - 1)
        except Exception:
            return None
        return None

    def remove_empty_data_directories(self):
        """Remove empty subdirectories while always preserving the skycov root."""
        data_root = os.path.normcase(os.path.abspath(self.dir_path))
        removed_count = 0
        removal_errors = 0
        for current_root, _, _ in os.walk(self.dir_path, topdown=False):
            absolute_root = os.path.normcase(os.path.abspath(current_root))
            if absolute_root == data_root:
                continue
            try:
                if os.path.commonpath([data_root, absolute_root]) != data_root:
                    raise ValueError("Adresář leží mimo datovou složku")
                if not os.listdir(current_root):
                    os.rmdir(current_root)
                    removed_count += 1
            except FileNotFoundError:
                # Jiná spuštěná instance mohla stejnou prázdnou složku
                # odstranit mezi os.walk() a os.rmdir().
                continue
            except (OSError, ValueError):
                removal_errors += 1
        return removed_count, removal_errors

    def run(self):
        self.progress.emit("Indexuji místní archiv skycov...")
        if not self.dir_path or not os.path.exists(self.dir_path):
            self.progress.emit("Chyba: Složka neexistuje!")
            self.finished.emit(([], None))
            return

        data_root = os.path.normcase(os.path.abspath(self.dir_path))
        all_files = []
        invalid_files = []
        for root, _, files in os.walk(self.dir_path):
            stn_code = os.path.basename(root)
            for file in files:
                if file.lower().endswith('.dat'):
                    filepath = os.path.join(root, file)
                    dt = self.parse_mpc_date(file)
                    if dt:
                        all_files.append((dt, filepath, stn_code))
                    else:
                        invalid_files.append(filepath)

        removed_count = 0
        removed_bytes = 0
        removal_errors = 0
        for filepath in invalid_files:
            try:
                absolute_path = os.path.normcase(os.path.abspath(filepath))
                if os.path.commonpath([data_root, absolute_path]) != data_root:
                    raise ValueError("Soubor leží mimo datovou složku")
                file_size = os.path.getsize(filepath)
                os.remove(filepath)
                removed_count += 1
                removed_bytes += file_size
            except FileNotFoundError:
                continue
            except (OSError, ValueError):
                removal_errors += 1

        if not all_files:
            removed_dirs, directory_errors = self.remove_empty_data_directories()
            self.progress.emit(
                "Nenalezeny žádné platně datované .dat soubory! "
                f"Odstraněno neplatných souborů: {removed_count}. "
                f"Odstraněno prázdných složek: {removed_dirs}. "
                f"Chyby: {removal_errors + directory_errors}."
            )
            self.finished.emit(([], None))
            return

        max_date = max(pair[0] for pair in all_files)
        retention_date = max_date - timedelta(days=SKYCOV_RETENTION_DAYS)
        retained_files = []
        for record in all_files:
            dt, filepath, _ = record
            if dt >= retention_date:
                retained_files.append(record)
                continue
            try:
                absolute_path = os.path.normcase(os.path.abspath(filepath))
                if os.path.commonpath([data_root, absolute_path]) != data_root:
                    raise ValueError("Soubor leží mimo datovou složku")
                file_size = os.path.getsize(filepath)
                os.remove(filepath)
                removed_count += 1
                removed_bytes += file_size
            except FileNotFoundError:
                # Souběžné načtení v jiné instanci mohlo soubor odstranit dříve.
                continue
            except (OSError, ValueError):
                removal_errors += 1
                retained_files.append(record)

        all_files = retained_files
        removed_dirs, directory_errors = self.remove_empty_data_directories()
        min_date = self.target_date - timedelta(days=self.days_back)
        target_files = self.files_in_date_window(
            all_files, self.target_date, self.days_back
        )

        removed_mb = removed_bytes / (1024 * 1024)
        cleanup_status = f" Odstraněno: {removed_count} souborů ({removed_mb:.1f} MiB)."
        cleanup_status += f" Prázdné složky: {removed_dirs}."
        if removal_errors:
            cleanup_status += f" Nešlo odstranit souborů: {removal_errors}."
        if directory_errors:
            cleanup_status += f" Nešlo odstranit složek: {directory_errors}."
        self.progress.emit(
            f"Načítám {len(target_files)} souborů pro období "
            f"{min_date:%Y-%m-%d} až {self.target_date:%Y-%m-%d}.{cleanup_status}"
        )

        polygons = []
        for dt, filepath, stn in target_files:
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 9:
                            try:
                                coords = [float(x) for x in parts[:8]]
                                ra1, dec1 = np.degrees(coords[0]), np.degrees(coords[1])
                                ra2, dec2 = np.degrees(coords[2]), np.degrees(coords[3])
                                ra3, dec3 = np.degrees(coords[4]), np.degrees(coords[5])
                                ra4, dec4 = np.degrees(coords[6]), np.degrees(coords[7])

                                if any(abs(dec) > 90 for dec in (dec1, dec2, dec3, dec4)):
                                    continue

                                ra_list = [r % 360 for r in [ra1, ra2, ra3, ra4]]
                                dec_list = [dec1, dec2, dec3, dec4]

                                ra = ra_list + [ra_list[0]]
                                dec = dec_list + [dec_list[0]]
                                mag = float(parts[8])

                                polygons.append((ra, dec, mag))
                            except ValueError:
                                continue
            except Exception:
                continue

        ref_time = Time(max_date)
        if target_files:
            result_status = (
                f"Připraveno. Zobrazeno období {min_date:%Y-%m-%d} až "
                f"{self.target_date:%Y-%m-%d}."
            )
        else:
            result_status = (
                f"Pro období {min_date:%Y-%m-%d} až {self.target_date:%Y-%m-%d} "
                "nejsou k dispozici žádná data pokrytí."
            )
        self.progress.emit(
            f"{result_status} Nejnovější MPC data: {max_date:%Y-%m-%d}."
            f"{cleanup_status}"
        )
        self.finished.emit((polygons, ref_time))


class StarCatalogReader(QThread):
    finished = pyqtSignal(object)
    progress = pyqtSignal(str)

    def __init__(self, cache_path):
        super().__init__()
        self.cache_path = cache_path

    def run(self):
        url = "https://raw.githubusercontent.com/astronexus/HYG-Database/master/hyg/v3/hyg_v38.csv.gz"
        try:
            if not os.path.exists(self.cache_path):
                self.progress.emit("Stahuji katalog hvězd HYG...")
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=30) as response:
                    with open(self.cache_path, 'wb') as output:
                        output.write(response.read())

            stars = []
            with gzip.open(self.cache_path, 'rt', encoding='utf-8', errors='ignore', newline='') as source:
                for row in csv.DictReader(source):
                    try:
                        magnitude = float(row['mag'])
                        if magnitude <= 7.0:
                            try:
                                color_index = float(row.get('ci', 'nan'))
                            except (TypeError, ValueError):
                                color_index = float('nan')
                            stars.append((
                                float(row['ra']) * 15.0,
                                float(row['dec']),
                                magnitude,
                                greek_bayer_label(row.get('bayer', '')),
                                row.get('con', ''),
                                star_color_from_index(color_index),
                            ))
                    except (KeyError, TypeError, ValueError):
                        continue
            self.finished.emit((stars, ""))
        except Exception as error:
            self.finished.emit(([], f"Katalog hvězd: {error}"))


class SupplementCatalogReader(QThread):
    """Load cached OpenNGC objects and refresh MPC comet elements."""

    finished = pyqtSignal(object)
    progress = pyqtSignal(str)

    def __init__(self, cache_dir):
        super().__init__()
        self.cache_dir = cache_dir

    @staticmethod
    def download(url, path):
        temporary_path = path + ".tmp"
        request = urllib.request.Request(
            url, headers={"User-Agent": "SkyCoverage planner/1.0"}
        )
        with urllib.request.urlopen(request, timeout=60) as response, open(
                temporary_path, "wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        os.replace(temporary_path, path)

    @staticmethod
    def sexagesimal(value, hours=False):
        parts = value.strip().split(":")
        if len(parts) != 3:
            raise ValueError("Neplatná šedesátková souřadnice")
        sign = -1.0 if parts[0].startswith("-") else 1.0
        first = abs(float(parts[0]))
        result = first + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0
        return sign * result * (15.0 if hours else 1.0)

    @classmethod
    def load_deep_sky(cls, path):
        objects = []
        excluded_types = {"*", "**", "Nova", "NonEx", "Dup", "Other"}
        with open(path, "r", encoding="utf-8", errors="ignore", newline="") as source:
            for row in csv.DictReader(source, delimiter=";"):
                try:
                    object_type = row.get("Type", "")
                    if object_type in excluded_types:
                        continue
                    magnitude_text = row.get("V-Mag") or row.get("B-Mag")
                    if not magnitude_text:
                        continue
                    magnitude = float(magnitude_text)
                    if magnitude > 11.0:
                        continue
                    objects.append({
                        "ra": cls.sexagesimal(row["RA"], hours=True) % 360.0,
                        "dec": cls.sexagesimal(row["Dec"]),
                        "mag": magnitude,
                        "type": object_type,
                        "major": float(row.get("MajAx") or 0.0),
                        "minor": float(row.get("MinAx") or row.get("MajAx") or 0.0),
                        "angle": float(row.get("PosAng") or 0.0),
                    })
                except (KeyError, TypeError, ValueError):
                    continue
        return objects

    @staticmethod
    def load_comets(path):
        elements = []
        with open(path, "r", encoding="ascii", errors="ignore") as source:
            for line in source:
                if len(line) < 110:
                    continue
                try:
                    year = int(line[14:18])
                    month = int(line[19:21])
                    perihelion_day = float(line[22:29])
                    perihelion = datetime(year, month, 1) + timedelta(
                        days=perihelion_day - 1.0
                    )
                    elements.append({
                        "tp_jd": float(Time(perihelion, scale="tt").jd),
                        "q": float(line[30:39]),
                        "e": float(line[41:49]),
                        "arg_peri": float(line[51:59]),
                        "node": float(line[61:69]),
                        "inclination": float(line[71:79]),
                        "h": float(line[91:95]),
                        "slope": float(line[96:100]),
                        "name": line[102:158].strip() or line[5:12].strip(),
                    })
                except (TypeError, ValueError):
                    continue
        return elements

    def run(self):
        os.makedirs(self.cache_dir, exist_ok=True)
        dso_path = os.path.join(self.cache_dir, "openngc.csv")
        comet_path = os.path.join(self.cache_dir, "CometEls.txt")
        warnings = []
        try:
            if not os.path.exists(dso_path):
                self.progress.emit("Stahuji lehký katalog DSO OpenNGC…")
                self.download(OPENNGC_URL, dso_path)
        except Exception as error:
            warnings.append(f"OpenNGC: {error}")
        try:
            comet_age = (
                datetime.now().timestamp() - os.path.getmtime(comet_path)
                if os.path.exists(comet_path) else float("inf")
            )
            if comet_age > 24 * 3600:
                self.progress.emit("Aktualizuji kometární elementy z MPC…")
                self.download(COMET_ELEMENTS_URL, comet_path)
        except Exception as error:
            warnings.append(f"MPC komety: {error}")

        try:
            deep_sky = self.load_deep_sky(dso_path) if os.path.exists(dso_path) else []
        except Exception as error:
            deep_sky = []
            warnings.append(f"OpenNGC parser: {error}")
        try:
            comets = self.load_comets(comet_path) if os.path.exists(comet_path) else []
        except Exception as error:
            comets = []
            warnings.append(f"MPC parser komet: {error}")
        self.finished.emit((deep_sky, comets, "; ".join(warnings)))


class SkyCoverageUpdater(QThread):
    """Download and transactionally replace the local MPC Sky Coverage archive."""

    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, target_dir):
        super().__init__()
        self.target_dir = os.path.abspath(target_dir)

    @staticmethod
    def member_date(member_name):
        base = os.path.basename(member_name).rsplit('.', 1)[0]
        if len(base) != 7 or not base.isdigit():
            return None
        try:
            year, day = int(base[:4]), int(base[4:])
            if not 1990 <= year <= 2100 or not 1 <= day <= 366:
                return None
            return datetime(year, 1, 1) + timedelta(days=day - 1)
        except ValueError:
            return None

    @classmethod
    def newest_local_date(cls, directory):
        dates = []
        if os.path.isdir(directory):
            for _, _, files in os.walk(directory):
                for filename in files:
                    if filename.lower().endswith(".dat"):
                        date = cls.member_date(filename)
                        if date is not None:
                            dates.append(date)
        return max(dates) if dates else None

    def run(self):
        parent_dir = os.path.dirname(self.target_dir)
        work_dir = tempfile.mkdtemp(prefix="skycov_update_", dir=parent_dir)
        archive_path = os.path.join(work_dir, "skycov.tgz")
        staging_dir = os.path.join(work_dir, "new_data")
        backup_dir = os.path.join(work_dir, "old_data")
        try:
            previous_newest = self.newest_local_date(self.target_dir)
            self.progress.emit("Stahuji aktuální balíček Sky Coverage z MPC…")
            request = urllib.request.Request(
                SKYCOV_URL, headers={"User-Agent": "SkyCoverage planner/1.0"}
            )
            with urllib.request.urlopen(request, timeout=60) as response, open(archive_path, "wb") as output:
                total_size = int(response.headers.get("Content-Length") or 0)
                downloaded = 0
                last_percent = -1
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
                    downloaded += len(chunk)
                    if total_size:
                        percent = int(downloaded * 100 / total_size)
                        if percent >= last_percent + 2:
                            self.progress.emit(
                                f"Stahuji Sky Coverage z MPC… {percent} % "
                                f"({downloaded / 1024 ** 2:.0f} MiB)"
                            )
                            last_percent = percent

            self.progress.emit("Kontroluji archiv a vybírám poslední tři měsíce…")
            with tarfile.open(archive_path, "r:gz") as archive:
                dated_members = []
                for member in archive.getmembers():
                    if not member.isfile() or not member.name.lower().endswith(".dat"):
                        continue
                    date = self.member_date(member.name)
                    station = os.path.basename(os.path.dirname(member.name))
                    if date is not None and len(station) == 3 and station.isalnum():
                        dated_members.append((date, station, member))
                if not dated_members:
                    raise ValueError("Archiv neobsahuje platné datované soubory .DAT")
                newest = max(item[0] for item in dated_members)
                cutoff = newest - timedelta(days=SKYCOV_RETENTION_DAYS)
                selected = [item for item in dated_members if item[0] >= cutoff]
                os.makedirs(staging_dir)
                for _, station, member in selected:
                    station_dir = os.path.join(staging_dir, station)
                    os.makedirs(station_dir, exist_ok=True)
                    source = archive.extractfile(member)
                    if source is None:
                        continue
                    with source, open(os.path.join(station_dir, os.path.basename(member.name)), "wb") as output:
                        shutil.copyfileobj(source, output)

            self.progress.emit("Nahrazuji místní data ověřenou aktualizací…")
            if os.path.exists(self.target_dir):
                os.replace(self.target_dir, backup_dir)
            try:
                os.replace(staging_dir, self.target_dir)
            except Exception:
                if os.path.exists(backup_dir) and not os.path.exists(self.target_dir):
                    os.replace(backup_dir, self.target_dir)
                raise
            if os.path.exists(backup_dir):
                shutil.rmtree(backup_dir)
            if previous_newest is not None and newest <= previous_newest:
                detail = (
                    f"Balíček MPC neobsahuje novější data než "
                    f"{previous_newest:%Y-%m-%d}."
                )
            else:
                previous_text = (
                    f" (původně {previous_newest:%Y-%m-%d})"
                    if previous_newest is not None else ""
                )
                detail = f"Nejnovější data: {newest:%Y-%m-%d}{previous_text}."
            self.finished.emit(
                True,
                f"Aktualizace dokončena: {len(selected)} souborů. {detail}",
            )
        except Exception as error:
            self.finished.emit(False, f"Aktualizace MPC selhala: {error}")
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SkyCoverage — MPC sky map")
        self.resize(1600, 920)
        self.setMinimumSize(1180, 720)
        self.apply_modern_style()

        app_dir = os.path.dirname(os.path.abspath(__file__))
        self.app_dir = app_dir
        self.selected_dir = os.path.join(app_dir, "skycov")
        self.planner_store = PlannerStore(os.path.join(app_dir, "planner_data.json"))
        self.planned_fields = []
        self.preview_fields = []
        self.preview_items = []
        self.mercator_center_ra = 180.0
        self.mercator_center_dec = 0.0
        self.mercator_center_azimuth = 180.0
        self.mercator_center_altitude = 0.0
        self.local_mercator_center_azimuth = 180.0
        self.local_mercator_center_altitude = 45.0
        self.local_mercator_observatory = None
        self.last_cursor_icrs = None
        self.mer_polygon_cache = None
        self.reprojecting_mercator = False
        self.deep_sky_objects = []
        self.comet_elements = []
        self.comet_position_cache = None
        self.refreshing_planner_ui = False
        self.const_lines = self.load_constellation_lines()
        self.const_boundaries = self.load_constellation_boundaries()

        central_widget = QWidget()
        central_widget.setObjectName("centralWidget")
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(14, 12, 14, 10)
        main_layout.setSpacing(9)

        header_frame = QFrame()
        header_frame.setObjectName("headerCard")
        header_layout = QHBoxLayout(header_frame)
        header_layout.setContentsMargins(16, 9, 16, 9)
        brand_layout = QVBoxLayout()
        brand_layout.setSpacing(0)
        title_label = QLabel("SKYCOVERAGE")
        title_label.setObjectName("appTitle")
        subtitle_label = QLabel("MPC · mapa oblohy · plánování")
        subtitle_label.setObjectName("appSubtitle")
        brand_layout.addWidget(title_label)
        brand_layout.addWidget(subtitle_label)
        header_layout.addLayout(brand_layout)
        main_layout.addWidget(header_frame)

        # Panel 1: Pohled, Observatoř, Čas, Mag
        primary_controls = QFrame()
        primary_controls.setObjectName("controlCard")
        self.primary_controls_layout = QGridLayout(primary_controls)
        self.primary_controls_layout.setContentsMargins(12, 7, 12, 7)
        self.primary_controls_layout.setHorizontalSpacing(12)
        self.primary_controls_layout.setVerticalSpacing(5)
        control_panel1 = QHBoxLayout()
        control_panel1.setContentsMargins(0, 0, 0, 0)
        control_panel1.setSpacing(8)
        observing_controls = QHBoxLayout()
        observing_controls.setContentsMargins(0, 0, 0, 0)
        observing_controls.setSpacing(8)
        data_controls = QHBoxLayout()
        data_controls.setContentsMargins(0, 0, 0, 0)
        data_controls.setSpacing(8)
        extra_controls = QHBoxLayout()
        extra_controls.setContentsMargins(0, 0, 0, 0)
        extra_controls.setSpacing(8)
        # Na širokých obrazovkách využijeme prostor vedle ovládacích prvků,
        # místo abychom čtyřmi krátkými řádky zbytečně zvyšovali hlavičku.
        self.primary_control_rows = (
            control_panel1, observing_controls, data_controls, extra_controls,
        )
        self.primary_controls_wide = None
        self.reflow_primary_controls(self.width() >= 1550)

        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItems(["Globální", "Lokální obzor"])
        self.view_mode_combo.setMinimumWidth(145)
        self.view_mode_combo.setToolTip(
            "Globální rovníkové souřadnice RA/Dec nebo místní obzorové souřadnice Alt/Az."
        )
        self.view_mode_combo.setObjectName("accentControl")
        self.view_mode_combo.currentIndexChanged.connect(self.on_view_change)
        control_panel1.addWidget(self.view_mode_combo)

        control_panel1.addWidget(QLabel("Projekce:"))
        self.projection_combo = QComboBox()
        self.projection_combo.addItems([
            "Obdélníková",
            "Hammer-Aitoff",
            "Válcová plochojevná",
            "Stereografická",
            "MER Mercatorova",
        ])
        self.projection_combo.setCurrentIndex(4)
        self.projection_combo.setMinimumWidth(195)
        self.projection_combo.setToolTip(
            "MER je výchozí širokoúhlý pohled; stereografická projekce "
            "zachovává místní tvary a úhly souhvězdí"
        )
        self.projection_combo.currentIndexChanged.connect(self.on_projection_change)
        control_panel1.addWidget(self.projection_combo)

        self.zoom_out_button = QPushButton("−")
        self.zoom_out_button.setToolTip("Oddálit mapu")
        self.zoom_out_button.setFixedWidth(34)
        self.zoom_out_button.clicked.connect(lambda: self.zoom_map(1.35))
        extra_controls.addWidget(self.zoom_out_button)
        self.zoom_in_button = QPushButton("+")
        self.zoom_in_button.setToolTip("Přiblížit mapu")
        self.zoom_in_button.setFixedWidth(34)
        self.zoom_in_button.clicked.connect(lambda: self.zoom_map(0.74))
        extra_controls.addWidget(self.zoom_in_button)
        self.zoom_reset_button = QPushButton("Celá mapa")
        self.zoom_reset_button.setToolTip("Obnovit výchozí rozsah mapy")
        self.zoom_reset_button.clicked.connect(self.reset_map_view)
        extra_controls.addWidget(self.zoom_reset_button)

        observing_controls.addWidget(QLabel("Místo:"))
        self.obs_combo = QComboBox()
        self.obs_combo.addItems(list(OBSERVATORIES.keys()))
        self.obs_combo.setMinimumWidth(170)
        self.obs_combo.currentIndexChanged.connect(self.request_redraw)
        observing_controls.addWidget(self.obs_combo)

        observing_controls.addWidget(QLabel("UTC:"))
        self.time_edit = MinuteStepDateTimeEdit(QDateTime.currentDateTimeUtc())
        self.time_edit.setTimeSpec(Qt.TimeSpec.UTC)
        self.time_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.time_edit.setMinimumWidth(165)
        self.time_edit.setCurrentSection(QDateTimeEdit.Section.MinuteSection)
        self.time_edit.setToolTip(
            "Šipka nahoru přidá 1 minutu, šipka dolů 1 minutu ubere"
        )
        self.time_edit.dateTimeChanged.connect(self.request_redraw)
        observing_controls.addWidget(self.time_edit)
        observing_controls.addStretch()

        data_controls.addWidget(QLabel("Dny:"))
        self.days_spin = QSpinBox()
        self.days_spin.setRange(1, SKYCOV_RETENTION_DAYS)
        self.days_spin.setSingleStep(1)
        self.days_spin.setValue(30)
        self.days_spin.setMinimumWidth(78)
        self.days_spin.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.days_spin.setToolTip("Šipka nahoru přidá 1 den, šipka dolů 1 den ubere")
        data_controls.addWidget(self.days_spin)

        data_controls.addWidget(QLabel("Mag.:"))
        self.mag_spin = QDoubleSpinBox()
        self.mag_spin.setRange(5.0, 25.0)
        self.mag_spin.setDecimals(1)
        self.mag_spin.setSingleStep(0.5)
        self.mag_spin.setValue(15.0)
        self.mag_spin.setMinimumWidth(88)
        self.mag_spin.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.mag_spin.setToolTip("Šipka nahoru přidá 0,5 mag, šipka dolů 0,5 mag ubere")
        data_controls.addWidget(self.mag_spin)

        self.btn_load = QPushButton("Načíst")
        self.btn_load.setObjectName("primaryButton")
        self.btn_load.setToolTip("Načíst data pokrytí pro zvolený časový rozsah")
        self.btn_load.clicked.connect(self.start_processing)
        data_controls.addWidget(self.btn_load)

        self.kreutz_check = QCheckBox("Kreutzův koridor I/II")
        self.kreutz_check.setToolTip(
            "Model koridorů Kreutz I a II pro 5, 10, 15, 20, 30, 40 a 45 dní před přísluním; "
            "značky ukazují 5 dní a potom intervaly po 10 dnech"
        )
        self.kreutz_check.setObjectName("warmCheck")
        self.kreutz_check.stateChanged.connect(self.request_redraw)
        extra_controls.insertWidget(0, self.kreutz_check)
        control_panel1.addStretch()
        data_controls.addStretch()
        extra_controls.addStretch()

        self.coord_label = QLabel(
            "Az: -- | Alt: --\n"
            "RA: -- | Dec: --\n"
            "Elongace: --\n"
            "Slunce: --"
        )
        self.coord_label.setTextFormat(Qt.TextFormat.PlainText)
        self.coord_label.setWordWrap(False)
        self.coord_label.setFixedSize(270, 90)
        self.coord_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self.coord_label.setObjectName("coordReadout")
        header_layout.addWidget(primary_controls, 1)
        header_layout.addWidget(self.coord_label)

        # Panel 2: Tlačítka soumraků
        secondary_controls = QFrame()
        secondary_controls.setObjectName("controlCard")
        secondary_layout = QVBoxLayout(secondary_controls)
        secondary_layout.setContentsMargins(12, 8, 12, 8)
        secondary_layout.setSpacing(5)
        control_panel2 = QHBoxLayout()
        control_panel2.setContentsMargins(0, 0, 0, 0)
        control_panel2.setSpacing(7)
        layer_controls = QHBoxLayout()
        layer_controls.setContentsMargins(0, 0, 0, 0)
        layer_controls.setSpacing(12)
        secondary_layout.addLayout(control_panel2)
        secondary_layout.addLayout(layer_controls)
        control_panel2.addWidget(QLabel("Soumrak:"))

        self.btn_naut_eve = QPushButton("Večer −12°")
        self.btn_naut_eve.setToolTip("Večerní nautický soumrak: Slunce 12° pod obzorem")
        self.btn_naut_eve.clicked.connect(lambda: self.set_twilight(-12.0, evening=True))
        control_panel2.addWidget(self.btn_naut_eve)

        self.btn_astro_eve = QPushButton("Večer −18°")
        self.btn_astro_eve.setToolTip("Začátek večerní astronomické noci: Slunce 18° pod obzorem")
        self.btn_astro_eve.clicked.connect(lambda: self.set_twilight(-18.0, evening=True))
        control_panel2.addWidget(self.btn_astro_eve)

        self.btn_astro_morn = QPushButton("Ráno −18°")
        self.btn_astro_morn.setToolTip("Konec ranní astronomické noci: Slunce 18° pod obzorem")
        self.btn_astro_morn.clicked.connect(lambda: self.set_twilight(-18.0, evening=False))
        control_panel2.addWidget(self.btn_astro_morn)

        self.btn_naut_morn = QPushButton("Ráno −12°")
        self.btn_naut_morn.setToolTip("Ranní nautický soumrak: Slunce 12° pod obzorem")
        self.btn_naut_morn.clicked.connect(lambda: self.set_twilight(-12.0, evening=False))
        control_panel2.addWidget(self.btn_naut_morn)

        for label, minutes in [("+30 min", 30), ("+1 h", 60), ("-1 h", -60), ("-30 min", -30)]:
            button = QPushButton(label)
            button.clicked.connect(lambda checked=False, delta=minutes: self.shift_time(delta))
            control_panel2.addWidget(button)

        self.const_boundaries_check = QCheckBox("Hranice souhvězdí")
        self.const_boundaries_check.setChecked(True)
        self.const_boundaries_check.setToolTip(
            "Zobrazit nebo skrýt oficiální hranice 88 souhvězdí"
        )
        self.const_boundaries_check.setObjectName("coolCheck")
        self.const_boundaries_check.stateChanged.connect(self.request_redraw)
        layer_controls.addWidget(self.const_boundaries_check)

        self.bayer_labels_check = QCheckBox("Bayer αβγ")
        self.bayer_labels_check.setChecked(True)
        self.bayer_labels_check.setToolTip(
            "Zobrazit řecká Bayerova označení jasných hvězd; hustota se přizpůsobuje přiblížení"
        )
        self.bayer_labels_check.setObjectName("coolCheck")
        self.bayer_labels_check.stateChanged.connect(self.request_redraw)
        layer_controls.addWidget(self.bayer_labels_check)

        self.planets_check = QCheckBox("Planety")
        self.planets_check.setChecked(True)
        self.planets_check.setToolTip(
            "Zobrazit nebo skrýt Merkur, Venuši, Mars, Jupiter, Saturn, Uran a Neptun"
        )
        self.planets_check.setObjectName("coolCheck")
        self.planets_check.stateChanged.connect(self.request_redraw)
        layer_controls.addWidget(self.planets_check)

        self.dso_check = QCheckBox("DSO ≤11m")
        self.dso_check.setChecked(True)
        self.dso_check.setToolTip(
            "Při větším přiblížení zobrazit objekty OpenNGC bez popisků"
        )
        self.dso_check.setObjectName("coolCheck")
        self.dso_check.stateChanged.connect(self.request_redraw)
        layer_controls.addWidget(self.dso_check)

        self.comets_check = QCheckBox("Komety ≤17m")
        self.comets_check.setChecked(True)
        self.comets_check.setToolTip(
            "Zobrazit komety nad obzorem podle aktuálních elementů MPC"
        )
        self.comets_check.setObjectName("warmCheck")
        self.comets_check.stateChanged.connect(self.request_redraw)
        layer_controls.addWidget(self.comets_check)

        control_panel2.addStretch()
        layer_controls.addStretch()
        main_layout.addWidget(secondary_controls)

        self.workspace_tabs = QTabWidget()
        self.workspace_tabs.setObjectName("workspaceTabs")
        self.map_workspace = QWidget()
        self.map_workspace.setObjectName("workspacePage")
        self.map_workspace_layout = QVBoxLayout(self.map_workspace)
        self.map_workspace_layout.setContentsMargins(0, 4, 0, 0)
        self.planning_workspace = QWidget()
        self.planning_workspace.setObjectName("workspacePage")
        self.planning_workspace_layout = QVBoxLayout(self.planning_workspace)
        self.planning_workspace_layout.setContentsMargins(0, 4, 0, 0)
        self.workspace_tabs.addTab(self.map_workspace, "Mapa pokrytí")
        self.workspace_tabs.addTab(self.planning_workspace, "Plánování pozorování")
        self.workspace_tabs.currentChanged.connect(self.on_workspace_tab_changed)
        main_layout.addWidget(self.workspace_tabs, 1)

        self.planner_tabs = QTabWidget()
        self.planner_tabs.setObjectName("plannerTabs")
        self.planner_tabs.setMinimumWidth(430)
        self.planner_tabs.setMaximumWidth(560)
        self.planning_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.planning_map_panel = QWidget()
        self.planning_map_layout = QVBoxLayout(self.planning_map_panel)
        self.planning_map_layout.setContentsMargins(0, 0, 0, 0)
        self.planning_splitter.addWidget(self.planning_map_panel)
        self.planning_splitter.addWidget(self.planner_tabs)
        self.planning_splitter.setStretchFactor(0, 5)
        self.planning_splitter.setStretchFactor(1, 0)
        self.planning_splitter.setSizes([1100, 460])
        self.planning_workspace_layout.addWidget(self.planning_splitter)

        planner_page = QWidget()
        planner_layout = QVBoxLayout(planner_page)
        planner_layout.setContentsMargins(10, 8, 10, 8)
        planner_controls = QGridLayout()
        planner_controls.setSpacing(6)

        self.placement_check = QCheckBox("Kliknutím v mapě přidávat pole")
        self.placement_check.setObjectName("warmCheck")
        self.placement_check.stateChanged.connect(self.on_placement_mode_changed)
        planner_controls.addWidget(self.placement_check, 0, 0, 1, 4)
        planner_controls.addWidget(QLabel("Sestava:"), 1, 0)
        self.plan_profile_combo = QComboBox()
        planner_controls.addWidget(self.plan_profile_combo, 1, 1, 1, 3)
        planner_controls.addWidget(QLabel("Mozaika:"), 2, 0, 1, 4)
        mosaic_controls = QGridLayout()
        mosaic_controls.setContentsMargins(0, 0, 0, 0)
        mosaic_controls.setHorizontalSpacing(8)
        mosaic_controls.addWidget(QLabel("Řádky"), 0, 0)
        mosaic_controls.addWidget(QLabel("Sloupce"), 0, 1)
        self.mosaic_overlap_label = QLabel("Překryv")
        mosaic_controls.addWidget(self.mosaic_overlap_label, 0, 2)
        self.mosaic_rows_spin = QSpinBox()
        self.mosaic_rows_spin.setRange(1, 20)
        self.mosaic_rows_spin.setValue(1)
        self.mosaic_rows_spin.setMinimumWidth(78)
        mosaic_controls.addWidget(self.mosaic_rows_spin, 1, 0)
        self.mosaic_columns_spin = QSpinBox()
        self.mosaic_columns_spin.setRange(1, 20)
        self.mosaic_columns_spin.setValue(1)
        self.mosaic_columns_spin.setMinimumWidth(78)
        mosaic_controls.addWidget(self.mosaic_columns_spin, 1, 1)
        self.mosaic_overlap_spin = QDoubleSpinBox()
        self.mosaic_overlap_spin.setRange(0, 90)
        self.mosaic_overlap_spin.setDecimals(1)
        self.mosaic_overlap_spin.setValue(10)
        self.mosaic_overlap_spin.setMinimumWidth(145)
        self.mosaic_overlap_spin.setSuffix(" %")
        mosaic_controls.addWidget(self.mosaic_overlap_spin, 1, 2)
        planner_controls.addLayout(mosaic_controls, 3, 0, 1, 4)
        orientation_label = QLabel("Orientace:")
        orientation_label.setMinimumWidth(120)
        planner_controls.addWidget(orientation_label, 4, 0)
        self.mosaic_orientation_combo = QComboBox()
        self.mosaic_orientation_combo.addItems([
            "Podle EQ sestavy",
            "Vodorovně s obzorem (EQ)",
        ])
        self.mosaic_orientation_combo.setToolTip(
            "Vodorovný režim rozloží středy polí podél místního obzoru pro čas mapy. "
            "Kamera zůstává orientovaná podle EQ montáže; rozteč se automaticky zmenší, "
            "aby zůstal zachován zadaný minimální překryv i u nakloněných polí."
        )
        self.mosaic_orientation_combo.currentIndexChanged.connect(
            self.on_mosaic_orientation_changed
        )
        planner_controls.addWidget(self.mosaic_orientation_combo, 4, 1, 1, 3)
        capture_controls = QGridLayout()
        capture_controls.setContentsMargins(0, 0, 0, 0)
        capture_controls.setHorizontalSpacing(8)
        capture_controls.addWidget(QLabel("Filtr"), 0, 0)
        capture_controls.addWidget(QLabel("Expozice"), 0, 1)
        capture_controls.addWidget(QLabel("Cykly"), 0, 2)
        self.filter_edit = QLineEdit("N")
        self.filter_edit.setMinimumWidth(70)
        self.filter_edit.setMaximumWidth(90)
        capture_controls.addWidget(self.filter_edit, 1, 0)
        self.exposure_spin = QDoubleSpinBox()
        self.exposure_spin.setRange(0.01, 86400)
        self.exposure_spin.setDecimals(2)
        self.exposure_spin.setValue(30)
        self.exposure_spin.setSuffix(" s")
        self.exposure_spin.setMinimumWidth(155)
        capture_controls.addWidget(self.exposure_spin, 1, 1)
        self.loops_spin = QSpinBox()
        self.loops_spin.setRange(1, 100000)
        self.loops_spin.setValue(15)
        self.loops_spin.setMinimumWidth(145)
        capture_controls.addWidget(self.loops_spin, 1, 2)
        planner_controls.addLayout(capture_controls, 5, 0, 1, 4)
        button_row = QHBoxLayout()
        self.undo_field_button = QPushButton("Zpět")
        self.undo_field_button.clicked.connect(self.undo_planned_field)
        button_row.addWidget(self.undo_field_button)
        self.clear_fields_button = QPushButton("Vymazat plán")
        self.clear_fields_button.clicked.connect(self.clear_planned_fields)
        button_row.addWidget(self.clear_fields_button)
        self.archive_plan_button = QPushButton("Uložit plán…")
        self.archive_plan_button.setObjectName("primaryButton")
        self.archive_plan_button.setToolTip("Uložit plán do archivu a exportovat textový soubor")
        self.archive_plan_button.clicked.connect(self.archive_current_plan)
        button_row.addWidget(self.archive_plan_button)
        planner_controls.addLayout(button_row, 6, 0, 1, 4)
        planner_layout.addLayout(planner_controls, 2)

        self.fields_table = QTableWidget(0, 4)
        self.fields_table.setHorizontalHeaderLabels(["#", "RA", "Dec", "Sestava"])
        self.fields_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.fields_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.fields_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        planner_layout.addWidget(self.fields_table, 2)

        self.sequence_duration_label = QLabel("Odhad délky sekvence: --")
        self.sequence_duration_label.setWordWrap(True)
        self.sequence_duration_label.setToolTip(
            "Expozice a vyčítání se počítají pro každý snímek. Přejezd se "
            "počítá mezi poli a při dalším cyklu také z posledního pole zpět "
            "na první. Úvodní najetí z neznámé výchozí polohy zahrnuto není."
        )
        planner_layout.addWidget(self.sequence_duration_label)

        self.script_output = QPlainTextEdit()
        self.script_output.setReadOnly(True)
        self.script_output.setPlaceholderText("Expoziční skript se vytvoří z polí v mapě.")
        planner_layout.addWidget(self.script_output, 3)
        self.planner_tabs.addTab(planner_page, "Plán")

        equipment_page = QWidget()
        equipment_layout = QGridLayout(equipment_page)
        equipment_layout.setContentsMargins(10, 8, 10, 8)
        equipment_layout.addWidget(QLabel("Profil sestavy:"), 0, 0, 1, 2)
        self.equipment_combo = QComboBox()
        self.equipment_combo.currentIndexChanged.connect(self.load_profile_into_editor)
        equipment_layout.addWidget(self.equipment_combo, 0, 2, 1, 2)
        equipment_layout.addWidget(QLabel("Název:"), 1, 0)
        self.profile_name_edit = QLineEdit()
        equipment_layout.addWidget(self.profile_name_edit, 1, 1, 1, 3)
        equipment_layout.addWidget(QLabel("Šířka pole"), 2, 0, 1, 2)
        equipment_layout.addWidget(QLabel("Výška pole"), 2, 2, 1, 2)
        self.profile_width_spin = QDoubleSpinBox()
        self.profile_width_spin.setRange(0.01, 180)
        self.profile_width_spin.setDecimals(3)
        self.profile_width_spin.setSuffix("°")
        self.profile_width_spin.setMinimumWidth(90)
        equipment_layout.addWidget(self.profile_width_spin, 3, 0, 1, 2)
        self.profile_height_spin = QDoubleSpinBox()
        self.profile_height_spin.setRange(0.01, 180)
        self.profile_height_spin.setDecimals(3)
        self.profile_height_spin.setSuffix("°")
        self.profile_height_spin.setMinimumWidth(90)
        equipment_layout.addWidget(self.profile_height_spin, 3, 2, 1, 2)
        equipment_layout.addWidget(QLabel("Rotace"), 4, 0, 1, 2)
        equipment_layout.addWidget(QLabel("Dithering"), 4, 2, 1, 2)
        self.profile_rotation_spin = QDoubleSpinBox()
        self.profile_rotation_spin.setRange(-180, 180)
        self.profile_rotation_spin.setDecimals(2)
        self.profile_rotation_spin.setSuffix("°")
        self.profile_rotation_spin.setMinimumWidth(90)
        equipment_layout.addWidget(self.profile_rotation_spin, 5, 0, 1, 2)
        self.profile_dither_spin = QDoubleSpinBox()
        self.profile_dither_spin.setRange(0, 20)
        self.profile_dither_spin.setDecimals(4)
        self.profile_dither_spin.setSuffix("°")
        self.profile_dither_spin.setMinimumWidth(90)
        equipment_layout.addWidget(self.profile_dither_spin, 5, 2, 1, 2)
        readout_label = QLabel("Vyčtení jednoho snímku")
        readout_label.setWordWrap(True)
        equipment_layout.addWidget(readout_label, 6, 0, 1, 2)
        slew_label = QLabel("Přejezd mezi poli")
        slew_label.setWordWrap(True)
        equipment_layout.addWidget(slew_label, 6, 2, 1, 2)
        self.profile_readout_spin = QDoubleSpinBox()
        self.profile_readout_spin.setRange(0, 3600)
        self.profile_readout_spin.setDecimals(1)
        self.profile_readout_spin.setSuffix(" s")
        self.profile_readout_spin.setMinimumWidth(130)
        self.profile_readout_spin.setToolTip(
            "Doba vyčtení a uložení jednoho CCD/CMOS snímku."
        )
        equipment_layout.addWidget(self.profile_readout_spin, 7, 0, 1, 2)
        self.profile_slew_spin = QDoubleSpinBox()
        self.profile_slew_spin.setRange(0, 3600)
        self.profile_slew_spin.setDecimals(1)
        self.profile_slew_spin.setSuffix(" s")
        self.profile_slew_spin.setMinimumWidth(130)
        self.profile_slew_spin.setToolTip(
            "Průměrný mrtvý čas jednoho přejezdu montáže mezi poli."
        )
        equipment_layout.addWidget(self.profile_slew_spin, 7, 2, 1, 2)
        self.save_profile_button = QPushButton("Uložit")
        self.save_profile_button.setObjectName("primaryButton")
        self.save_profile_button.setToolTip("Uložit změny profilu sestavy")
        self.save_profile_button.clicked.connect(self.save_equipment_profile)
        equipment_layout.addWidget(self.save_profile_button, 8, 0)
        self.new_profile_button = QPushButton("Nový")
        self.new_profile_button.setToolTip("Vytvořit nový profil sestavy")
        self.new_profile_button.clicked.connect(self.new_equipment_profile)
        equipment_layout.addWidget(self.new_profile_button, 8, 1)
        self.delete_profile_button = QPushButton("Odstranit")
        self.delete_profile_button.setToolTip("Odstranit vybraný profil sestavy")
        self.delete_profile_button.clicked.connect(self.delete_equipment_profile)
        equipment_layout.addWidget(self.delete_profile_button, 8, 2, 1, 2)
        equipment_layout.setColumnStretch(3, 1)
        self.planner_tabs.addTab(equipment_page, "Vybavení")

        archive_page = QWidget()
        archive_layout = QVBoxLayout(archive_page)
        archive_layout.setContentsMargins(10, 8, 10, 8)
        self.archive_table = QTableWidget(0, 5)
        self.archive_table.setHorizontalHeaderLabels(["Zobrazit", "Plán", "Uloženo UTC", "Polí", "Sestava"])
        self.archive_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.archive_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.archive_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.archive_table.itemChanged.connect(self.archive_visibility_changed)
        archive_layout.addWidget(self.archive_table)
        archive_buttons = QHBoxLayout()
        load_archive_button = QPushButton("Načíst plán")
        load_archive_button.setToolTip("Načíst vybraný archiv jako aktivní plán")
        load_archive_button.clicked.connect(self.load_selected_archive)
        archive_buttons.addWidget(load_archive_button)
        delete_archive_button = QPushButton("Odstranit")
        delete_archive_button.setToolTip("Odstranit vybraný plán z archivu")
        delete_archive_button.clicked.connect(self.delete_selected_archive)
        archive_buttons.addWidget(delete_archive_button)
        archive_buttons.addStretch()
        archive_layout.addLayout(archive_buttons)
        self.planner_tabs.addTab(archive_page, "Archiv")

        data_page = QWidget()
        data_layout = QHBoxLayout(data_page)
        data_layout.setContentsMargins(12, 12, 12, 12)
        self.data_age_label = QLabel("Nejnovější MPC data: zjišťuji…")
        self.data_age_label.setObjectName("dataAgeLabel")
        self.data_age_label.setWordWrap(True)
        data_layout.addWidget(self.data_age_label)
        data_layout.addStretch()
        self.update_data_button = QPushButton("Aktualizovat")
        self.update_data_button.setObjectName("primaryButton")
        self.update_data_button.setToolTip("Stáhnout a ověřit aktuální data Sky Coverage z MPC")
        self.update_data_button.clicked.connect(self.update_skycoverage_data)
        data_layout.addWidget(self.update_data_button)
        self.planner_tabs.addTab(data_page, "Data MPC")

        # Plot
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setObjectName("skyPlot")
        self.plot_widget.setBackground('#060a18')
        self.plot_widget.getViewBox().invertX(True)
        self.plot_widget.getViewBox().setAspectLocked(True, ratio=1.0)
        axis_font = QFont("Segoe UI", 9)
        for axis_name in ("left", "bottom"):
            axis = self.plot_widget.getAxis(axis_name)
            axis.setPen(pg.mkPen((70, 89, 118, 180), width=1))
            axis.setTextPen(pg.mkPen((142, 158, 184, 220)))
            axis.setTickFont(axis_font)
        self.plot_widget.scene().sigMouseMoved.connect(self.on_mouse_move)
        self.plot_widget.scene().sigMouseClicked.connect(self.on_map_click)
        self.map_workspace_layout.addWidget(self.plot_widget)

        self.status_label = QLabel("Připraveno.")
        self.status_label.setObjectName("statusLabel")
        main_layout.addWidget(self.status_label)

        self.polygons = []
        self.stars = []
        self.ref_time = None
        self.local_sun_azimuth = None
        self.current_sun_altitude = None
        self.uncovered_cache = {}
        self.bayer_label_data = []
        self.bayer_label_items = []
        self.redraw_timer = QTimer(self)
        self.redraw_timer.setSingleShot(True)
        self.redraw_timer.timeout.connect(self.redraw)
        self.reload_timer = QTimer(self)
        self.reload_timer.setSingleShot(True)
        self.reload_timer.timeout.connect(self.start_processing)
        self.axis_tick_timer = QTimer(self)
        self.axis_tick_timer.setSingleShot(True)
        self.axis_tick_timer.timeout.connect(self.update_axis_ticks)
        self.bayer_label_timer = QTimer(self)
        self.bayer_label_timer.setSingleShot(True)
        self.bayer_label_timer.timeout.connect(self.refresh_bayer_labels)
        self.mer_pan_timer = QTimer(self)
        self.mer_pan_timer.setSingleShot(True)
        self.mer_pan_timer.timeout.connect(self.refresh_mercator_view)
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_axis_tick_update
        )
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_bayer_label_update
        )
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_mercator_refresh
        )

        self.days_spin.valueChanged.connect(self.request_data_reload)
        self.mag_spin.valueChanged.connect(self.request_redraw)
        self.obs_combo.currentIndexChanged.connect(self.request_redraw)
        self.time_edit.dateTimeChanged.connect(self.request_redraw)
        self.time_edit.dateTimeChanged.connect(self.request_data_reload)
        self.filter_edit.textChanged.connect(self.refresh_exposure_script)
        self.exposure_spin.valueChanged.connect(self.refresh_exposure_script)
        self.loops_spin.valueChanged.connect(self.refresh_exposure_script)
        self.refresh_equipment_ui()
        self.refresh_archive_table()
        self.refresh_plan_table()
        star_cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "astro_cache", "hyg_v38.csv.gz")
        self.star_reader = StarCatalogReader(star_cache)
        self.star_reader.finished.connect(self.on_stars_loaded)
        self.star_reader.progress.connect(self.status_label.setText)
        self.star_reader.start()
        self.supplement_reader = SupplementCatalogReader(os.path.dirname(star_cache))
        self.supplement_reader.finished.connect(self.on_supplement_catalogs_loaded)
        self.supplement_reader.progress.connect(self.status_label.setText)
        self.supplement_reader.start()
        self.set_default_mer_range()
        self.update_axis_ticks()
        self.start_processing()

    def reflow_primary_controls(self, wide):
        """Use two compact rows when there is room, four rows in narrow windows."""
        wide = bool(wide)
        if self.primary_controls_wide == wide:
            return

        for row_layout in self.primary_control_rows:
            self.primary_controls_layout.removeItem(row_layout)

        if wide:
            for index, row_layout in enumerate(self.primary_control_rows):
                self.primary_controls_layout.addLayout(
                    row_layout, index // 2, index % 2
                )
            self.primary_controls_layout.setColumnStretch(1, 1)
        else:
            for index, row_layout in enumerate(self.primary_control_rows):
                self.primary_controls_layout.addLayout(row_layout, index, 0, 1, 2)
            self.primary_controls_layout.setColumnStretch(1, 0)

        self.primary_controls_wide = wide
        self.primary_controls_layout.invalidate()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "primary_control_rows"):
            self.reflow_primary_controls(event.size().width() >= 1550)

    def apply_modern_style(self):
        self.setStyleSheet("""
            QMainWindow, QWidget#centralWidget {
                background-color: #070b16;
                color: #dbe7f7;
                font-family: "Segoe UI";
                font-size: 12px;
            }
            QFrame#headerCard {
                background-color: #0d1525;
                border: 1px solid #1d2a42;
                border-radius: 12px;
            }
            QFrame#controlCard {
                background-color: #0d1422;
                border: 1px solid #1b2940;
                border-radius: 10px;
            }
            QLabel {
                color: #9cabc1;
                background: transparent;
                border: none;
            }
            QLabel#appTitle {
                color: #f2f7ff;
                font-size: 18px;
                font-weight: 700;
                letter-spacing: 2px;
            }
            QLabel#appSubtitle {
                color: #6f829e;
                font-size: 11px;
            }
            QLabel#coordReadout {
                color: #72e6d1;
                background-color: #09111e;
                border: 1px solid #1d3c48;
                border-radius: 8px;
                padding: 5px 10px;
                font-family: "Cascadia Mono", "Consolas";
                font-size: 11px;
                font-weight: 600;
            }
            QLabel#statusLabel {
                color: #71829a;
                padding: 1px 6px;
                font-size: 11px;
            }
            QComboBox, QSpinBox, QDoubleSpinBox, QDateTimeEdit, QLineEdit, QPlainTextEdit,
            QTableWidget {
                color: #e8f0fb;
                background-color: #152137;
                border: 1px solid #2a3a56;
                border-radius: 7px;
                padding: 6px 9px;
                min-height: 20px;
                selection-background-color: #245e72;
            }
            QSpinBox, QDoubleSpinBox, QDateTimeEdit {
                color: #f4f8ff;
                background-color: #111c30;
                padding-left: 8px;
                padding-right: 31px;
                font-weight: 600;
            }
            QSpinBox::up-button, QDoubleSpinBox::up-button, QDateTimeEdit::up-button {
                subcontrol-origin: border;
                subcontrol-position: top right;
                width: 27px;
                background-color: #20314d;
                border-left: 1px solid #3a4e6d;
                border-bottom: 1px solid #3a4e6d;
                border-top-right-radius: 6px;
            }
            QSpinBox::down-button, QDoubleSpinBox::down-button, QDateTimeEdit::down-button {
                subcontrol-origin: border;
                subcontrol-position: bottom right;
                width: 27px;
                background-color: #20314d;
                border-left: 1px solid #3a4e6d;
                border-bottom-right-radius: 6px;
            }
            QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
            QDateTimeEdit::up-button:hover, QSpinBox::down-button:hover,
            QDoubleSpinBox::down-button:hover, QDateTimeEdit::down-button:hover {
                background-color: #2b4567;
            }
            QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover, QDateTimeEdit:hover,
            QLineEdit:hover, QPlainTextEdit:hover, QTableWidget:hover {
                border-color: #3e6380;
                background-color: #182740;
            }
            QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus, QDateTimeEdit:focus,
            QLineEdit:focus, QPlainTextEdit:focus, QTableWidget:focus {
                border-color: #42cbb8;
            }
            QComboBox#accentControl {
                color: #75ead5;
                font-weight: 600;
            }
            QComboBox QAbstractItemView {
                color: #e8f0fb;
                background-color: #111c30;
                border: 1px solid #2a3a56;
                selection-background-color: #22445d;
                outline: none;
            }
            QPushButton {
                color: #cbd8e9;
                background-color: #17243a;
                border: 1px solid #2a3b57;
                border-radius: 7px;
                padding: 6px 10px;
                min-height: 20px;
            }
            QPushButton:hover {
                color: #ffffff;
                background-color: #20314d;
                border-color: #46617f;
            }
            QPushButton:pressed {
                background-color: #101a2c;
            }
            QPushButton#primaryButton {
                color: #061512;
                background-color: #52d5b9;
                border-color: #69e4ca;
                font-weight: 700;
                padding-left: 14px;
                padding-right: 14px;
            }
            QPushButton#primaryButton:hover {
                background-color: #6ae5ca;
            }
            QPushButton:disabled {
                color: #59677b;
                background-color: #111a2a;
                border-color: #1e2a3d;
            }
            QCheckBox {
                color: #b9c7d9;
                spacing: 6px;
                background: transparent;
            }
            QCheckBox#coolCheck { color: #8fd2bd; font-weight: 600; }
            QCheckBox#warmCheck { color: #ffad82; font-weight: 600; }
            QCheckBox::indicator {
                width: 15px;
                height: 15px;
                border: 1px solid #42536e;
                border-radius: 4px;
                background-color: #111b2d;
            }
            QCheckBox::indicator:checked {
                background-color: #49cbb3;
                border-color: #70e2cf;
            }
            QTabWidget#plannerTabs::pane {
                background-color: #0d1422;
                border: 1px solid #1b2940;
                border-radius: 8px;
            }
            QTabWidget#workspaceTabs::pane, QWidget#workspacePage {
                background-color: #070b16;
                border: none;
            }
            QTabBar::tab {
                color: #8fa2bb;
                background-color: #101a2b;
                border: 1px solid #1b2940;
                padding: 6px 16px;
            }
            QTabBar::tab:selected { color: #72e6d1; background-color: #17243a; }
            QHeaderView::section {
                color: #9fb0c6;
                background-color: #101b2d;
                border: none;
                padding: 4px;
            }
            QLabel#dataAgeLabel { color: #72e6d1; font-size: 14px; font-weight: 600; }
            QToolTip {
                color: #edf5ff;
                background-color: #111c2e;
                border: 1px solid #35506c;
                padding: 5px;
            }
            QGraphicsView#skyPlot {
                border: 1px solid #1e2c45;
                border-radius: 10px;
            }
        """)

    def load_cached_geojson(self, filename, url):
        cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "astro_cache")
        os.makedirs(cache_dir, exist_ok=True)
        json_path = os.path.join(cache_dir, filename)

        if not os.path.exists(json_path):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=10) as resp:
                    with open(json_path, 'wb') as f:
                        f.write(resp.read())
            except Exception:
                return []

        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception:
            return {}

    def load_constellation_lines(self):
        data = self.load_cached_geojson(
            "constellations.lines.json",
            "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json",
        )
        lines_coords = []
        for feature in data.get("features", []):
            geom = feature.get("geometry", {})
            gtype = geom.get("type", "")
            coords = geom.get("coordinates", [])
            if gtype == "MultiLineString":
                lines_coords.extend(coords)
            elif gtype == "LineString":
                lines_coords.append(coords)
        return lines_coords

    def load_constellation_boundaries(self):
        data = self.load_cached_geojson(
            "constellations.bounds.json",
            "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.bounds.json",
        )
        boundaries = []
        for feature in data.get("features", []):
            geom = feature.get("geometry", {})
            gtype = geom.get("type", "")
            coords = geom.get("coordinates", [])
            if gtype == "Polygon":
                boundaries.extend(coords)
            elif gtype == "MultiPolygon":
                for polygon in coords:
                    boundaries.extend(polygon)
        return boundaries

    def refresh_equipment_ui(self):
        profiles = self.planner_store.data["profiles"]
        selected_name = self.plan_profile_combo.currentText() if self.plan_profile_combo.count() else ""
        self.plan_profile_combo.blockSignals(True)
        self.equipment_combo.blockSignals(True)
        self.plan_profile_combo.clear()
        self.equipment_combo.clear()
        for profile in profiles:
            self.plan_profile_combo.addItem(profile["name"])
            self.equipment_combo.addItem(profile["name"])
        index = max(0, self.plan_profile_combo.findText(selected_name))
        self.plan_profile_combo.setCurrentIndex(index)
        self.equipment_combo.setCurrentIndex(index)
        self.plan_profile_combo.blockSignals(False)
        self.equipment_combo.blockSignals(False)
        self.load_profile_into_editor(index)

    def current_profile(self):
        index = self.plan_profile_combo.currentIndex()
        profiles = self.planner_store.data["profiles"]
        return profiles[index if 0 <= index < len(profiles) else 0]

    def load_profile_into_editor(self, index):
        profiles = self.planner_store.data["profiles"]
        if not 0 <= index < len(profiles):
            return
        profile = profiles[index]
        self.profile_name_edit.setText(str(profile["name"]))
        self.profile_width_spin.setValue(float(profile["width_deg"]))
        self.profile_height_spin.setValue(float(profile["height_deg"]))
        self.profile_rotation_spin.setValue(float(profile["rotation_deg"]))
        self.profile_dither_spin.setValue(float(profile.get("dither_deg", 0.05)))
        self.profile_readout_spin.setValue(float(profile.get("readout_seconds", 0.0)))
        self.profile_slew_spin.setValue(float(profile.get("slew_seconds", 0.0)))
        self.plan_profile_combo.setCurrentIndex(index)

    def save_equipment_profile(self):
        name = self.profile_name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Profil sestavy", "Zadejte název profilu.")
            return
        profile = {
            "name": name,
            "width_deg": self.profile_width_spin.value(),
            "height_deg": self.profile_height_spin.value(),
            "rotation_deg": self.profile_rotation_spin.value(),
            "dither_deg": self.profile_dither_spin.value(),
            "readout_seconds": self.profile_readout_spin.value(),
            "slew_seconds": self.profile_slew_spin.value(),
        }
        index = self.equipment_combo.currentIndex()
        previous_name = self.planner_store.data["profiles"][index]["name"]
        duplicate = next((i for i, item in enumerate(self.planner_store.data["profiles"])
                          if item["name"].casefold() == name.casefold() and i != index), None)
        if duplicate is not None:
            QMessageBox.warning(self, "Profil sestavy", "Profil s tímto názvem už existuje.")
            return
        self.planner_store.data["profiles"][index] = profile
        for field in self.planned_fields:
            field_profile = field.get("profile", {})
            if field_profile.get("name") == previous_name:
                field_profile["name"] = name
                field_profile["readout_seconds"] = profile["readout_seconds"]
                field_profile["slew_seconds"] = profile["slew_seconds"]
        self.planner_store.save()
        self.refresh_equipment_ui()
        self.equipment_combo.setCurrentIndex(index)
        self.refresh_plan_table()
        self.redraw()
        self.status_label.setText(f"Profil „{name}“ byl uložen.")

    def new_equipment_profile(self):
        number = len(self.planner_store.data["profiles"]) + 1
        self.planner_store.data["profiles"].append({
            "name": f"Nová sestava {number}",
            "width_deg": 10.0,
            "height_deg": 7.0,
            "rotation_deg": 0.0,
            "dither_deg": 0.05,
            "readout_seconds": 0.0,
            "slew_seconds": 0.0,
        })
        self.planner_store.save()
        self.refresh_equipment_ui()
        last = len(self.planner_store.data["profiles"]) - 1
        self.equipment_combo.setCurrentIndex(last)
        self.profile_name_edit.selectAll()
        self.profile_name_edit.setFocus()

    def delete_equipment_profile(self):
        profiles = self.planner_store.data["profiles"]
        if len(profiles) == 1:
            QMessageBox.information(self, "Profil sestavy", "Alespoň jeden profil musí zůstat zachován.")
            return
        index = self.equipment_combo.currentIndex()
        name = profiles[index]["name"]
        if QMessageBox.question(
                self, "Odstranit profil", f"Opravdu odstranit profil „{name}“?") != QMessageBox.StandardButton.Yes:
            return
        del profiles[index]
        self.planner_store.save()
        self.refresh_equipment_ui()

    def map_position_to_icrs(self, x, y):
        qdt = self.time_edit.dateTime().toPyDateTime()
        frame_altaz = AltAz(
            obstime=Time(qdt), location=OBSERVATORIES[self.obs_combo.currentText()]
        )
        if self.view_mode_combo.currentIndex() == 1:
            projection = self.projection_combo.currentIndex()
            if projection == 1 and x ** 2 / 8.0 + y ** 2 / 2.0 > 1.0:
                return None
            local_x, altitude = self.inverse_local(x, y)
            if projection == 1:
                azimuth = (local_x + (self.local_sun_azimuth or 0)) % 360
            else:
                azimuth = (local_x + (self.local_sun_azimuth or 0) - 180) % 360
            azimuth, altitude = float(np.asarray(azimuth)), float(np.asarray(altitude))
            if not np.isfinite(azimuth) or not np.isfinite(altitude) or not -90 <= altitude <= 90:
                return None
            return SkyCoord(az=azimuth * u.deg, alt=altitude * u.deg,
                            frame=frame_altaz).transform_to("icrs")
        if self.projection_combo.currentIndex() == 1 and x ** 2 / 8.0 + y ** 2 / 2.0 > 1.0:
            return None
        ra, dec = self.inverse_global(x, y)
        ra, dec = float(np.asarray(ra)), float(np.asarray(dec))
        if not np.isfinite(ra) or not np.isfinite(dec) or not -90 <= dec <= 90:
            return None
        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    def on_mosaic_orientation_changed(self, index):
        if index == 1:
            self.mosaic_overlap_label.setText("Min. překryv")
            self.status_label.setText(
                "Mozaika bude vodorovná s obzorem v čase mapy; rozteč zohlední "
                "naklonění EQ kamery."
            )
        else:
            self.mosaic_overlap_label.setText("Překryv")
        self.mosaic_overlap_spin.setSuffix(" %")
        self.clear_field_preview()

    @staticmethod
    def horizon_offset_coordinate(origin_altaz, horizontal_deg, vertical_deg):
        """Offset in local horizon axes while keeping a row at constant altitude."""
        altitude = float(origin_altaz.alt.deg) + float(vertical_deg)
        altitude = float(np.clip(altitude, -89.8, 89.8))
        cos_altitude = max(0.03, abs(math.cos(math.radians(altitude))))
        azimuth = (float(origin_altaz.az.deg) + float(horizontal_deg) / cos_altitude) % 360.0
        return SkyCoord(az=azimuth * u.deg, alt=altitude * u.deg,
                        frame=origin_altaz.frame)

    def eq_camera_angle_to_horizon(self, coordinate, profile, frame_altaz):
        """Return sensor X-axis angle from the local horizontal direction."""
        east, north = rotated_offset(0.05, 0.0, float(profile["rotation_deg"]))
        endpoint_ra, endpoint_dec = offset_coordinate(
            float(coordinate.icrs.ra.deg), float(coordinate.icrs.dec.deg), east, north
        )
        centre = coordinate.transform_to(frame_altaz)
        endpoint = SkyCoord(
            ra=endpoint_ra * u.deg, dec=endpoint_dec * u.deg, frame="icrs"
        ).transform_to(frame_altaz)
        delta_azimuth = (
            float(endpoint.az.deg) - float(centre.az.deg) + 180.0
        ) % 360.0 - 180.0
        horizontal = delta_azimuth * math.cos(math.radians(float(centre.alt.deg)))
        vertical = float(endpoint.alt.deg) - float(centre.alt.deg)
        return math.degrees(math.atan2(vertical, horizontal))

    def make_planner_mosaic(self, coordinate, frame_altaz):
        profile = self.current_profile()
        rows = self.mosaic_rows_spin.value()
        columns = self.mosaic_columns_spin.value()
        overlap = self.mosaic_overlap_spin.value()
        if self.mosaic_orientation_combo.currentIndex() == 0:
            return make_mosaic(
                float(coordinate.icrs.ra.deg), float(coordinate.icrs.dec.deg), profile,
                rows, columns, overlap,
            )

        origin_altaz = coordinate.transform_to(frame_altaz)
        width = float(profile["width_deg"])
        height = float(profile["height_deg"])
        centre_angle = self.eq_camera_angle_to_horizon(coordinate, profile, frame_altaz)
        horizontal_step, vertical_step = safe_horizon_grid_steps(
            width, height, centre_angle, overlap
        )

        # A wide mosaic sees a changing parallactic angle.  Sample the full
        # footprint and use the smallest safe spacing found there.
        for _ in range(2):
            sampled_horizontal = [horizontal_step]
            sampled_vertical = [vertical_step]
            for row_index in range(rows):
                y = ((rows - 1) / 2.0 - row_index) * vertical_step
                for column_index in range(columns):
                    x = (column_index - (columns - 1) / 2.0) * horizontal_step
                    sample = self.horizon_offset_coordinate(origin_altaz, x, y).transform_to("icrs")
                    angle = self.eq_camera_angle_to_horizon(sample, profile, frame_altaz)
                    step_x, step_y = safe_horizon_grid_steps(width, height, angle, overlap)
                    sampled_horizontal.append(step_x)
                    sampled_vertical.append(step_y)
            horizontal_step = min(sampled_horizontal)
            vertical_step = min(sampled_vertical)

        result = []
        for row_index in range(rows):
            y = ((rows - 1) / 2.0 - row_index) * vertical_step
            for column_index in range(columns):
                x = (column_index - (columns - 1) / 2.0) * horizontal_step
                field_coordinate = self.horizon_offset_coordinate(
                    origin_altaz, x, y
                ).transform_to("icrs")
                result.append({
                    "ra_deg": float(field_coordinate.ra.deg),
                    "dec_deg": float(field_coordinate.dec.deg),
                    "profile": dict(profile),
                    "mosaic_orientation": "horizon_eq",
                    "requested_overlap_percent": float(overlap),
                    "horizontal_step_deg": float(horizontal_step),
                    "vertical_step_deg": float(vertical_step),
                })
        return result

    def on_map_click(self, event):
        if (not self.placement_check.isChecked()
                or self.workspace_tabs.currentIndex() != 1
                or event.button() != Qt.MouseButton.LeftButton):
            return
        view_box = self.plot_widget.plotItem.vb
        if not view_box.sceneBoundingRect().contains(event.scenePos()):
            return
        point = view_box.mapSceneToView(event.scenePos())
        coordinate = self.map_position_to_icrs(point.x(), point.y())
        if coordinate is None:
            self.status_label.setText("Pole nelze umístit mimo platnou oblast projekce.")
            return
        qdt = self.time_edit.dateTime().toPyDateTime()
        frame_altaz = AltAz(
            obstime=Time(qdt), location=OBSERVATORIES[self.obs_combo.currentText()]
        )
        new_fields = self.make_planner_mosaic(coordinate, frame_altaz)
        self.planned_fields.extend(new_fields)
        self.refresh_plan_table()
        self.redraw()
        event.accept()
        detail = ""
        if new_fields and new_fields[0].get("mosaic_orientation") == "horizon_eq":
            detail = (
                f" Vodorovná rozteč {new_fields[0]['horizontal_step_deg']:.2f}°, "
                f"svislá {new_fields[0]['vertical_step_deg']:.2f}°."
            )
        self.status_label.setText(
            f"Přidáno {len(new_fields)} polí; plán nyní obsahuje "
            f"{len(self.planned_fields)} polí.{detail}"
        )

    def on_placement_mode_changed(self, state):
        if state != Qt.CheckState.Checked.value:
            self.clear_field_preview()

    def clear_field_preview(self):
        for item in getattr(self, "preview_items", []):
            try:
                self.plot_widget.removeItem(item)
            except (RuntimeError, ValueError):
                pass
        self.preview_items = []
        self.preview_fields = []

    def update_field_preview(self, coordinate, frame_altaz, is_local):
        self.clear_field_preview()
        self.preview_fields = self.make_planner_mosaic(coordinate, frame_altaz)
        line_x, line_y = [], []
        for field in self.preview_fields:
            vertices = field_outline(field)
            ra = np.asarray([point[0] for point in vertices])
            dec = np.asarray([point[1] for point in vertices])
            if is_local:
                azimuth, altitude = self.transform_coords(ra, dec, frame_altaz)
                px, py = self.project_local(self.local_azimuth(azimuth), altitude)
            else:
                px, py = self.project_global(ra, dec)
            for edge_index in range(len(px) - 1):
                values = (px[edge_index], py[edge_index], px[edge_index + 1], py[edge_index + 1])
                if not all(np.isfinite(value) for value in values):
                    continue
                if self.projection_combo.currentIndex() in (1, 3, 4):
                    append_projected_segment(
                        line_x, line_y, *values,
                        max_step=self.projected_segment_limit(0.8),
                    )
                else:
                    append_sky_segment(line_x, line_y, *values)
        if line_x:
            preview = pg.PlotCurveItem(
                x=np.asarray(line_x), y=np.asarray(line_y), connect="finite",
                pen=pg.mkPen(
                    color=(255, 190, 80, 220), width=1.8,
                    style=Qt.PenStyle.DashLine,
                ),
            )
            preview.setZValue(40)
            self.plot_widget.addItem(preview)
            self.preview_items.append(preview)

    def refresh_plan_table(self):
        self.fields_table.setRowCount(len(self.planned_fields))
        for row, field in enumerate(self.planned_fields):
            coordinate = SkyCoord(ra=field["ra_deg"] * u.deg, dec=field["dec_deg"] * u.deg)
            values = [
                str(row + 1),
                coordinate.ra.to_string(unit=u.hour, sep=":", precision=1, pad=True),
                coordinate.dec.to_string(unit=u.deg, sep=":", precision=0, alwayssign=True, pad=True),
                field["profile"]["name"],
            ]
            for column, value in enumerate(values):
                self.fields_table.setItem(row, column, QTableWidgetItem(value))
        self.refresh_exposure_script()

    def refresh_exposure_script(self, *_):
        if not hasattr(self, "script_output"):
            return
        duration = estimate_sequence_duration(
            self.planned_fields,
            self.exposure_spin.value(),
            self.loops_spin.value(),
        )
        if self.planned_fields:
            self.sequence_duration_label.setText(
                f"Odhad délky sekvence: {format_duration(duration['total_seconds'])} "
                f"· {duration['frame_count']} snímků\n"
                f"Expozice {format_duration(duration['exposure_seconds'])} · "
                f"vyčítání {format_duration(duration['readout_seconds'])} · "
                f"přejezdy {format_duration(duration['slew_seconds'])}"
            )
        else:
            self.sequence_duration_label.setText("Odhad délky sekvence: --")
        self.script_output.setPlainText(generate_exposure_script(
            self.planned_fields, self.filter_edit.text(),
            self.exposure_spin.value(), self.loops_spin.value(),
        ))

    def undo_planned_field(self):
        if self.planned_fields:
            self.planned_fields.pop()
            self.refresh_plan_table()
            self.redraw()

    def clear_planned_fields(self):
        self.planned_fields.clear()
        self.refresh_plan_table()
        self.redraw()

    def archive_current_plan(self):
        if not self.planned_fields:
            QMessageBox.information(self, "Archiv polí", "Nejprve přidejte pole do mapy.")
            return
        qdt = self.time_edit.dateTime().toPyDateTime()
        eval_time = Time(qdt)
        location = OBSERVATORIES[self.obs_combo.currentText()]
        before_time = eval_time - 5 * u.min
        after_time = eval_time + 5 * u.min
        sun_before = get_sun(before_time).transform_to(
            AltAz(obstime=before_time, location=location)
        )
        sun_after = get_sun(after_time).transform_to(
            AltAz(obstime=after_time, location=location)
        )
        is_morning = float(sun_after.alt.deg) > float(sun_before.alt.deg)
        profile_names = [
            field.get("profile", {}).get("name", "Sestava")
            for field in self.planned_fields
        ]
        default_name = suggested_archive_name(qdt, profile_names, is_morning)

        existing_names = {
            str(item.get("name", "")).casefold()
            for item in self.planner_store.data["archives"]
        }
        base_name = default_name
        suffix = 2
        while default_name.casefold() in existing_names or os.path.exists(
            os.path.join(self.app_dir, "observation_plans", default_name + ".txt")
        ):
            default_name = f"{base_name}_{suffix}"
            suffix += 1
        name, accepted = QInputDialog.getText(self, "Uložit plán", "Název pozorování:", text=default_name)
        if not accepted or not name.strip():
            return
        if name.strip().casefold() in existing_names:
            QMessageBox.warning(
                self, "Uložit plán",
                "Pozorování s tímto názvem už v archivu existuje. Zvolte jiný název."
            )
            return
        archive = new_archive(
            name.strip(), self.planned_fields, self.filter_edit.text(),
            self.exposure_spin.value(), self.loops_spin.value(),
        )
        observation_time = qdt.strftime("%Y-%m-%d %H:%M")
        try:
            text_path = export_plan_text(
                os.path.join(self.app_dir, "observation_plans"), archive,
                observation_time, self.obs_combo.currentText(),
            )
        except FileExistsError:
            QMessageBox.warning(
                self, "Uložit plán",
                "Textový soubor s tímto názvem už existuje. Zvolte jiný název."
            )
            return
        except OSError as error:
            QMessageBox.warning(
                self, "Uložit plán", f"Textový soubor se nepodařilo uložit:\n{error}"
            )
            return
        self.planner_store.data["archives"].append(archive)
        self.planner_store.save()
        self.refresh_archive_table()
        self.redraw()
        self.status_label.setText(
            f"Plán „{name.strip()}“ byl uložen do archivu i do {text_path}."
        )

    def refresh_archive_table(self):
        self.refreshing_planner_ui = True
        archives = self.planner_store.data["archives"]
        self.archive_table.setRowCount(len(archives))
        for row, archive in enumerate(archives):
            visible_item = QTableWidgetItem("")
            visible_item.setFlags(visible_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            visible_item.setCheckState(Qt.CheckState.Checked if archive.get("visible", True)
                                       else Qt.CheckState.Unchecked)
            self.archive_table.setItem(row, 0, visible_item)
            profiles = sorted({field.get("profile", {}).get("name", "?") for field in archive.get("fields", [])})
            values = [archive.get("name", "Bez názvu"), archive.get("created_utc", ""),
                      str(len(archive.get("fields", []))), ", ".join(profiles)]
            for column, value in enumerate(values, start=1):
                self.archive_table.setItem(row, column, QTableWidgetItem(value))
        self.refreshing_planner_ui = False

    def archive_visibility_changed(self, item):
        if self.refreshing_planner_ui or item.column() != 0:
            return
        row = item.row()
        archives = self.planner_store.data["archives"]
        if 0 <= row < len(archives):
            archives[row]["visible"] = item.checkState() == Qt.CheckState.Checked
            self.planner_store.save()
            self.redraw()

    def selected_archive_index(self):
        rows = self.archive_table.selectionModel().selectedRows()
        return rows[0].row() if rows else -1

    def load_selected_archive(self):
        index = self.selected_archive_index()
        if index < 0:
            QMessageBox.information(self, "Archiv polí", "Vyberte řádek archivu.")
            return
        archive = self.planner_store.data["archives"][index]
        self.planned_fields = json.loads(json.dumps(archive.get("fields", [])))
        self.filter_edit.setText(str(archive.get("filter", "N")))
        self.exposure_spin.setValue(float(archive.get("exposure_seconds", 30)))
        self.loops_spin.setValue(int(archive.get("loops", 15)))
        self.refresh_plan_table()
        self.workspace_tabs.setCurrentIndex(1)
        self.planner_tabs.setCurrentIndex(0)
        self.redraw()

    def delete_selected_archive(self):
        index = self.selected_archive_index()
        if index < 0:
            QMessageBox.information(self, "Archiv polí", "Vyberte řádek archivu.")
            return
        archive = self.planner_store.data["archives"][index]
        if QMessageBox.question(
                self, "Odstranit z archivu",
                f"Opravdu odstranit plán „{archive.get('name', 'Bez názvu')}“?") != QMessageBox.StandardButton.Yes:
            return
        del self.planner_store.data["archives"][index]
        self.planner_store.save()
        self.refresh_archive_table()
        self.redraw()

    def update_skycoverage_data(self):
        self.update_data_button.setEnabled(False)
        self.btn_load.setEnabled(False)
        self.skycov_updater = SkyCoverageUpdater(self.selected_dir)
        self.skycov_updater.progress.connect(self.status_label.setText)
        self.skycov_updater.finished.connect(self.on_skycoverage_updated)
        self.skycov_updater.start()

    def on_skycoverage_updated(self, successful, message):
        self.update_data_button.setEnabled(True)
        self.btn_load.setEnabled(True)
        self.status_label.setText(message)
        if successful:
            QMessageBox.information(self, "Aktualizace MPC", message)
            self.start_processing()
        else:
            QMessageBox.warning(self, "Aktualizace MPC", message)

    def set_twilight(self, target_alt_deg, evening=True):
        """Nastaví dnešní večerní nebo následující ranní soumrak v UTC."""
        obs_location = OBSERVATORIES[self.obs_combo.currentText()]

        base_date = datetime.now(timezone.utc).replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        if evening:
            search_date = base_date
            start_hour, end_hour = 12, 30
        else:
            search_date = base_date + timedelta(days=1)
            start_hour, end_hour = 0, 18
        times = [search_date + timedelta(minutes=m) for m in range(start_hour * 60, end_hour * 60 + 1)]

        astropy_times = Time(times)
        frame = AltAz(obstime=astropy_times, location=obs_location)
        suns = get_sun(astropy_times).transform_to(frame)
        alts = suns.alt.deg

        if evening:
            crossings = np.where(
                (alts[:-1] > target_alt_deg) & (alts[1:] <= target_alt_deg)
            )[0]
        else:
            crossings = np.where(
                (alts[:-1] < target_alt_deg) & (alts[1:] >= target_alt_deg)
            )[0]

        if len(crossings):
            best_idx = int(crossings[0])
            if abs(alts[best_idx + 1] - target_alt_deg) < abs(alts[best_idx] - target_alt_deg):
                best_idx += 1
            best_dt = times[best_idx]
        else:
            best_idx = int(np.argmin(np.abs(alts - target_alt_deg)))
            best_dt = times[best_idx]
        best_dt = best_dt.replace(tzinfo=timezone.utc)

        self.time_edit.blockSignals(True)
        self.time_edit.setDateTime(
            QDateTime.fromSecsSinceEpoch(int(best_dt.timestamp()), Qt.TimeSpec.UTC)
        )
        self.time_edit.blockSignals(False)
        self.view_mode_combo.setCurrentIndex(1)
        self.projection_combo.setCurrentIndex(4)
        self.local_mercator_center_azimuth = float(suns[best_idx].az.deg)
        self.local_mercator_center_altitude = 35.0
        self.local_mercator_observatory = self.obs_combo.currentText()
        self.local_sun_azimuth = float(suns[best_idx].az.deg)
        self.mercator_center_azimuth = 180.0
        self.mercator_center_altitude = 35.0
        self.set_default_mer_range()
        self.redraw()
        self.request_data_reload()

    def shift_time(self, minutes):
        self.time_edit.setDateTime(
            self.time_edit.dateTime().addSecs(int(minutes * 60))
        )

    def on_workspace_tab_changed(self, index):
        """Move the single live map between the coverage and planning workspaces."""
        self.clear_field_preview()
        if index == 0:
            self.map_workspace_layout.addWidget(self.plot_widget)
        else:
            self.planning_map_layout.addWidget(self.plot_widget)

    def zoom_map(self, factor):
        self.plot_widget.getViewBox().scaleBy((factor, factor))
        self.request_axis_tick_update()
        self.request_bayer_label_update()

    def reset_map_view(self):
        if (self.view_mode_combo.currentIndex() == 1
                and self.projection_combo.currentIndex() == 4):
            self.local_mercator_observatory = None
        self.on_view_change()

    def set_default_mer_range(self):
        half_width = MER_DEFAULT_HORIZONTAL_FOV_DEG / 2.0
        view_box = self.plot_widget.getViewBox()
        bounds = view_box.sceneBoundingRect()
        aspect = bounds.width() / max(bounds.height(), 1.0)
        half_height = half_width / max(aspect, 0.1)
        view_box.setRange(
            xRange=(-half_width, half_width),
            yRange=(-half_height, half_height),
            padding=0,
            disableAutoRange=True,
        )

    def set_mercator_center(self, coordinate):
        coordinate = coordinate.transform_to("icrs")
        self.mercator_center_ra = float(coordinate.ra.deg)
        self.mercator_center_dec = float(coordinate.dec.deg)
        qdt = self.time_edit.dateTime().toPyDateTime()
        frame_altaz = AltAz(
            obstime=Time(qdt), location=OBSERVATORIES[self.obs_combo.currentText()]
        )
        sun_altaz = get_sun(Time(qdt)).transform_to(frame_altaz)
        center_altaz = coordinate.transform_to(frame_altaz)
        self.mercator_center_azimuth = (
            float(center_altaz.az.deg) - float(sun_altaz.az.deg) + 180.0
        ) % 360.0
        self.mercator_center_altitude = float(center_altaz.alt.deg)

    def on_view_change(self):
        is_local = self.view_mode_combo.currentIndex() == 1
        projection = self.projection_combo.currentIndex()
        self.plot_widget.clear()
        self.plot_widget.getViewBox().setAspectLocked(True, ratio=1.0)
        if is_local:
            self.plot_widget.getViewBox().invertX(False)
            if projection == 1:
                self.plot_widget.setXRange(-3, 3)
                self.plot_widget.setYRange(-1.5, 1.5)
                self.plot_widget.setLabel('bottom', 'Relativní azimut (střed = Slunce)')
                self.plot_widget.setLabel('left', 'Výška Alt (stupně)')
            elif projection == 3:
                self.plot_widget.setXRange(-3.7, 3.7)
                self.plot_widget.setYRange(-3.7, 3.7)
                self.plot_widget.setLabel('bottom', 'Stereografická obloha · zenit uprostřed')
                self.plot_widget.setLabel('left', '')
            elif projection == 4:
                self.set_default_mer_range()
                self.plot_widget.setLabel('bottom', 'MER Mercatorova · oblique')
                self.plot_widget.setLabel('left', '')
            else:
                self.plot_widget.setXRange(0, 360)
                upper_y = LAMBERT_Y_SCALE if projection == 2 else 90.0
                self.plot_widget.setYRange(0, upper_y)
                self.plot_widget.setLabel('bottom', 'Azimut vůči Slunci (střed = směr Slunce)')
                self.plot_widget.setLabel('left', 'Výška nad obzorem Alt (deg)')
        else:
            self.plot_widget.getViewBox().invertX(True)
            if projection == 1:
                self.plot_widget.setXRange(-3, 3)
                self.plot_widget.setYRange(-1.5, 1.5)
                self.plot_widget.setLabel('bottom', 'Rektascenze RA (hodiny)')
                self.plot_widget.setLabel('left', 'Deklinace Dec (stupně)')
            elif projection == 3:
                self.plot_widget.setXRange(-STEREO_VIEW_RADIUS, STEREO_VIEW_RADIUS)
                self.plot_widget.setYRange(-STEREO_VIEW_RADIUS, STEREO_VIEW_RADIUS)
                self.plot_widget.setLabel('bottom', 'Stereografická · střed RA 12h, Dec 0°')
                self.plot_widget.setLabel('left', '')
            elif projection == 4:
                self.set_default_mer_range()
                self.plot_widget.setLabel('bottom', 'MER Mercatorova · oblique')
                self.plot_widget.setLabel('left', '')
            else:
                self.plot_widget.setXRange(0, 360)
                upper_y = LAMBERT_Y_SCALE if projection == 2 else 90.0
                self.plot_widget.setYRange(-upper_y, upper_y)
                self.plot_widget.setLabel('bottom', 'Rektascenze RA (hodiny)')
                self.plot_widget.setLabel('left', 'Deklinace Dec (stupně)')
        self.redraw()
        self.update_axis_ticks()

    def on_projection_change(self, index=None):
        if self.projection_combo.currentIndex() == 4:
            if self.last_cursor_icrs is not None:
                self.set_mercator_center(self.last_cursor_icrs)
            else:
                self.set_mercator_center(get_sun(Time(self.time_edit.dateTime().toPyDateTime())))
        self.on_view_change()

    def update_axis_ticks(self):
        is_local = self.view_mode_combo.currentIndex() == 1
        projection = self.projection_combo.currentIndex()
        hammer = projection == 1
        lambert = projection == 2
        stereo = projection == 3
        mercator = projection == 4
        bottom = self.plot_widget.getAxis('bottom')
        left = self.plot_widget.getAxis('left')
        if stereo or mercator:
            if is_local:
                bottom.setLabel(
                    'Stereografická obloha · zenit uprostřed' if stereo
                    else 'MER Mercatorova · oblique'
                )
            else:
                bottom.setLabel(
                    'Stereografická · střed RA 12h, Dec 0°' if stereo
                    else 'MER Mercatorova · oblique'
                )
            left.setLabel('')
            bottom.setTicks([[]])
            left.setTicks([[]])
            return
        if is_local:
            if hammer:
                bottom_label = 'Relativní azimut (Hammer-Aitoff)'
            elif lambert:
                bottom_label = 'Relativní azimut (válcová plochojevná)'
            else:
                bottom_label = 'Relativní azimut (střed = Slunce)'
            bottom.setLabel(bottom_label)
            left.setLabel('Výška Alt (stupně)')
            x_range, y_range = self.plot_widget.getViewBox().viewRange()
            x_width = abs(x_range[1] - x_range[0])
            y_height = abs(y_range[1] - y_range[0])
            if hammer:
                tick_step = 10 if x_width < 4.5 or y_height < 2.25 else 30
                center_x = (x_range[0] + x_range[1]) / 2.0
                center_y = (y_range[0] + y_range[1]) / 2.0
                if center_x ** 2 / 8.0 + center_y ** 2 / 2.0 <= 1.0:
                    center_longitude, center_altitude = self.inverse_local(
                        center_x, center_y
                    )
                    reference_altitude = float(np.clip(center_altitude, -80.0, 80.0))
                    reference_azimuth = (float(center_longitude) + 180.0) % 360.0
                else:
                    reference_altitude = 0.0
                    reference_azimuth = 180.0

                az_values = np.arange(0, 360, tick_step)
                az_pos, _ = self.project_local(
                    self.local_azimuth(az_values),
                    np.full_like(az_values, reference_altitude, dtype=float),
                )
                bottom.setTicks([[
                    (float(position), f'{value}°')
                    for position, value in zip(az_pos, az_values)
                ]])

                alt_values = np.arange(-90, 91, tick_step)
                _, alt_pos = self.project_local(
                    np.full_like(alt_values, reference_azimuth, dtype=float),
                    alt_values,
                )
                left.setTicks([[
                    (float(position), f'{value}°')
                    for position, value in zip(alt_pos, alt_values)
                ]])
            else:
                y_zoom_threshold = 50.0 if lambert else 80.0
                tick_step = 10 if x_width < 240.0 or y_height < y_zoom_threshold else 30
                az_values = np.arange(0, 360, tick_step)
                az_positions = self.local_azimuth(az_values)
                bottom.setTicks([[
                    (float(position), f'{value}°')
                    for position, value in zip(az_positions, az_values)
                ]])
                alt_values = np.arange(-90, 91, tick_step)
                _, alt_positions = self.project_local(
                    np.full_like(alt_values, 180, dtype=float),
                    alt_values,
                )
                left.setTicks([[
                    (float(position), f'{value}°')
                    for position, value in zip(alt_positions, alt_values)
                ]])
            return

        bottom.setLabel('Rektascenze RA (hodiny)')
        left.setLabel('Deklinace Dec (stupně)')
        ra_values = np.arange(0, 360, 15)
        ra_minor = np.arange(0, 360, 7.5)
        ra_positions = self.project_global(ra_values, np.zeros_like(ra_values))[0]
        ra_minor_positions = self.project_global(ra_minor, np.zeros_like(ra_minor))[0]
        dec_values = np.arange(-80, 81, 20)
        dec_minor = np.arange(-90, 91, 5)
        dec_positions = self.project_global(np.full_like(dec_values, 180), dec_values)[1]
        dec_minor_positions = self.project_global(np.full_like(dec_minor, 180), dec_minor)[1]
        bottom.setTicks([
            [(float(position), '') for position in ra_minor_positions],
            [(float(position), f'{int(ra / 15)}h') for position, ra in zip(ra_positions, ra_values)],
        ])
        left.setTicks([
            [(float(position), '') for position in dec_minor_positions],
            [(float(position), f'{int(dec)}°') for position, dec in zip(dec_positions, dec_values)],
        ])

    def on_mouse_move(self, pos):
        view_box = self.plot_widget.plotItem.vb
        if not view_box.sceneBoundingRect().contains(pos):
            return
        mouse_point = view_box.mapSceneToView(pos)
        x, y = mouse_point.x(), mouse_point.y()
        qdt = self.time_edit.dateTime().toPyDateTime()
        eval_time = Time(qdt)
        location = OBSERVATORIES[self.obs_combo.currentText()]
        frame_altaz = AltAz(obstime=eval_time, location=location)
        if self.view_mode_combo.currentIndex() == 1:
            projection = self.projection_combo.currentIndex()
            if projection == 1:
                if x ** 2 / 8.0 + y ** 2 / 2.0 > 1.0:
                    self.show_invalid_coordinates()
                    return
            local_x, altitude = self.inverse_local(x, y)
            if projection == 1:
                azimuth = (local_x + (self.local_sun_azimuth or 0)) % 360
            else:
                azimuth = (local_x + (self.local_sun_azimuth or 0) - 180) % 360
            azimuth = float(np.asarray(azimuth))
            altitude = float(np.asarray(altitude))
            if (not np.isfinite(azimuth) or not np.isfinite(altitude) or
                    altitude < -90.0 or altitude > 90.0):
                self.show_invalid_coordinates()
                return
            target_altaz = SkyCoord(az=azimuth * u.deg, alt=altitude * u.deg, frame=frame_altaz)
            target = target_altaz.transform_to('icrs')
        else:
            azimuth = None
            if self.projection_combo.currentIndex() == 1 and (x ** 2 / 8.0 + y ** 2 / 2.0 > 1.0):
                self.show_invalid_coordinates()
                return
            ra, dec = self.inverse_global(x, y)
            ra = float(np.asarray(ra))
            dec = float(np.asarray(dec))
            if (not np.isfinite(ra) or not np.isfinite(dec) or
                    dec < -90.0 or dec > 90.0):
                self.show_invalid_coordinates()
                return
            target = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

        sun_altaz = get_sun(eval_time).transform_to(frame_altaz)
        self.current_sun_altitude = float(sun_altaz.alt.deg)
        target_altaz = target.transform_to(frame_altaz)
        elongation = target_altaz.separation(sun_altaz).deg
        target = target.transform_to('icrs')
        self.last_cursor_icrs = target
        ra_hours = target.ra.deg / 15.0
        display_azimuth = float(target_altaz.az.deg)
        display_altitude = float(target_altaz.alt.deg)
        self.coord_label.setText(
            f"Az: {display_azimuth:5.1f}° | Alt: {display_altitude:+5.1f}°\n"
            f"RA: {int(ra_hours):02d}h {int((ra_hours % 1) * 60):02d}m "
            f"| Dec: {target.dec.deg:+5.1f}°\n"
            f"Elongace: {elongation:.1f}°\n"
            f"{format_sun_altitude(self.current_sun_altitude)}"
        )
        if self.placement_check.isChecked() and self.workspace_tabs.currentIndex() == 1:
            self.update_field_preview(target, frame_altaz, self.view_mode_combo.currentIndex() == 1)
        else:
            self.clear_field_preview()

    def show_invalid_coordinates(self):
        self.clear_field_preview()
        self.coord_label.setText(
            "Az: -- | Alt: --\n"
            "RA: -- | Dec: --\n"
            "Mimo platnou oblast projekce\n"
            f"{format_sun_altitude(self.current_sun_altitude)}"
        )

    def update_sun_altitude_readout(self, altitude):
        """Refresh the persistent Sun line without discarding cursor details."""
        self.current_sun_altitude = float(altitude)
        lines = [
            line for line in self.coord_label.text().splitlines()
            if not line.startswith("Slunce:")
        ]
        lines.append(format_sun_altitude(self.current_sun_altitude))
        self.coord_label.setText("\n".join(lines))

    def request_redraw(self):
        self.redraw_timer.start(80)

    def request_axis_tick_update(self, *_):
        self.axis_tick_timer.start(50)

    def request_bayer_label_update(self, *_):
        self.bayer_label_timer.start(120)

    def request_mercator_refresh(self, *_):
        if (self.projection_combo.currentIndex() == 4
                and not self.reprojecting_mercator):
            # Qt moves the existing scene immediately while dragging. Refill
            # and reproject it only after the view has settled.
            self.mer_pan_timer.start(300)

    def refresh_mercator_view(self):
        if self.projection_combo.currentIndex() != 4:
            return
        view_box = self.plot_widget.getViewBox()
        x_range, y_range = view_box.viewRange()
        center_x = (x_range[0] + x_range[1]) / 2.0
        center_y = (y_range[0] + y_range[1]) / 2.0
        width = max(float(x_range[1] - x_range[0]), 1e-6)
        height = max(float(y_range[1] - y_range[0]), 1e-6)

        self.reprojecting_mercator = True
        try:
            if self.view_mode_combo.currentIndex() == 1:
                projected_azimuth, altitude = self.inverse_local(center_x, center_y)
                self.local_mercator_center_azimuth = (
                    float(np.asarray(projected_azimuth))
                    + (self.local_sun_azimuth or 0.0) - 180.0
                ) % 360.0
                self.local_mercator_center_altitude = float(np.clip(
                    np.asarray(altitude), -85.0, 85.0
                ))
            else:
                right_ascension, declination = self.inverse_global(center_x, center_y)
                self.mercator_center_ra = float(np.asarray(right_ascension)) % 360.0
                self.mercator_center_dec = float(np.clip(
                    np.asarray(declination), -85.0, 85.0
                ))

            view_box.setRange(
                xRange=(-width / 2.0, width / 2.0),
                yRange=(-height / 2.0, height / 2.0),
                padding=0,
                disableAutoRange=True,
            )
            self.mer_polygon_cache = None
            self.redraw()
        finally:
            self.reprojecting_mercator = False

    def request_data_reload(self):
        self.reload_timer.start(250)

    def start_processing(self):
        self.btn_load.setEnabled(False)
        self.update_data_button.setEnabled(False)
        displayed_time = self.time_edit.dateTime().toPyDateTime()
        self.reader = PointingDataReader(
            self.selected_dir,
            self.days_spin.value(),
            displayed_time,
        )
        self.reader.progress.connect(self.status_label.setText)
        self.reader.finished.connect(self.on_data_loaded)
        self.reader.start()

    def on_data_loaded(self, result):
        self.polygons, self.ref_time = result
        self.uncovered_cache.clear()
        self.btn_load.setEnabled(True)
        self.update_data_button.setEnabled(True)
        if self.ref_time is None:
            self.data_age_label.setText("Nejnovější MPC data: nejsou k dispozici")
        else:
            newest = self.ref_time.to_datetime(timezone=timezone.utc)
            age_days = max(0, (datetime.now(timezone.utc).date() - newest.date()).days)
            warning = " — AKTUALIZUJTE" if age_days > 7 else ""
            self.data_age_label.setText(
                f"Nejnovější MPC data: {newest:%Y-%m-%d} (stáří {age_days} dní){warning}"
            )
        self.redraw()

    def on_stars_loaded(self, result):
        self.stars, error = result
        if error:
            self.status_label.setText(error)
        self.request_redraw()

    def on_supplement_catalogs_loaded(self, result):
        self.deep_sky_objects, self.comet_elements, warning = result
        self.comet_position_cache = None
        if warning:
            self.status_label.setText(warning)
        else:
            self.status_label.setText(
                f"Katalogy: {len(self.deep_sky_objects)} DSO, "
                f"{len(self.comet_elements)} kometárních drah."
            )
        self.request_redraw()

    def transform_coords(self, ra_arr, dec_arr, frame_altaz):
        sc = SkyCoord(ra=np.atleast_1d(ra_arr)*u.deg, dec=np.atleast_1d(dec_arr)*u.deg)
        altaz = sc.transform_to(frame_altaz)
        return altaz.az.deg, altaz.alt.deg

    def local_azimuth(self, azimuth):
        if self.local_sun_azimuth is None:
            return np.asarray(azimuth)
        return (np.asarray(azimuth) - self.local_sun_azimuth + 180.0) % 360.0

    def projected_segment_limit(self, default=0.35):
        """Return a seam-rejection limit in the active projection's units."""
        return 45.0 if self.projection_combo.currentIndex() == 4 else default

    def approximate_horizontal_fov(self):
        x_range = self.plot_widget.getViewBox().viewRange()[0]
        span = abs(float(x_range[1] - x_range[0]))
        projection = self.projection_combo.currentIndex()
        if projection in (0, 2, 4):
            return span
        if projection == 1:
            return span * 60.0
        return span * 180.0 / (2.0 * STEREO_VIEW_RADIUS)

    def visible_star_magnitude(self):
        field_width = self.approximate_horizontal_fov()
        if field_width <= 55.0:
            return 7.0
        if field_width <= 85.0:
            return 6.0
        return 5.0

    @staticmethod
    def preferred_local_mercator_center(location):
        """Face the celestial equator while keeping the local horizon below it."""
        latitude = float(location.lat.to_value(u.deg))
        facing_azimuth = 180.0 if latitude >= 0.0 else 0.0
        return facing_azimuth, 45.0

    def project_global(self, ra, dec):
        ra = np.asarray(ra, dtype=float)
        dec = np.asarray(dec, dtype=float)
        projection = self.projection_combo.currentIndex()
        if projection == 0:
            return ra % 360.0, dec
        if projection == 2:
            return ra % 360.0, LAMBERT_Y_SCALE * np.sin(np.radians(dec))
        if projection == 4:
            return self.project_oblique_mercator(
                ra, dec, self.mercator_center_ra, self.mercator_center_dec
            )

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
            y = np.where(safe, 2.0 * np.sin(latitude) / safe_denominator, np.nan)
            visible = np.hypot(x, y) <= STEREO_VIEW_RADIUS * 1.02
            return np.where(visible, x, np.nan), np.where(visible, y, np.nan)

        denominator = np.sqrt(1.0 + np.cos(latitude) * np.cos(longitude / 2.0))
        return (
            2.0 * np.sqrt(2.0) * np.cos(latitude) * np.sin(longitude / 2.0) / denominator,
            np.sqrt(2.0) * np.sin(latitude) / denominator,
        )

    def inverse_global(self, x, y):
        projection = self.projection_combo.currentIndex()
        if projection == 0:
            return np.asarray(x), np.asarray(y)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if projection == 2:
            valid = np.abs(y) <= LAMBERT_Y_SCALE
            dec = np.where(
                valid,
                np.degrees(np.arcsin(np.clip(y / LAMBERT_Y_SCALE, -1.0, 1.0))),
                np.nan,
            )
            return x % 360.0, dec
        if projection == 4:
            return self.inverse_oblique_mercator(
                x, y, self.mercator_center_ra, self.mercator_center_dec
            )
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
                np.where(valid, (np.degrees(longitude) + 180.0) % 360.0, np.nan),
                np.where(valid, np.degrees(dec), np.nan),
            )
        z = np.sqrt(np.maximum(0.0, 1.0 - (x ** 2 / 16.0) - (y ** 2 / 4.0)))
        longitude = 2.0 * np.arctan2(z * x, 2.0 * (2.0 * z ** 2 - 1.0))
        latitude = np.arcsin(np.clip(z * y, -1.0, 1.0))
        return (np.degrees(longitude) + 180.0) % 360.0, np.degrees(latitude)

    def project_local(self, azimuth, altitude):
        azimuth = np.asarray(azimuth, dtype=float)
        altitude = np.asarray(altitude, dtype=float)
        projection = self.projection_combo.currentIndex()
        if projection == 0:
            return azimuth, altitude
        if projection == 2:
            return azimuth, LAMBERT_Y_SCALE * np.sin(np.radians(altitude))
        if projection == 4:
            return self.project_oblique_mercator(
                azimuth, altitude,
                self.mercator_center_azimuth, self.mercator_center_altitude,
            )
        if projection == 3:
            angle = np.radians((azimuth - 180.0 + 180.0) % 360.0 - 180.0)
            radius = 2.0 * np.tan(np.radians((90.0 - altitude) / 2.0))
            return radius * np.sin(angle), -radius * np.cos(angle)
        longitude = np.radians((azimuth - 180.0 + 180.0) % 360.0 - 180.0)
        latitude = np.radians(altitude)
        denominator = np.sqrt(1.0 + np.cos(latitude) * np.cos(longitude / 2.0))
        return (
            2.0 * np.sqrt(2.0) * np.cos(latitude) * np.sin(longitude / 2.0) / denominator,
            np.sqrt(2.0) * np.sin(latitude) / denominator,
        )

    def inverse_local(self, x, y):
        projection = self.projection_combo.currentIndex()
        if projection == 0:
            return np.asarray(x), np.asarray(y)
        if projection == 2:
            x = np.asarray(x, dtype=float)
            y = np.asarray(y, dtype=float)
            valid = np.abs(y) <= LAMBERT_Y_SCALE
            altitude = np.where(
                valid,
                np.degrees(np.arcsin(np.clip(y / LAMBERT_Y_SCALE, -1.0, 1.0))),
                np.nan,
            )
            return x, altitude

        if projection == 4:
            return self.inverse_oblique_mercator(
                x, y,
                self.mercator_center_azimuth, self.mercator_center_altitude,
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

        z = np.sqrt(np.maximum(0.0, 1.0 - (x ** 2 / 16.0) - (y ** 2 / 4.0)))
        longitude = 2.0 * np.arctan2(z * x, 2.0 * (2.0 * z ** 2 - 1.0))
        latitude = np.arcsin(np.clip(z * y, -1.0, 1.0))
        return np.degrees(longitude), np.degrees(latitude)

    @staticmethod
    def project_oblique_mercator(longitude_deg, latitude_deg,
                                 center_longitude_deg, center_latitude_deg):
        longitude = np.radians(np.asarray(longitude_deg, dtype=float))
        latitude = np.radians(np.asarray(latitude_deg, dtype=float))
        center_longitude = np.radians(float(center_longitude_deg))
        center_latitude = np.radians(float(center_latitude_deg))
        delta_longitude = (longitude - center_longitude + np.pi) % (2 * np.pi) - np.pi
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
        safe_latitude = np.clip(rotated_latitude, np.radians(-89.5), np.radians(89.5))
        return (
            np.degrees(rotated_longitude),
            np.degrees(np.arcsinh(np.tan(safe_latitude))),
        )

    @staticmethod
    def inverse_oblique_mercator(x, y, center_longitude_deg,
                                 center_latitude_deg):
        rotated_longitude = np.radians(np.asarray(x, dtype=float))
        rotated_latitude = np.arctan(np.sinh(np.radians(np.asarray(y, dtype=float))))
        center_longitude = np.radians(float(center_longitude_deg))
        center_latitude = np.radians(float(center_latitude_deg))
        east = np.cos(rotated_latitude) * np.sin(rotated_longitude)
        north = np.sin(rotated_latitude)
        forward = np.cos(rotated_latitude) * np.cos(rotated_longitude)
        latitude = np.arcsin(np.clip(
            np.cos(center_latitude) * north
            + np.sin(center_latitude) * forward,
            -1.0, 1.0,
        ))
        longitude = center_longitude + np.arctan2(
            east,
            np.cos(center_latitude) * forward
            - np.sin(center_latitude) * north,
        )
        return np.degrees(longitude) % 360.0, np.degrees(latitude)

    @staticmethod
    def comet_heliocentric_ecliptic(element, jd):
        """Approximate a comet's heliocentric J2000 ecliptic vector in AU."""
        gaussian_constant = 0.01720209895
        q = element["q"]
        eccentricity = element["e"]
        elapsed_days = jd - element["tp_jd"]
        if q <= 0.0 or eccentricity < 0.0:
            return None

        if abs(eccentricity - 1.0) < 1e-4:
            barker = gaussian_constant * elapsed_days / np.sqrt(2.0 * q ** 3)
            anomaly = 2.0 * np.sinh(np.arcsinh(1.5 * barker) / 3.0)
            true_anomaly = 2.0 * np.arctan(anomaly)
            radius = q * (1.0 + anomaly ** 2)
        elif eccentricity < 1.0:
            semi_major = q / (1.0 - eccentricity)
            mean_anomaly = gaussian_constant * elapsed_days / semi_major ** 1.5
            mean_anomaly = (mean_anomaly + np.pi) % (2.0 * np.pi) - np.pi
            eccentric_anomaly = mean_anomaly
            if eccentricity > 0.8:
                eccentric_anomaly = np.pi if mean_anomaly >= 0.0 else -np.pi
            for _ in range(30):
                correction = (
                    eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)
                    - mean_anomaly
                ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
                eccentric_anomaly -= correction
                if abs(correction) < 1e-12:
                    break
            radius = semi_major * (1.0 - eccentricity * np.cos(eccentric_anomaly))
            true_anomaly = 2.0 * np.arctan2(
                np.sqrt(1.0 + eccentricity) * np.sin(eccentric_anomaly / 2.0),
                np.sqrt(1.0 - eccentricity) * np.cos(eccentric_anomaly / 2.0),
            )
        else:
            semi_major = q / (eccentricity - 1.0)
            mean_anomaly = gaussian_constant * elapsed_days / semi_major ** 1.5
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
                np.sqrt(eccentricity + 1.0) * np.sinh(hyperbolic_anomaly / 2.0),
                np.sqrt(eccentricity - 1.0) * np.cosh(hyperbolic_anomaly / 2.0),
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

    def current_comet_positions(self, eval_time):
        cache_key = round(float(eval_time.tt.jd), 6)
        if self.comet_position_cache is not None and self.comet_position_cache[0] == cache_key:
            return self.comet_position_cache[1]
        if not self.comet_elements:
            return []

        earth_heliocentric = (
            get_body_barycentric("earth", eval_time).xyz.to_value(u.au)
            - get_body_barycentric("sun", eval_time).xyz.to_value(u.au)
        )
        obliquity = np.radians(23.4392911)
        positions = []
        for element in self.comet_elements:
            result = self.comet_heliocentric_ecliptic(element, eval_time.tt.jd)
            if result is None:
                continue
            ecliptic, radius = result
            equatorial = np.array([
                ecliptic[0],
                ecliptic[1] * np.cos(obliquity) - ecliptic[2] * np.sin(obliquity),
                ecliptic[1] * np.sin(obliquity) + ecliptic[2] * np.cos(obliquity),
            ])
            geocentric = equatorial - earth_heliocentric
            distance = float(np.linalg.norm(geocentric))
            if distance <= 0.0 or radius <= 0.0:
                continue
            magnitude = (
                element["h"] + 5.0 * np.log10(distance)
                + 2.5 * element["slope"] * np.log10(radius)
            )
            if not np.isfinite(magnitude) or magnitude > 17.0:
                continue
            ra = np.degrees(np.arctan2(geocentric[1], geocentric[0])) % 360.0
            dec = np.degrees(np.arcsin(np.clip(geocentric[2] / distance, -1.0, 1.0)))
            positions.append((ra, dec, float(magnitude), element["name"]))
        self.comet_position_cache = (cache_key, positions)
        return positions

    def kreutz_position(self, eval_time, days_before, population, variant):
        """Return a geocentric inbound position for a seasonal Kreutz orbit."""
        elements = {
            # C/2024 S1 (ATLAS), JPL solution 13.
            "I": (0.007953356, 0.999916164, 141.90213, 347.15413, 69.16762),
            # C/1882 R1-A, representative of the second subgroup.
            "II": (0.007750, 0.999899, 142.0112, 347.6563, 69.5843),
        }
        perihelion_distance, eccentricity, inclination, node, argument = elements[population]
        node += variant * 1.8
        argument += variant * 0.9
        perihelion_distance *= 1.0 + variant * 0.025

        reference_perihelion = eval_time + float(days_before) * u.day
        days_from_perihelion = -float(days_before)
        gravitational_parameter = 0.0002959122083
        semi_latus_rectum = perihelion_distance * (1.0 + eccentricity)
        scale = np.sqrt(semi_latus_rectum ** 3 / (2.0 * gravitational_parameter))
        eccentric_anomaly = -max(abs(days_from_perihelion) / max(scale, 1e-12), 0.0) ** (1.0 / 3.0)
        for _ in range(10):
            function = scale * (eccentric_anomaly + eccentric_anomaly ** 3 / 3.0) - days_from_perihelion
            derivative = scale * (1.0 + eccentric_anomaly ** 2)
            eccentric_anomaly -= function / derivative

        true_anomaly = 2.0 * np.arctan(eccentric_anomaly)
        radius = semi_latus_rectum / (1.0 + eccentricity * np.cos(true_anomaly))
        x_orbit = radius * np.cos(true_anomaly)
        y_orbit = radius * np.sin(true_anomaly)

        node_rad, inclination_rad, argument_rad = np.radians([node, inclination, argument])
        cos_node, sin_node = np.cos(node_rad), np.sin(node_rad)
        cos_inc, sin_inc = np.cos(inclination_rad), np.sin(inclination_rad)
        cos_arg, sin_arg = np.cos(argument_rad), np.sin(argument_rad)
        x = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit + (-cos_node * sin_arg - sin_node * cos_arg * cos_inc) * y_orbit
        y = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit + (-sin_node * sin_arg + cos_node * cos_arg * cos_inc) * y_orbit
        z = (sin_arg * sin_inc) * x_orbit + (cos_arg * sin_inc) * y_orbit

        ecliptic_frame = HeliocentricMeanEcliptic(equinox='J2000')
        earth_barycentric = get_body_barycentric('earth', eval_time) - get_body_barycentric('sun', eval_time)
        earth = SkyCoord(
            x=earth_barycentric.x, y=earth_barycentric.y, z=earth_barycentric.z,
            representation_type='cartesian', frame='icrs'
        ).transform_to(ecliptic_frame)
        geo_x = x - earth.cartesian.x.to_value(u.au)
        geo_y = y - earth.cartesian.y.to_value(u.au)
        geo_z = z - earth.cartesian.z.to_value(u.au)
        comet = SkyCoord(
            x=geo_x * u.au, y=geo_y * u.au, z=geo_z * u.au,
            representation_type='cartesian',
            frame=ecliptic_frame
        ).transform_to('icrs')
        return float(comet.ra.deg), float(comet.dec.deg)

    def add_kreutz_corridors(self, eval_time, frame_altaz, is_local):
        if not self.kreutz_check.isChecked():
            return

        styles = {
            "I": ((255, 100, 80, 210), Qt.PenStyle.DashLine),
            "II": ((190, 100, 255, 210), Qt.PenStyle.DotLine),
        }
        days_values = [5, 10, 15, 20, 30, 40, 45]
        marker_days = {5, 10, 20, 30, 40}
        for population, (color, line_style) in styles.items():
            corridor_x, corridor_y = [], []
            central_markers = []
            for variant in np.linspace(-2.0, 2.0, 9):
                points_x, points_y = [], []
                for days_before in days_values:
                    ra, dec = self.kreutz_position(eval_time, days_before, population, variant)
                    if is_local:
                        az, alt = self.transform_coords(ra, dec, frame_altaz)
                        x, y = self.project_local(self.local_azimuth(az), alt)
                    else:
                        x, y = self.project_global([ra], [dec])
                    points_x.append(float(x[0]))
                    points_y.append(float(y[0]))
                    if np.isclose(variant, 0.0) and days_before in marker_days:
                        central_markers.append((days_before, points_x[-1], points_y[-1]))
                for index in range(len(points_x) - 1):
                    if is_local and self.projection_combo.currentIndex() in (1, 3, 4):
                        append_projected_segment(
                            corridor_x, corridor_y,
                            points_x[index], points_y[index],
                            points_x[index + 1], points_y[index + 1],
                            self.projected_segment_limit(0.5),
                        )
                    elif is_local:
                        append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
                    else:
                        append_projected_segment(
                            corridor_x, corridor_y,
                            points_x[index], points_y[index],
                            points_x[index + 1], points_y[index + 1],
                            self.projected_segment_limit(0.5),
                        ) if self.projection_combo.currentIndex() in (1, 3, 4) else append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
            self.add_curve(corridor_x, corridor_y, pg.mkPen(color=color, width=1.5, style=line_style), 6)

            central_markers = [
                point for point in central_markers
                if np.isfinite(point[1]) and np.isfinite(point[2])
            ]
            if central_markers:
                marker_item = pg.ScatterPlotItem(
                    x=[point[1] for point in central_markers],
                    y=[point[2] for point in central_markers],
                    symbol="o",
                    size=8,
                    pen=pg.mkPen((235, 235, 245, 230), width=0.8),
                    brush=pg.mkBrush(color),
                )
                marker_item.setZValue(8)
                self.plot_widget.addItem(marker_item)

            label_anchor = (0.5, 1.15) if population == "I" else (0.5, -0.15)
            for days_before, x, y in central_markers:
                label = pg.TextItem(
                    text=f"{days_before} d",
                    color=color,
                    anchor=label_anchor,
                    border=None,
                    fill=None,
                )
                label.setZValue(9)
                label.setPos(x, y)
                self.plot_widget.addItem(label)

    def add_constellation_boundaries(self, frame_altaz, is_local):
        if not self.const_boundaries_check.isChecked() or not self.const_boundaries:
            return

        boundary_arrays = []
        for boundary in self.const_boundaries:
            coordinates = np.asarray(boundary, dtype=float)
            if coordinates.ndim == 2 and coordinates.shape[0] >= 2 and coordinates.shape[1] >= 2:
                boundary_arrays.append(coordinates[:, :2])

        if not boundary_arrays:
            return

        transformed = None
        if is_local:
            all_ra = np.concatenate([coordinates[:, 0] % 360.0 for coordinates in boundary_arrays])
            all_dec = np.concatenate([coordinates[:, 1] for coordinates in boundary_arrays])
            transformed = self.transform_coords(all_ra, all_dec, frame_altaz)

        boundary_x, boundary_y = [], []
        offset = 0
        projected = self.projection_combo.currentIndex() in (1, 3, 4)
        for coordinates in boundary_arrays:
            point_count = len(coordinates)
            if is_local:
                azimuth = transformed[0][offset:offset + point_count]
                altitude = transformed[1][offset:offset + point_count]
                offset += point_count
                x_values, y_values = self.project_local(
                    self.local_azimuth(azimuth), altitude
                )
            else:
                altitude = None
                x_values, y_values = self.project_global(
                    coordinates[:, 0], coordinates[:, 1]
                )

            for index in range(point_count - 1):
                if projected:
                    append_projected_segment(
                        boundary_x, boundary_y,
                        x_values[index], y_values[index],
                        x_values[index + 1], y_values[index + 1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(
                        boundary_x, boundary_y,
                        x_values[index], y_values[index],
                        x_values[index + 1], y_values[index + 1],
                    )

        self.add_curve(
            boundary_x,
            boundary_y,
            pg.mkPen(
                color=(188, 137, 72, 105),
                width=0.75,
                style=Qt.PenStyle.DashLine,
            ),
            3,
        )

    def add_curve(self, x_values, y_values, pen, z_value=0):
        if len(x_values) == 0:
            return
        item = pg.PlotCurveItem(x=np.array(x_values), y=np.array(y_values), pen=pen, connect="finite")
        item.setZValue(z_value)
        self.plot_widget.addItem(item)

    def add_below_horizon_shading(self):
        """Shade the local hemisphere below altitude 0 without hiding objects."""
        x_range, y_range = self.plot_widget.getViewBox().viewRange()
        if self.projection_combo.currentIndex() == 3:
            # In the local stereographic projection altitude 0 is exactly the
            # circle r=2.  Keep that boundary as a vector path: rasterising the
            # mask made the horizon visibly stair-stepped at shallow slopes.
            outer_rect = QRectF(
                float(x_range[0]), float(y_range[0]),
                float(x_range[1] - x_range[0]),
                float(y_range[1] - y_range[0]),
            )
            ground_path = QPainterPath()
            ground_path.setFillRule(Qt.FillRule.OddEvenFill)
            ground_path.addRect(outer_rect)
            ground_path.addEllipse(QRectF(-2.0, -2.0, 4.0, 4.0))
            ground = AntialiasedPathItem(ground_path)
            ground.setBrush(QBrush(QColor(125, 82, 45, 38)))
            ground.setPen(QPen(Qt.PenStyle.NoPen))
            ground.setZValue(-20)
            self.plot_widget.addItem(ground)
            return

        columns = int(np.clip(self.plot_widget.width() // 4, 180, 420))
        rows = int(np.clip(self.plot_widget.height() // 4, 100, 260))
        x_values = np.linspace(x_range[0], x_range[1], columns)
        y_values = np.linspace(y_range[0], y_range[1], rows)
        grid_x, grid_y = np.meshgrid(x_values, y_values)
        _, altitude = self.inverse_local(grid_x, grid_y)
        below = np.isfinite(altitude) & (np.asarray(altitude) < 0.0)
        if self.projection_combo.currentIndex() == 1:
            below &= grid_x ** 2 / 8.0 + grid_y ** 2 / 2.0 <= 1.0
        image_data = np.zeros((rows, columns, 4), dtype=np.uint8)
        image_data[below] = (125, 82, 45, 38)
        ground = pg.ImageItem(image=image_data, axisOrder="row-major")
        ground.setRect(QRectF(
            float(x_range[0]), float(y_range[0]),
            float(x_range[1] - x_range[0]), float(y_range[1] - y_range[0]),
        ))
        ground.setZValue(-20)
        self.plot_widget.addItem(ground)

    def add_deep_sky_objects(self, frame_altaz, is_local):
        if (not self.dso_check.isChecked() or not self.deep_sky_objects
                or self.approximate_horizontal_fov() > 85.0):
            return
        ra = np.asarray([item["ra"] for item in self.deep_sky_objects])
        dec = np.asarray([item["dec"] for item in self.deep_sky_objects])
        if is_local:
            azimuth, altitude = self.transform_coords(ra, dec, frame_altaz)
            x_values, y_values = self.project_local(
                self.local_azimuth(azimuth), altitude
            )
        else:
            x_values, y_values = self.project_global(ra, dec)
        x_values = np.asarray(x_values)
        y_values = np.asarray(y_values)
        x_range, y_range = self.plot_widget.getViewBox().viewRange()
        visible = (
            np.isfinite(x_values) & np.isfinite(y_values)
            & (x_values >= x_range[0]) & (x_values <= x_range[1])
            & (y_values >= y_range[0]) & (y_values <= y_range[1])
        )

        ellipse = QPainterPath()
        ellipse.addEllipse(-0.75, -0.38, 1.5, 0.76)
        circle = QPainterPath()
        circle.addEllipse(-0.5, -0.5, 1.0, 1.0)
        groups = {
            "galaxy": (ellipse, (175, 145, 255, 150)),
            "cluster": (circle, (255, 205, 105, 155)),
            "nebula": (circle, (105, 225, 185, 155)),
        }
        for group_name, (symbol, color) in groups.items():
            indices = []
            sizes = []
            for index in np.where(visible)[0]:
                item = self.deep_sky_objects[index]
                object_type = item["type"]
                if object_type.startswith("G"):
                    category = "galaxy"
                elif object_type in {"OCl", "GCl", "*Ass", "Cl+N"}:
                    category = "cluster"
                else:
                    category = "nebula"
                if category != group_name:
                    continue
                indices.append(index)
                sizes.append(float(np.clip(5.0 + np.sqrt(item["major"]) * 1.5, 5.0, 18.0)))
            if not indices:
                continue
            marker = pg.ScatterPlotItem(
                x=x_values[indices], y=y_values[indices], size=sizes,
                symbol=symbol, brush=None,
                pen=pg.mkPen(color=color, width=1.0),
            )
            marker.setZValue(16)
            self.plot_widget.addItem(marker)

    def add_visible_comets(self, eval_time, frame_altaz, is_local):
        if not self.comets_check.isChecked():
            return
        comets = self.current_comet_positions(eval_time)
        if not comets:
            return
        ra = np.asarray([item[0] for item in comets])
        dec = np.asarray([item[1] for item in comets])
        azimuth, altitude = self.transform_coords(ra, dec, frame_altaz)
        above_horizon = np.asarray(altitude) >= 0.0
        if is_local:
            x_values, y_values = self.project_local(
                self.local_azimuth(azimuth), altitude
            )
        else:
            x_values, y_values = self.project_global(ra, dec)
        x_values = np.asarray(x_values)
        y_values = np.asarray(y_values)
        x_range, y_range = self.plot_widget.getViewBox().viewRange()
        visible = (
            above_horizon & np.isfinite(x_values) & np.isfinite(y_values)
            & (x_values >= x_range[0]) & (x_values <= x_range[1])
            & (y_values >= y_range[0]) & (y_values <= y_range[1])
        )
        indices = np.where(visible)[0]
        if len(indices) == 0:
            return
        magnitudes = np.asarray([item[2] for item in comets])
        marker = pg.ScatterPlotItem(
            x=x_values[indices], y=y_values[indices],
            size=np.clip(17.0 - 0.55 * magnitudes[indices], 6.0, 13.0),
            symbol="d", brush=pg.mkBrush(100, 235, 225, 210),
            pen=pg.mkPen((220, 255, 250, 240), width=1.0),
        )
        marker.setZValue(26)
        self.plot_widget.addItem(marker)

        if self.approximate_horizontal_fov() <= 85.0:
            for index in indices[np.argsort(magnitudes[indices])[:30]]:
                label = pg.TextItem(
                    text=f'☄ {comets[index][3]}  {magnitudes[index]:.1f}m',
                    color=(125, 245, 235, 225), anchor=(-0.08, 1.15),
                )
                label.setZValue(27)
                label.setPos(float(x_values[index]), float(y_values[index]))
                self.plot_widget.addItem(label)

    def add_planets(self, eval_time, frame_altaz, is_local):
        """Draw the seven major planets, omitting locally hidden bodies."""
        if not self.planets_check.isChecked():
            return

        label_font = QFont("Segoe UI", 9, QFont.Weight.DemiBold)
        for body_name, (label_text, color, marker_size) in PLANET_STYLES.items():
            body = get_body(body_name, eval_time)
            if is_local:
                body_altaz = body.transform_to(frame_altaz)
                altitude = float(body_altaz.alt.deg)
                if altitude < 0.0:
                    continue
                x_values, y_values = self.project_local(
                    self.local_azimuth([float(body_altaz.az.deg)]), [altitude]
                )
            else:
                x_values, y_values = self.project_global(
                    [float(body.ra.deg)], [float(body.dec.deg)]
                )

            x_value = float(np.asarray(x_values)[0])
            y_value = float(np.asarray(y_values)[0])
            if not np.isfinite(x_value) or not np.isfinite(y_value):
                continue

            halo = pg.ScatterPlotItem(
                x=[x_value], y=[y_value], size=marker_size * 2.4,
                brush=pg.mkBrush(color[0], color[1], color[2], 42), pen=None,
            )
            halo.setZValue(23)
            self.plot_widget.addItem(halo)

            marker = pg.ScatterPlotItem(
                x=[x_value], y=[y_value], size=marker_size,
                brush=pg.mkBrush(color),
                pen=pg.mkPen((245, 250, 255, 210), width=0.8),
            )
            marker.setZValue(24)
            self.plot_widget.addItem(marker)

            label = pg.TextItem(
                text=label_text, color=color, anchor=(-0.22, 0.68)
            )
            label.setFont(label_font)
            label.setZValue(25)
            label.setPos(x_value, y_value)
            self.plot_widget.addItem(label)

    def add_planner_overlays(self, frame_altaz, is_local):
        groups = [(self.planned_fields, (255, 70, 85, 235), 2.0, 31)]
        archive_colors = [
            (235, 95, 255, 220), (255, 165, 55, 220),
            (95, 235, 130, 220), (255, 225, 75, 220),
        ]
        for index, archive in enumerate(self.planner_store.data["archives"]):
            if archive.get("visible", True):
                groups.append((archive.get("fields", []), archive_colors[index % len(archive_colors)], 1.8, 29))

        for fields, color, width, z_value in groups:
            line_x, line_y = [], []
            center_x, center_y = [], []
            for field in fields:
                try:
                    vertices = field_outline(field)
                except (KeyError, TypeError, ValueError):
                    continue
                ra = np.asarray([point[0] for point in vertices])
                dec = np.asarray([point[1] for point in vertices])
                if is_local:
                    azimuth, altitude = self.transform_coords(ra, dec, frame_altaz)
                    px, py = self.project_local(self.local_azimuth(azimuth), altitude)
                    center_az, center_alt = self.transform_coords(
                        [field["ra_deg"]], [field["dec_deg"]], frame_altaz
                    )
                    cx, cy = self.project_local(self.local_azimuth(center_az), center_alt)
                else:
                    px, py = self.project_global(ra, dec)
                    cx, cy = self.project_global([field["ra_deg"]], [field["dec_deg"]])
                if np.isfinite(cx[0]) and np.isfinite(cy[0]):
                    center_x.append(float(cx[0]))
                    center_y.append(float(cy[0]))
                for edge_index in range(len(px) - 1):
                    values = (px[edge_index], py[edge_index], px[edge_index + 1], py[edge_index + 1])
                    if not all(np.isfinite(value) for value in values):
                        continue
                    if self.projection_combo.currentIndex() in (1, 3, 4):
                        append_projected_segment(
                            line_x, line_y, *values,
                            max_step=self.projected_segment_limit(0.8),
                        )
                    else:
                        append_sky_segment(line_x, line_y, *values)
            if line_x:
                outline = pg.PlotCurveItem(
                    x=np.asarray(line_x), y=np.asarray(line_y),
                    pen=pg.mkPen(color=color, width=width), connect="finite",
                )
                outline.setZValue(z_value)
                self.plot_widget.addItem(outline)
            if center_x:
                centres = pg.ScatterPlotItem(
                    x=center_x, y=center_y, size=5,
                    brush=pg.mkBrush(color), pen=pg.mkPen(color=color, width=0.8),
                )
                centres.setZValue(z_value + 1)
                self.plot_widget.addItem(centres)

    def refresh_bayer_labels(self):
        for item in self.bayer_label_items:
            try:
                self.plot_widget.removeItem(item)
            except RuntimeError:
                pass
        self.bayer_label_items = []

        if not self.bayer_labels_check.isChecked() or not self.bayer_label_data:
            return

        star_x, star_y, star_mag, star_bayer = self.bayer_label_data
        x_range, y_range = self.plot_widget.getViewBox().viewRange()
        x_min, x_max = sorted(x_range)
        y_min, y_max = sorted(y_range)
        x_span = x_max - x_min
        y_span = y_max - y_min
        projection = self.projection_combo.currentIndex()
        is_local = self.view_mode_combo.currentIndex() == 1
        if projection in (0, 2):
            full_width = 360.0
        elif projection == 1:
            full_width = 6.0
        else:
            full_width = 7.4 if is_local else 2.0 * STEREO_VIEW_RADIUS
        zoom_ratio = x_span / full_width
        if zoom_ratio > 0.80:
            magnitude_limit = 3.0
        elif zoom_ratio > 0.50:
            magnitude_limit = 3.5
        else:
            magnitude_limit = 4.0

        in_view = (
            (star_x >= x_min) & (star_x <= x_max) &
            (star_y >= y_min) & (star_y <= y_max)
        )
        candidates = np.where(
            in_view & (star_mag <= magnitude_limit) & (star_bayer != "")
        )[0]
        min_dx = x_span * 27.0 / max(self.plot_widget.width(), 1)
        min_dy = y_span * 17.0 / max(self.plot_widget.height(), 1)
        accepted = []
        label_font = QFont("Segoe UI Symbol", 9, QFont.Weight.DemiBold)
        for index in candidates[np.argsort(star_mag[candidates])]:
            x_value = float(star_x[index])
            y_value = float(star_y[index])
            if any(
                abs(x_value - other_x) < min_dx and abs(y_value - other_y) < min_dy
                for other_x, other_y in accepted
            ):
                continue
            accepted.append((x_value, y_value))
            label = pg.TextItem(
                text=str(star_bayer[index]),
                color=(218, 230, 255, 225),
                anchor=(-0.28, 0.62),
            )
            label.setFont(label_font)
            label.setZValue(22)
            label.setPos(x_value, y_value)
            self.plot_widget.addItem(label)
            self.bayer_label_items.append(label)

    def add_reference_grid(self, frame_altaz, is_local):
        projected = self.projection_combo.currentIndex() in (1, 3, 4)
        if is_local:
            major_x, major_y, minor_x, minor_y = [], [], [], []
            for azimuth in range(0, 360, 10):
                plot_azimuth = self.local_azimuth(azimuth)
                if projected:
                    alts = np.linspace(-30, 90, 61)
                    lx, ly = self.project_local(np.full_like(alts, plot_azimuth), alts)
                else:
                    lx, ly = self.project_local(
                        [plot_azimuth, plot_azimuth], [-30, 90]
                    )
                target_x, target_y = (
                    (major_x, major_y) if azimuth % 30 == 0
                    else (minor_x, minor_y)
                )
                target_x.extend(lx.tolist() + [np.nan])
                target_y.extend(ly.tolist() + [np.nan])

            self.add_curve(
                minor_x, minor_y,
                pg.mkPen(color=(100, 120, 140, 35), width=0.4),
            )
            self.add_curve(
                major_x, major_y,
                pg.mkPen(color=(120, 145, 165, 80), width=0.7),
            )

            az_samples = np.linspace(0, 360, 361)
            for altitude in range(-30, 90, 10):
                if altitude == 0:
                    color, width, z_value = (105, 62, 30, 255), 4.0, 20
                elif altitude < 0:
                    color, width, z_value = (165, 120, 80, 75), 0.65, 9
                elif altitude % 30 == 0:
                    color, width, z_value = (235, 165, 80, 160), 1.05, 11
                else:
                    color, width, z_value = (220, 150, 75, 135), 0.9, 10
                circle_az = self.local_azimuth(az_samples)
                cx, cy = self.project_local(circle_az, np.full_like(circle_az, altitude))

                res_x, res_y = [], []
                for i in range(len(cx) - 1):
                    if projected:
                        append_projected_segment(
                            res_x, res_y, cx[i], cy[i], cx[i+1], cy[i+1],
                            self.projected_segment_limit(),
                        )
                    else:
                        append_sky_segment(res_x, res_y, cx[i], cy[i], cx[i+1], cy[i+1])
                line_style = (
                    Qt.PenStyle.DotLine if altitude < 0
                    else Qt.PenStyle.SolidLine
                )
                self.add_curve(
                    res_x, res_y,
                    pg.mkPen(color=color, width=width, style=line_style),
                    z_value,
                )
            return

        major_x, major_y, minor_x, minor_y = [], [], [], []
        for ra in np.arange(0, 360, 7.5):
            samples = np.linspace(-90, 90, 91)
            x, y = self.project_global(np.full_like(samples, ra), samples)
            target_x, target_y = (major_x, major_y) if ra % 15 == 0 else (minor_x, minor_y)
            target_x.extend(x.tolist() + [np.nan])
            target_y.extend(y.tolist() + [np.nan])
        for dec in range(-90, 91, 5):
            samples = np.linspace(0, 360, 361)
            x, y = self.project_global(samples, np.full_like(samples, dec))
            target_x, target_y = (major_x, major_y) if dec % 20 == 0 else (minor_x, minor_y)
            for i in range(len(x)-1):
                if projected:
                    append_projected_segment(
                        target_x, target_y, x[i], y[i], x[i+1], y[i+1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(target_x, target_y, x[i], y[i], x[i+1], y[i+1])
        self.add_curve(minor_x, minor_y, pg.mkPen(color=(100, 120, 140, 45), width=0.45))
        self.add_curve(major_x, major_y, pg.mkPen(color=(150, 170, 190, 100), width=0.9))

    def add_reference_circles(self, frame_altaz, is_local):
        equator_pen = pg.mkPen(color=(80, 220, 255, 210), width=1.3)
        ecliptic_pen = pg.mkPen(
            color=(255, 155, 35, 240),
            width=2.0,
            style=Qt.PenStyle.DashLine,
        )

        if is_local:
            ra = np.linspace(0, 360, 361)
            dec = np.zeros_like(ra)
            az, alt = self.transform_coords(ra, dec, frame_altaz)
            az = self.local_azimuth(az)
            equator_x, equator_y = [], []
            for i in range(len(az) - 1):
                line_x, line_y = self.project_local([az[i], az[i + 1]], [alt[i], alt[i + 1]])
                if self.projection_combo.currentIndex() in (1, 3, 4):
                    append_projected_segment(
                        equator_x, equator_y,
                        line_x[0], line_y[0], line_x[1], line_y[1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(equator_x, equator_y, line_x[0], line_y[0], line_x[1], line_y[1])
            self.add_curve(equator_x, equator_y, equator_pen, 4)
        else:
            equator_ra = np.linspace(0, 360, 361)
            equator_x, equator_y = self.project_global(equator_ra, np.zeros_like(equator_ra))
            self.add_curve(
                equator_x.tolist() + [np.nan],
                equator_y.tolist() + [np.nan],
                equator_pen,
                4,
            )

        longitude = np.linspace(0, 360, 361)
        obliquity = np.radians(23.4393)
        ecliptic_ra = np.degrees(np.arctan2(
            np.sin(np.radians(longitude)) * np.cos(obliquity),
            np.cos(np.radians(longitude))
        )) % 360
        ecliptic_dec = np.degrees(np.arcsin(
            np.sin(obliquity) * np.sin(np.radians(longitude))
        ))

        if is_local:
            az, alt = self.transform_coords(ecliptic_ra, ecliptic_dec, frame_altaz)
            az = self.local_azimuth(az)
            ecliptic_x, ecliptic_y = [], []
            for i in range(len(az) - 1):
                line_x, line_y = self.project_local([az[i], az[i + 1]], [alt[i], alt[i + 1]])
                if self.projection_combo.currentIndex() in (1, 3, 4):
                    append_projected_segment(
                        ecliptic_x, ecliptic_y,
                        line_x[0], line_y[0], line_x[1], line_y[1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(ecliptic_x, ecliptic_y, line_x[0], line_y[0], line_x[1], line_y[1])
        else:
            projected_x, projected_y = self.project_global(ecliptic_ra, ecliptic_dec)
            ecliptic_x, ecliptic_y = [], []
            for index in range(len(projected_x) - 1):
                if self.projection_combo.currentIndex() in (1, 3, 4):
                    append_projected_segment(
                        ecliptic_x, ecliptic_y,
                        projected_x[index], projected_y[index],
                        projected_x[index + 1], projected_y[index + 1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(
                        ecliptic_x, ecliptic_y,
                        projected_x[index], projected_y[index],
                        projected_x[index + 1], projected_y[index + 1],
                    )

        self.add_curve(ecliptic_x, ecliptic_y, ecliptic_pen, 4)

        galactic = SkyCoord(l=np.linspace(0, 360, 361) * u.deg, b=np.zeros(361) * u.deg, frame=Galactic).icrs
        if is_local:
            gal_x, gal_y = self.transform_coords(galactic.ra.deg, galactic.dec.deg, frame_altaz)
            gal_x = self.local_azimuth(gal_x)
            galactic_x, galactic_y = [], []
            for i in range(len(gal_x) - 1):
                line_x, line_y = self.project_local(
                    [gal_x[i], gal_x[i + 1]], [gal_y[i], gal_y[i + 1]])
                if self.projection_combo.currentIndex() in (1, 3, 4):
                    append_projected_segment(
                        galactic_x, galactic_y,
                        line_x[0], line_y[0], line_x[1], line_y[1],
                        self.projected_segment_limit(),
                    )
                else:
                    append_sky_segment(
                        galactic_x, galactic_y,
                        line_x[0], line_y[0], line_x[1], line_y[1],
                    )
        else:
            galactic_x, galactic_y = self.project_global(galactic.ra.deg, galactic.dec.deg)
            if self.projection_combo.currentIndex() in (1, 3, 4):
                split_x, split_y = [], []
                for i in range(len(galactic_x) - 1):
                    if (abs(galactic_x[i + 1] - galactic_x[i]) < 0.3 and
                            abs(galactic_y[i + 1] - galactic_y[i]) < 0.3):
                        split_x.extend([galactic_x[i], galactic_x[i + 1], np.nan])
                        split_y.extend([galactic_y[i], galactic_y[i + 1], np.nan])
                galactic_x, galactic_y = split_x, split_y
            else:
                galactic_x, galactic_y = galactic_x.tolist(), galactic_y.tolist()
            if self.projection_combo.currentIndex() not in (1, 3, 4):
                wrapped_x, wrapped_y = [], []
                for i in range(len(galactic_x) - 1):
                    append_sky_segment(wrapped_x, wrapped_y, galactic_x[i], galactic_y[i], galactic_x[i + 1], galactic_y[i + 1])
                galactic_x, galactic_y = wrapped_x, wrapped_y
            else:
                galactic_x.append(np.nan)
                galactic_y.append(np.nan)
        self.add_curve(galactic_x, galactic_y, pg.mkPen(color=(220, 80, 220, 180), width=1.2, style=Qt.PenStyle.DashLine), 4)

    def redraw(self):
        self.plot_widget.clear()
        self.preview_items = []
        self.preview_fields = []
        self.bayer_label_items = []
        self.bayer_label_data = []
        
        qdt = self.time_edit.dateTime().toPyDateTime()
        eval_time = Time(qdt)
        obs_location = OBSERVATORIES[self.obs_combo.currentText()]
        frame_altaz = AltAz(obstime=eval_time, location=obs_location)
        is_local = self.view_mode_combo.currentIndex() == 1

        # 1. Slunce a Měsíc
        sun = get_sun(eval_time)
        moon = get_body("moon", eval_time)
        self.update_sun_altitude_readout(
            float(sun.transform_to(frame_altaz).alt.deg)
        )

        if is_local:
            sun_x, sun_y = self.transform_coords(sun.ra.deg, sun.dec.deg, frame_altaz)
            moon_x, moon_y = self.transform_coords(moon.ra.deg, moon.dec.deg, frame_altaz)
            sun_altitude = float(sun_y[0])
            moon_altitude = float(moon_y[0])
            self.local_sun_azimuth = float(sun_x[0])
            if self.projection_combo.currentIndex() == 4:
                observatory = self.obs_combo.currentText()
                if self.local_mercator_observatory != observatory:
                    facing_azimuth, center_altitude = (
                        self.preferred_local_mercator_center(obs_location)
                    )
                    self.local_mercator_center_azimuth = facing_azimuth
                    self.local_mercator_center_altitude = center_altitude
                    self.local_mercator_observatory = observatory
                self.mercator_center_azimuth = float(
                    self.local_azimuth(self.local_mercator_center_azimuth)
                )
                self.mercator_center_altitude = self.local_mercator_center_altitude
            sun_x = self.local_azimuth(sun_x)
            moon_x = self.local_azimuth(moon_x)
            sun_x, sun_y = self.project_local(sun_x, sun_y)
            moon_x, moon_y = self.project_local(moon_x, moon_y)
            lower_altitude = max(-30.0, sun_altitude - 10.0)
            local_lower_altitude = lower_altitude
            if self.projection_combo.currentIndex() == 1:
                self.plot_widget.setXRange(-3, 3, padding=0)
                self.plot_widget.setYRange(-1.5, 1.5, padding=0)
            elif self.projection_combo.currentIndex() == 3:
                radius = 2.0 * np.tan(np.radians((90.0 - lower_altitude) / 2.0))
                radius *= 1.04
                self.plot_widget.setXRange(-radius, radius, padding=0)
                self.plot_widget.setYRange(-radius, radius, padding=0)
            elif self.projection_combo.currentIndex() == 4:
                # MER si zachovává rozsah při běžném plynulém posouvání.
                pass
            else:
                _, projected_limits = self.project_local(
                    [180.0, 180.0], [lower_altitude, 90.0]
                )
                self.plot_widget.setYRange(
                    float(projected_limits[0]),
                    float(projected_limits[1]),
                    padding=0,
                )
        else:
            self.local_sun_azimuth = None
            sun_altitude = float(sun.dec.deg)
            moon_altitude = float(moon.dec.deg)
            sun_x, sun_y = self.project_global([sun.ra.deg], [sun.dec.deg])
            moon_x, moon_y = self.project_global([moon.ra.deg], [moon.dec.deg])

        if is_local:
            self.add_below_horizon_shading()
        self.add_reference_grid(frame_altaz, is_local)
        self.add_reference_circles(frame_altaz, is_local)
        self.add_constellation_boundaries(frame_altaz, is_local)
        self.add_deep_sky_objects(frame_altaz, is_local)

        if self.stars:
            star_ra = np.array([star[0] for star in self.stars])
            star_dec = np.array([star[1] for star in self.stars])
            star_mag = np.array([star[2] for star in self.stars])
            star_bayer = np.array([
                star[3] if len(star) > 3 else "" for star in self.stars
            ], dtype=object)
            star_colors = np.array([
                star[5] if len(star) > 5 else (225, 232, 255, 235)
                for star in self.stars
            ], dtype=object)
            catalog_visible = star_mag <= self.visible_star_magnitude()
            if is_local:
                star_az, star_alt = self.transform_coords(
                    star_ra[catalog_visible], star_dec[catalog_visible], frame_altaz
                )
                star_x, star_y = self.project_local(
                    self.local_azimuth(star_az), star_alt
                )
            else:
                star_x, star_y = self.project_global(
                    star_ra[catalog_visible], star_dec[catalog_visible]
                )
            visible_mag = star_mag[catalog_visible]
            visible_bayer = star_bayer[catalog_visible]
            visible_colors = star_colors[catalog_visible]
            finite = np.isfinite(star_x) & np.isfinite(star_y)
            star_x = np.asarray(star_x)[finite]
            star_y = np.asarray(star_y)[finite]
            visible_mag = visible_mag[finite]
            visible_bayer = visible_bayer[finite]
            visible_colors = visible_colors[finite]

            star_size = np.clip(9.2 - 1.55 * visible_mag, 1.6, 9.5)
            bright = visible_mag <= 2.5
            if np.any(bright):
                halo_item = pg.ScatterPlotItem(
                    x=star_x[bright],
                    y=star_y[bright],
                    size=star_size[bright] * 2.7,
                    brush=(135, 175, 255, 35),
                    pen=None,
                )
                halo_item.setZValue(17)
                self.plot_widget.addItem(halo_item)

            star_item = pg.ScatterPlotItem(
                x=star_x,
                y=star_y,
                size=star_size,
                brush=[pg.mkBrush(tuple(color)) for color in visible_colors],
                pen=pg.mkPen((245, 250, 255, 120), width=0.35),
            )
            star_item.setZValue(18)
            self.plot_widget.addItem(star_item)

            self.bayer_label_data = (
                star_x.copy(),
                star_y.copy(),
                visible_mag.copy(),
                visible_bayer.copy(),
            )
            self.refresh_bayer_labels()

        if np.all(np.isfinite(sun_x)) and np.all(np.isfinite(sun_y)):
            self.plot_widget.addItem(pg.ScatterPlotItem(x=sun_x, y=sun_y, size=18, brush='#ffcc00', pen=pg.mkPen('w', width=1.5)))
        if np.all(np.isfinite(moon_x)) and np.all(np.isfinite(moon_y)):
            self.plot_widget.addItem(pg.ScatterPlotItem(x=moon_x, y=moon_y, size=12, brush='#eeeeee', pen=pg.mkPen('b', width=1)))
        self.add_planets(eval_time, frame_altaz, is_local)
        self.add_visible_comets(eval_time, frame_altaz, is_local)

        # 2. Elongační kružnice kolem Slunce (20° až 70°)
        angles = np.linspace(0, 2*np.pi, 120)
        elongation_styles = {
            20: ((255, 235, 120, 190), Qt.PenStyle.SolidLine),
            30: ((255, 190, 80, 190), Qt.PenStyle.DashLine),
            40: ((100, 220, 255, 190), Qt.PenStyle.DotLine),
            50: ((150, 255, 150, 190), Qt.PenStyle.DashDotLine),
            60: ((255, 110, 180, 190), Qt.PenStyle.DashDotDotLine),
            70: ((180, 140, 255, 200), Qt.PenStyle.SolidLine),
        }
        for elongation in [20, 30, 40, 50, 60, 70]:
            e_rad = np.radians(elongation)
            s_ra = np.radians(sun.ra.deg)
            s_dec = np.radians(sun.dec.deg)

            dec_c = np.arcsin(np.sin(s_dec) * np.cos(e_rad) + np.cos(s_dec) * np.sin(e_rad) * np.cos(angles))
            # Ošetření přetečení deklinace mimo rozsah [-90, +90]
            dec_deg_c = np.clip(np.degrees(dec_c), -89.9, 89.9)
            
            ra_c = s_ra + np.arctan2(np.sin(angles) * np.sin(e_rad) * np.cos(s_dec), np.cos(e_rad) - np.sin(s_dec) * np.sin(dec_c))
            ra_deg_c = np.degrees(ra_c) % 360

            if is_local:
                x_c, y_c = self.transform_coords(ra_deg_c, dec_deg_c, frame_altaz)
                x_c = self.local_azimuth(x_c)
                projected_x, projected_y = [], []
                for i in range(len(x_c) - 1):
                    line_x, line_y = self.project_local(
                        [x_c[i], x_c[i + 1]], [y_c[i], y_c[i + 1]])
                    if self.projection_combo.currentIndex() in (1, 3, 4):
                        append_projected_segment(
                            projected_x, projected_y,
                            line_x[0], line_y[0], line_x[1], line_y[1],
                            self.projected_segment_limit(),
                        )
                    else:
                        append_sky_segment(
                            projected_x, projected_y,
                            line_x[0], line_y[0], line_x[1], line_y[1])
                x_c, y_c = np.asarray(projected_x), np.asarray(projected_y)
            else:
                raw_x, raw_y = self.project_global(ra_deg_c, dec_deg_c)
                projected_x, projected_y = [], []
                for index in range(len(raw_x) - 1):
                    if self.projection_combo.currentIndex() in (1, 3, 4):
                        append_projected_segment(
                            projected_x, projected_y,
                            raw_x[index], raw_y[index],
                            raw_x[index + 1], raw_y[index + 1],
                            self.projected_segment_limit(),
                        )
                    else:
                        append_sky_segment(
                            projected_x, projected_y,
                            raw_x[index], raw_y[index],
                            raw_x[index + 1], raw_y[index + 1],
                        )
                x_c, y_c = np.asarray(projected_x), np.asarray(projected_y)

            if len(x_c) > 0:
                color, style = elongation_styles[elongation]
                self.plot_widget.addItem(pg.PlotCurveItem(x=x_c, y=y_c, pen=pg.mkPen(color=color, width=1.4, style=style)))

        self.add_kreutz_corridors(eval_time, frame_altaz, is_local)

        # 3. Zorná pole z archivu
        if self.polygons:
            limit_mag = self.mag_spin.value()
            x_all, y_all = [], []
            visible_polygons = [polygon for polygon in self.polygons if polygon[2] >= limit_mag]

            if self.projection_combo.currentIndex() == 4 and visible_polygons:
                cache_key = (
                    id(self.polygons), len(visible_polygons), float(limit_mag), is_local,
                    self.obs_combo.currentText(), eval_time.utc.isot,
                    round(self.mercator_center_ra, 9), round(self.mercator_center_dec, 9),
                    round(self.mercator_center_azimuth, 9),
                    round(self.mercator_center_altitude, 9),
                )
                if self.mer_polygon_cache is not None and self.mer_polygon_cache[0] == cache_key:
                    px, py = self.mer_polygon_cache[1], self.mer_polygon_cache[2]
                else:
                    ra_matrix = np.asarray([polygon[0] for polygon in visible_polygons], dtype=float)
                    dec_matrix = np.asarray([polygon[1] for polygon in visible_polygons], dtype=float)
                    if is_local:
                        shape = ra_matrix.shape
                        azimuth, altitude = self.transform_coords(
                            ra_matrix.ravel(), dec_matrix.ravel(), frame_altaz
                        )
                        altitude = np.asarray(altitude).reshape(shape)
                        azimuth = np.asarray(azimuth).reshape(shape)
                        px, py = self.project_local(self.local_azimuth(azimuth), altitude)
                    else:
                        px, py = self.project_global(ra_matrix, dec_matrix)
                    self.mer_polygon_cache = (cache_key, px, py)
                shape = px.shape

                x_range, y_range = self.plot_widget.getViewBox().viewRange()
                x_margin = abs(x_range[1] - x_range[0]) * 0.12
                y_margin = abs(y_range[1] - y_range[0]) * 0.12
                outside = (
                    (px < x_range[0] - x_margin) | (px > x_range[1] + x_margin)
                    | (py < y_range[0] - y_margin) | (py > y_range[1] + y_margin)
                )
                px = np.where(outside, np.nan, px)
                py = np.where(outside, np.nan, py)
                visible_rows = np.any(np.isfinite(px) & np.isfinite(py), axis=1)
                px = px[visible_rows]
                py = py[visible_rows]
                edge_x = np.stack((px[:, :-1], px[:, 1:]), axis=2).reshape(-1, 2)
                edge_y = np.stack((py[:, :-1], py[:, 1:]), axis=2).reshape(-1, 2)
                valid_edges = np.all(np.isfinite(edge_x) & np.isfinite(edge_y), axis=1)
                x_all = edge_x[valid_edges].ravel()
                y_all = edge_y[valid_edges].ravel()

            transformed_vertices = None
            if self.projection_combo.currentIndex() != 4 and is_local and visible_polygons:
                all_ra = np.concatenate([np.asarray(polygon[0]) for polygon in visible_polygons])
                all_dec = np.concatenate([np.asarray(polygon[1]) for polygon in visible_polygons])
                transformed_vertices = self.transform_coords(all_ra, all_dec, frame_altaz)
                vertex_offset = 0

            for ra, dec, mag in ([] if self.projection_combo.currentIndex() == 4 else visible_polygons):
                ra_arr = np.array(ra)
                dec_arr = np.array(dec)

                if is_local:
                    vertex_count = len(ra_arr)
                    az = transformed_vertices[0][vertex_offset:vertex_offset + vertex_count]
                    alt = transformed_vertices[1][vertex_offset:vertex_offset + vertex_count]
                    vertex_offset += vertex_count
                    px, py = self.project_local(self.local_azimuth(az), alt)
                else:
                    px, py = self.project_global(ra_arr, dec_arr)

                for i in range(len(px) - 1):
                    if self.projection_combo.currentIndex() in (1, 3, 4):
                        if (abs(px[i + 1] - px[i]) < 0.35 and
                                abs(py[i + 1] - py[i]) < 0.35):
                            x_all.extend([px[i], px[i + 1], np.nan])
                            y_all.extend([py[i], py[i + 1], np.nan])
                    else:
                        append_sky_segment(x_all, y_all, px[i], py[i], px[i+1], py[i+1])

            if len(x_all) > 0:
                if self.projection_combo.currentIndex() == 4:
                    poly_item = pg.PlotCurveItem(
                        x=np.asarray(x_all), y=np.asarray(y_all),
                        pen=pg.mkPen(color=(0, 180, 220, 40), width=0.5),
                        connect="pairs", skipFiniteCheck=True,
                    )
                    poly_item.setSegmentedLineMode("on")
                    poly_item.setCacheMode(
                        QGraphicsItem.CacheMode.DeviceCoordinateCache,
                        QSize(4096, 4096),
                    )
                else:
                    poly_item = pg.PlotCurveItem(
                        x=np.array(x_all), y=np.array(y_all),
                        pen=pg.mkPen(color=(0, 180, 220, 40), width=0.5),
                        connect="finite",
                    )
                poly_item.setZValue(1)
                self.plot_widget.addItem(poly_item)

        # 4. Souhvězdí
        if self.const_lines:
            c_x, c_y = [], []
            segments = []
            for line in self.const_lines:
                for i in range(len(line) - 1):
                    p1, p2 = line[i], line[i+1]
                    segments.append((p1[0] % 360, p1[1], p2[0] % 360, p2[1]))

            if segments:
                segment_array = np.asarray(segments, dtype=float)
                ra_pairs = segment_array[:, (0, 2)]
                dec_pairs = segment_array[:, (1, 3)]
                if is_local:
                    azimuth, altitude = self.transform_coords(
                        ra_pairs.ravel(), dec_pairs.ravel(), frame_altaz
                    )
                    line_x, line_y = self.project_local(
                        self.local_azimuth(np.asarray(azimuth).reshape(-1, 2)),
                        np.asarray(altitude).reshape(-1, 2),
                    )
                else:
                    line_x, line_y = self.project_global(ra_pairs, dec_pairs)

                projected = self.projection_combo.currentIndex() in (1, 3, 4)
                for index in range(len(segments)):
                    if projected:
                        append_projected_segment(
                            c_x, c_y,
                            line_x[index, 0], line_y[index, 0],
                            line_x[index, 1], line_y[index, 1],
                            self.projected_segment_limit(),
                        )
                    else:
                        append_sky_segment(
                            c_x, c_y,
                            line_x[index, 0], line_y[index, 0],
                            line_x[index, 1], line_y[index, 1],
                        )

            if c_x:
                const_item = pg.PlotCurveItem(
                    x=np.array(c_x),
                    y=np.array(c_y),
                    pen=pg.mkPen(color=(105, 220, 145, 155), width=0.9),
                    connect="finite",
                )
                const_item.setZValue(15)
                self.plot_widget.addItem(const_item)

        # Aktivní plán je červený, zapnuté archivní plány mají tlumené barvy.
        self.add_planner_overlays(frame_altaz, is_local)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
