import sys
import os
import json
import csv
import gzip
import urllib.request
from datetime import datetime, timedelta, timezone
import numpy as np

from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QSpinBox, QDoubleSpinBox, 
                             QPushButton, QComboBox, QDateTimeEdit, QCheckBox,
                             QFrame)
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QDateTime, QTimer
from PyQt6.QtGui import QFont
import pyqtgraph as pg

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (get_sun, get_body, EarthLocation, SkyCoord,
                                 AltAz, Galactic, GeocentricTrueEcliptic,
                                 HeliocentricMeanEcliptic, get_body_barycentric)


OBSERVATORIES = {
    '048 - Kleť (CZ)': EarthLocation.from_geodetic(14.2858, 48.8639, 1070),
    '950 - La Palma (ES)': EarthLocation.from_geodetic(-17.8816, 28.7603, 2396),
    'I05 - Paranal (CL)': EarthLocation.from_geodetic(-70.4042, -24.6272, 2635),
    'I47 - Pierre Auger (AR)': EarthLocation.from_geodetic(-69.3158, -35.2060, 1400)
}

SKYCOV_RETENTION_DAYS = 62
LAMBERT_Y_SCALE = 180.0 / np.pi
STEREO_VIEW_RADIUS = 7.5

BAYER_GREEK = {
    "Alp": "α", "Bet": "β", "Gam": "γ", "Del": "δ", "Eps": "ε",
    "Zet": "ζ", "Eta": "η", "The": "θ", "Iot": "ι", "Kap": "κ",
    "Lam": "λ", "Mu": "μ", "Nu": "ν", "Xi": "ξ", "Omi": "ο",
    "Pi": "π", "Rho": "ρ", "Sig": "σ", "Tau": "τ", "Ups": "υ",
    "Phi": "φ", "Chi": "χ", "Psi": "ψ", "Ome": "ω",
}
SUPERSCRIPT_DIGITS = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


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

    def __init__(self, dir_path, days_back):
        super().__init__()
        self.dir_path = dir_path
        self.days_back = days_back

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
        min_date = max_date - timedelta(days=self.days_back)
        target_files = [p for p in all_files if min_date <= p[0] <= max_date]

        removed_mb = removed_bytes / (1024 * 1024)
        cleanup_status = f" Odstraněno: {removed_count} souborů ({removed_mb:.1f} MiB)."
        cleanup_status += f" Prázdné složky: {removed_dirs}."
        if removal_errors:
            cleanup_status += f" Nešlo odstranit souborů: {removal_errors}."
        if directory_errors:
            cleanup_status += f" Nešlo odstranit složek: {directory_errors}."
        self.progress.emit(f"Načítám {len(target_files)} souborů.{cleanup_status}")

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
        self.progress.emit(
            f"Připraveno. Nejnovější MPC data: {max_date:%Y-%m-%d}."
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
                        if magnitude <= 5.0:
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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SkyCoverage — MPC sky map")
        self.resize(1600, 920)
        self.setMinimumSize(1180, 720)
        self.apply_modern_style()

        self.selected_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "skycov")
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
        subtitle_label = QLabel("MPC coverage · planetárium · plánování pozorování")
        subtitle_label.setObjectName("appSubtitle")
        brand_layout.addWidget(title_label)
        brand_layout.addWidget(subtitle_label)
        header_layout.addLayout(brand_layout)
        header_layout.addStretch()
        main_layout.addWidget(header_frame)

        # Panel 1: Pohled, Observatoř, Čas, Mag
        primary_controls = QFrame()
        primary_controls.setObjectName("controlCard")
        control_panel1 = QHBoxLayout(primary_controls)
        control_panel1.setContentsMargins(12, 9, 12, 9)
        control_panel1.setSpacing(8)

        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItems(["Globální (RA/Dec)", "Lokální obzor (Alt/Az)"])
        self.view_mode_combo.setObjectName("accentControl")
        self.view_mode_combo.currentIndexChanged.connect(self.on_view_change)
        control_panel1.addWidget(self.view_mode_combo)

        control_panel1.addWidget(QLabel("Projekce:"))
        self.projection_combo = QComboBox()
        self.projection_combo.addItems([
            "Obdélníková",
            "Hammer-Aitoff",
            "Válcová plochojevná",
            "Stereografická (věrné tvary)",
        ])
        self.projection_combo.setToolTip(
            "Stereografická projekce zachovává místní tvary a úhly souhvězdí"
        )
        self.projection_combo.currentIndexChanged.connect(self.on_projection_change)
        control_panel1.addWidget(self.projection_combo)

        control_panel1.addWidget(QLabel("Observatoř:"))
        self.obs_combo = QComboBox()
        self.obs_combo.addItems(list(OBSERVATORIES.keys()))
        self.obs_combo.currentIndexChanged.connect(self.request_redraw)
        control_panel1.addWidget(self.obs_combo)

        control_panel1.addWidget(QLabel("Čas (UTC):"))
        self.time_edit = QDateTimeEdit(QDateTime.currentDateTimeUtc())
        self.time_edit.setTimeSpec(Qt.TimeSpec.UTC)
        self.time_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.time_edit.dateTimeChanged.connect(self.request_redraw)
        control_panel1.addWidget(self.time_edit)

        control_panel1.addWidget(QLabel("Dny zpět:"))
        self.days_spin = QSpinBox()
        self.days_spin.setRange(1, SKYCOV_RETENTION_DAYS)
        self.days_spin.setValue(30)
        control_panel1.addWidget(self.days_spin)

        control_panel1.addWidget(QLabel("Min. mag:"))
        self.mag_spin = QDoubleSpinBox()
        self.mag_spin.setRange(5.0, 25.0)
        self.mag_spin.setSingleStep(0.5)
        self.mag_spin.setValue(15.0)
        control_panel1.addWidget(self.mag_spin)

        self.btn_load = QPushButton("Načíst data")
        self.btn_load.setObjectName("primaryButton")
        self.btn_load.clicked.connect(self.start_processing)
        control_panel1.addWidget(self.btn_load)

        self.kreutz_check = QCheckBox("Kreutzův koridor I/II")
        self.kreutz_check.setToolTip(
            "Model koridorů Kreutz I a II pro 5, 10, 15, 20, 30, 40 a 45 dní před přísluním; "
            "značky ukazují 5 dní a potom intervaly po 10 dnech"
        )
        self.kreutz_check.setObjectName("warmCheck")
        self.kreutz_check.stateChanged.connect(self.request_redraw)
        control_panel1.addWidget(self.kreutz_check)

        self.coord_label = QLabel(
            "Az: -- | Alt: --\n"
            "RA: -- | Dec: --\n"
            "Elongace: --"
        )
        self.coord_label.setTextFormat(Qt.TextFormat.PlainText)
        self.coord_label.setWordWrap(False)
        self.coord_label.setFixedSize(360, 64)
        self.coord_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self.coord_label.setObjectName("coordReadout")
        control_panel1.addStretch()
        header_layout.addWidget(self.coord_label)

        main_layout.addWidget(primary_controls)

        # Panel 2: Tlačítka soumraků
        secondary_controls = QFrame()
        secondary_controls.setObjectName("controlCard")
        control_panel2 = QHBoxLayout(secondary_controls)
        control_panel2.setContentsMargins(12, 8, 12, 8)
        control_panel2.setSpacing(7)
        control_panel2.addWidget(QLabel("Soumrak:"))

        self.btn_naut_eve = QPushButton("Večerní Nautický (-12°)")
        self.btn_naut_eve.clicked.connect(lambda: self.set_twilight(-12.0, evening=True))
        control_panel2.addWidget(self.btn_naut_eve)

        self.btn_astro_eve = QPushButton("Večerní Astro Noc (-18°)")
        self.btn_astro_eve.clicked.connect(lambda: self.set_twilight(-18.0, evening=True))
        control_panel2.addWidget(self.btn_astro_eve)

        self.btn_astro_morn = QPushButton("Ranní Astro Noc (-18°)")
        self.btn_astro_morn.clicked.connect(lambda: self.set_twilight(-18.0, evening=False))
        control_panel2.addWidget(self.btn_astro_morn)

        self.btn_naut_morn = QPushButton("Ranní Nautický (-12°)")
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
        control_panel2.addWidget(self.const_boundaries_check)

        self.bayer_labels_check = QCheckBox("Bayer αβγ")
        self.bayer_labels_check.setChecked(True)
        self.bayer_labels_check.setToolTip(
            "Zobrazit řecká Bayerova označení jasných hvězd; hustota se přizpůsobuje přiblížení"
        )
        self.bayer_labels_check.setObjectName("coolCheck")
        self.bayer_labels_check.stateChanged.connect(self.request_redraw)
        control_panel2.addWidget(self.bayer_labels_check)

        control_panel2.addStretch()
        main_layout.addWidget(secondary_controls)

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
        main_layout.addWidget(self.plot_widget)

        self.status_label = QLabel("Připraveno.")
        self.status_label.setObjectName("statusLabel")
        main_layout.addWidget(self.status_label)

        self.polygons = []
        self.stars = []
        self.ref_time = None
        self.local_sun_azimuth = None
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
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_axis_tick_update
        )
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_bayer_label_update
        )

        self.days_spin.valueChanged.connect(self.request_data_reload)
        self.mag_spin.valueChanged.connect(self.request_redraw)
        self.obs_combo.currentIndexChanged.connect(self.request_redraw)
        self.time_edit.dateTimeChanged.connect(self.request_redraw)
        star_cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "astro_cache", "hyg_v38.csv.gz")
        self.star_reader = StarCatalogReader(star_cache)
        self.star_reader.finished.connect(self.on_stars_loaded)
        self.star_reader.progress.connect(self.status_label.setText)
        self.star_reader.start()
        self.update_axis_ticks()
        self.start_processing()

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
            QComboBox, QSpinBox, QDoubleSpinBox, QDateTimeEdit {
                color: #e8f0fb;
                background-color: #152137;
                border: 1px solid #2a3a56;
                border-radius: 7px;
                padding: 6px 9px;
                min-height: 20px;
                selection-background-color: #245e72;
            }
            QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover, QDateTimeEdit:hover {
                border-color: #3e6380;
                background-color: #182740;
            }
            QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus, QDateTimeEdit:focus {
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
        self.redraw()

    def shift_time(self, minutes):
        self.time_edit.setDateTime(
            self.time_edit.dateTime().addSecs(int(minutes * 60))
        )

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
            else:
                self.plot_widget.setXRange(0, 360)
                upper_y = LAMBERT_Y_SCALE if projection == 2 else 90.0
                self.plot_widget.setYRange(-upper_y, upper_y)
                self.plot_widget.setLabel('bottom', 'Rektascenze RA (hodiny)')
                self.plot_widget.setLabel('left', 'Deklinace Dec (stupně)')
        self.redraw()
        self.update_axis_ticks()

    def on_projection_change(self):
        self.on_view_change()

    def update_axis_ticks(self):
        is_local = self.view_mode_combo.currentIndex() == 1
        projection = self.projection_combo.currentIndex()
        hammer = projection == 1
        lambert = projection == 2
        stereo = projection == 3
        bottom = self.plot_widget.getAxis('bottom')
        left = self.plot_widget.getAxis('left')
        if stereo:
            if is_local:
                bottom.setLabel('Stereografická obloha · zenit uprostřed')
            else:
                bottom.setLabel('Stereografická · střed RA 12h, Dec 0°')
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
        target_altaz = target.transform_to(frame_altaz)
        elongation = target_altaz.separation(sun_altaz).deg
        target = target.transform_to('icrs')
        ra_hours = target.ra.deg / 15.0
        display_azimuth = float(target_altaz.az.deg)
        display_altitude = float(target_altaz.alt.deg)
        self.coord_label.setText(
            f"Az: {display_azimuth:5.1f}° | Alt: {display_altitude:+5.1f}°\n"
            f"RA: {int(ra_hours):02d}h {int((ra_hours % 1) * 60):02d}m "
            f"| Dec: {target.dec.deg:+5.1f}°\n"
            f"Elongace: {elongation:.1f}°"
        )

    def show_invalid_coordinates(self):
        self.coord_label.setText(
            "Az: -- | Alt: --\n"
            "RA: -- | Dec: --\n"
            "Mimo platnou oblast projekce"
        )

    def request_redraw(self):
        self.redraw_timer.start(80)

    def request_axis_tick_update(self, *_):
        self.axis_tick_timer.start(50)

    def request_bayer_label_update(self, *_):
        self.bayer_label_timer.start(120)

    def request_data_reload(self):
        self.reload_timer.start(250)

    def start_processing(self):
        self.btn_load.setEnabled(False)
        self.reader = PointingDataReader(self.selected_dir, self.days_spin.value())
        self.reader.progress.connect(self.status_label.setText)
        self.reader.finished.connect(self.on_data_loaded)
        self.reader.start()

    def on_data_loaded(self, result):
        self.polygons, self.ref_time = result
        self.uncovered_cache.clear()
        self.btn_load.setEnabled(True)
        self.redraw()

    def on_stars_loaded(self, result):
        self.stars, error = result
        if error:
            self.status_label.setText(error)
        self.request_redraw()

    def transform_coords(self, ra_arr, dec_arr, frame_altaz):
        sc = SkyCoord(ra=np.atleast_1d(ra_arr)*u.deg, dec=np.atleast_1d(dec_arr)*u.deg)
        altaz = sc.transform_to(frame_altaz)
        return altaz.az.deg, altaz.alt.deg

    def local_azimuth(self, azimuth):
        if self.local_sun_azimuth is None:
            return np.asarray(azimuth)
        return (np.asarray(azimuth) - self.local_sun_azimuth + 180.0) % 360.0

    def project_global(self, ra, dec):
        ra = np.asarray(ra, dtype=float)
        dec = np.asarray(dec, dtype=float)
        projection = self.projection_combo.currentIndex()
        if projection == 0:
            return ra % 360.0, dec
        if projection == 2:
            return ra % 360.0, LAMBERT_Y_SCALE * np.sin(np.radians(dec))

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
                    if is_local and self.projection_combo.currentIndex() in (1, 3):
                        append_projected_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1], 0.5)
                    elif is_local:
                        append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
                    else:
                        append_projected_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1], 0.5) if self.projection_combo.currentIndex() in (1, 3) else append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
            self.add_curve(corridor_x, corridor_y, pg.mkPen(color=color, width=1.5, style=line_style), 6)

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
        projected = self.projection_combo.currentIndex() in (1, 3)
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
                if is_local and altitude[index] < 0.0 and altitude[index + 1] < 0.0:
                    continue
                if projected:
                    append_projected_segment(
                        boundary_x, boundary_y,
                        x_values[index], y_values[index],
                        x_values[index + 1], y_values[index + 1],
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
        projected = self.projection_combo.currentIndex() in (1, 3)
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
                    color, width, z_value = (185, 100, 40, 255), 4.0, 20
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
                        append_projected_segment(res_x, res_y, cx[i], cy[i], cx[i+1], cy[i+1])
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
                    append_projected_segment(target_x, target_y, x[i], y[i], x[i+1], y[i+1])
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
                if alt[i] >= 0 and alt[i + 1] >= 0:
                    line_x, line_y = self.project_local([az[i], az[i + 1]], [alt[i], alt[i + 1]])
                    if self.projection_combo.currentIndex() in (1, 3):
                        append_projected_segment(equator_x, equator_y, line_x[0], line_y[0], line_x[1], line_y[1])
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
                if self.projection_combo.currentIndex() in (1, 3):
                    append_projected_segment(ecliptic_x, ecliptic_y, line_x[0], line_y[0], line_x[1], line_y[1])
                else:
                    append_sky_segment(ecliptic_x, ecliptic_y, line_x[0], line_y[0], line_x[1], line_y[1])
        else:
            projected_x, projected_y = self.project_global(ecliptic_ra, ecliptic_dec)
            ecliptic_x, ecliptic_y = [], []
            for index in range(len(projected_x) - 1):
                if self.projection_combo.currentIndex() in (1, 3):
                    append_projected_segment(
                        ecliptic_x, ecliptic_y,
                        projected_x[index], projected_y[index],
                        projected_x[index + 1], projected_y[index + 1],
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
                if gal_y[i] >= 0 and gal_y[i + 1] >= 0:
                    line_x, line_y = self.project_local(
                        [gal_x[i], gal_x[i + 1]], [gal_y[i], gal_y[i + 1]])
                    if self.projection_combo.currentIndex() in (1, 3):
                        append_projected_segment(
                            galactic_x, galactic_y,
                            line_x[0], line_y[0], line_x[1], line_y[1],
                        )
                    else:
                        append_sky_segment(
                            galactic_x, galactic_y,
                            line_x[0], line_y[0], line_x[1], line_y[1],
                        )
        else:
            galactic_x, galactic_y = self.project_global(galactic.ra.deg, galactic.dec.deg)
            if self.projection_combo.currentIndex() in (1, 3):
                split_x, split_y = [], []
                for i in range(len(galactic_x) - 1):
                    if (abs(galactic_x[i + 1] - galactic_x[i]) < 0.3 and
                            abs(galactic_y[i + 1] - galactic_y[i]) < 0.3):
                        split_x.extend([galactic_x[i], galactic_x[i + 1], np.nan])
                        split_y.extend([galactic_y[i], galactic_y[i + 1], np.nan])
                galactic_x, galactic_y = split_x, split_y
            else:
                galactic_x, galactic_y = galactic_x.tolist(), galactic_y.tolist()
            if self.projection_combo.currentIndex() not in (1, 3):
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

        if is_local:
            sun_x, sun_y = self.transform_coords(sun.ra.deg, sun.dec.deg, frame_altaz)
            moon_x, moon_y = self.transform_coords(moon.ra.deg, moon.dec.deg, frame_altaz)
            sun_altitude = float(sun_y[0])
            moon_altitude = float(moon_y[0])
            self.local_sun_azimuth = float(sun_x[0])
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

        self.add_reference_grid(frame_altaz, is_local)
        self.add_reference_circles(frame_altaz, is_local)
        self.add_constellation_boundaries(frame_altaz, is_local)

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
            if is_local:
                star_az, star_alt = self.transform_coords(star_ra, star_dec, frame_altaz)
                visible = star_alt >= local_lower_altitude
                star_x, star_y = self.project_local(self.local_azimuth(star_az[visible]), star_alt[visible])
            else:
                visible = np.ones(len(star_ra), dtype=bool)
                star_x, star_y = self.project_global(star_ra, star_dec)
            visible_mag = star_mag[visible]
            visible_bayer = star_bayer[visible]
            visible_colors = star_colors[visible]
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

        if not is_local or sun_altitude >= -30:
            self.plot_widget.addItem(pg.ScatterPlotItem(x=sun_x, y=sun_y, size=18, brush='#ffcc00', pen=pg.mkPen('w', width=1.5)))
        if not is_local or moon_altitude >= 0:
            self.plot_widget.addItem(pg.ScatterPlotItem(x=moon_x, y=moon_y, size=12, brush='#eeeeee', pen=pg.mkPen('b', width=1)))

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
                    if self.projection_combo.currentIndex() in (1, 3):
                        append_projected_segment(
                            projected_x, projected_y,
                            line_x[0], line_y[0], line_x[1], line_y[1])
                    else:
                        append_sky_segment(
                            projected_x, projected_y,
                            line_x[0], line_y[0], line_x[1], line_y[1])
                x_c, y_c = np.asarray(projected_x), np.asarray(projected_y)
            else:
                raw_x, raw_y = self.project_global(ra_deg_c, dec_deg_c)
                projected_x, projected_y = [], []
                for index in range(len(raw_x) - 1):
                    if self.projection_combo.currentIndex() in (1, 3):
                        append_projected_segment(
                            projected_x, projected_y,
                            raw_x[index], raw_y[index],
                            raw_x[index + 1], raw_y[index + 1],
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

            transformed_vertices = None
            if is_local and visible_polygons:
                all_ra = np.concatenate([np.asarray(polygon[0]) for polygon in visible_polygons])
                all_dec = np.concatenate([np.asarray(polygon[1]) for polygon in visible_polygons])
                transformed_vertices = self.transform_coords(all_ra, all_dec, frame_altaz)
                vertex_offset = 0

            for ra, dec, mag in visible_polygons:
                ra_arr = np.array(ra)
                dec_arr = np.array(dec)

                if is_local:
                    vertex_count = len(ra_arr)
                    az = transformed_vertices[0][vertex_offset:vertex_offset + vertex_count]
                    alt = transformed_vertices[1][vertex_offset:vertex_offset + vertex_count]
                    vertex_offset += vertex_count
                    if np.all(alt < 0):
                        continue
                    px, py = self.project_local(self.local_azimuth(az), alt)
                else:
                    px, py = self.project_global(ra_arr, dec_arr)

                for i in range(len(px) - 1):
                    if self.projection_combo.currentIndex() in (1, 3):
                        if (abs(px[i + 1] - px[i]) < 0.35 and
                                abs(py[i + 1] - py[i]) < 0.35):
                            x_all.extend([px[i], px[i + 1], np.nan])
                            y_all.extend([py[i], py[i + 1], np.nan])
                    else:
                        append_sky_segment(x_all, y_all, px[i], py[i], px[i+1], py[i+1])

            if x_all:
                poly_item = pg.PlotCurveItem(x=np.array(x_all), y=np.array(y_all), pen=pg.mkPen(color=(0, 180, 220, 40), width=0.5), connect="finite")
                poly_item.setZValue(1)
                self.plot_widget.addItem(poly_item)

        # 4. Souhvězdí
        if self.const_lines:
            c_x, c_y = [], []
            for line in self.const_lines:
                for i in range(len(line) - 1):
                    p1, p2 = line[i], line[i+1]
                    ra1 = p1[0] % 360
                    ra2 = p2[0] % 360
                    dec1, dec2 = p1[1], p2[1]

                    if is_local:
                        az, alt = self.transform_coords(np.array([ra1, ra2]), np.array([dec1, dec2]), frame_altaz)
                        if np.all(alt < 0):
                            continue
                        plot_az = self.local_azimuth(az)
                        line_x, line_y = self.project_local(plot_az, alt)
                        if self.projection_combo.currentIndex() in (1, 3):
                            append_projected_segment(c_x, c_y, line_x[0], line_y[0], line_x[1], line_y[1])
                        else:
                            append_sky_segment(c_x, c_y, line_x[0], line_y[0], line_x[1], line_y[1])
                    else:
                        x1, y1 = self.project_global([ra1], [dec1])
                        x2, y2 = self.project_global([ra2], [dec2])
                        if self.projection_combo.currentIndex() not in (1, 3):
                            append_sky_segment(c_x, c_y, x1[0], y1[0], x2[0], y2[0])
                        else:
                            append_projected_segment(c_x, c_y, x1[0], y1[0], x2[0], y2[0])

            if c_x:
                const_item = pg.PlotCurveItem(
                    x=np.array(c_x),
                    y=np.array(c_y),
                    pen=pg.mkPen(color=(105, 220, 145, 155), width=0.9),
                    connect="finite",
                )
                const_item.setZValue(15)
                self.plot_widget.addItem(const_item)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
