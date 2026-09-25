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
                             QPushButton, QComboBox, QDateTimeEdit, QCheckBox)
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QDateTime, QTimer
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
                            stars.append((float(row['ra']) * 15.0, float(row['dec']), magnitude))
                    except (KeyError, TypeError, ValueError):
                        continue
            self.finished.emit((stars, ""))
        except Exception as error:
            self.finished.emit(([], f"Katalog hvězd: {error}"))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MPC Sky Coverage Visualizer & Planetarium")
        self.resize(1400, 850)
        self.setStyleSheet("background-color: #111111; color: white;")

        self.selected_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "skycov")
        self.const_lines = self.load_constellation_lines()

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Panel 1: Pohled, Observatoř, Čas, Mag
        control_panel1 = QHBoxLayout()

        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItems(["Globální (RA/Dec)", "Lokální obzor (Alt/Az)"])
        self.view_mode_combo.setStyleSheet("background-color: #222; color: #00ffcc; font-weight: bold; padding: 4px;")
        self.view_mode_combo.currentIndexChanged.connect(self.on_view_change)
        control_panel1.addWidget(self.view_mode_combo)

        control_panel1.addWidget(QLabel("Projekce:"))
        self.projection_combo = QComboBox()
        self.projection_combo.addItems(["Obdélníková", "Hammer-Aitoff"])
        self.projection_combo.setStyleSheet("background-color: #222; color: white; padding: 4px;")
        self.projection_combo.currentIndexChanged.connect(self.on_projection_change)
        control_panel1.addWidget(self.projection_combo)

        control_panel1.addWidget(QLabel("Observatoř:"))
        self.obs_combo = QComboBox()
        self.obs_combo.addItems(list(OBSERVATORIES.keys()))
        self.obs_combo.setStyleSheet("background-color: #222; color: white; padding: 4px;")
        self.obs_combo.currentIndexChanged.connect(self.request_redraw)
        control_panel1.addWidget(self.obs_combo)

        control_panel1.addWidget(QLabel("Čas (UTC):"))
        self.time_edit = QDateTimeEdit(QDateTime.currentDateTimeUtc())
        self.time_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.time_edit.setStyleSheet("background-color: #222; color: white; padding: 4px;")
        self.time_edit.dateTimeChanged.connect(self.request_redraw)
        control_panel1.addWidget(self.time_edit)

        spin_style = "QSpinBox, QDoubleSpinBox { background-color: #222; color: white; border: 1px solid #555; padding-right: 10px; min-height: 25px; }"

        control_panel1.addWidget(QLabel("Dny zpět:"))
        self.days_spin = QSpinBox()
        self.days_spin.setRange(1, SKYCOV_RETENTION_DAYS)
        self.days_spin.setValue(30)
        self.days_spin.setStyleSheet(spin_style)
        control_panel1.addWidget(self.days_spin)

        control_panel1.addWidget(QLabel("Min. mag:"))
        self.mag_spin = QDoubleSpinBox()
        self.mag_spin.setRange(5.0, 25.0)
        self.mag_spin.setSingleStep(0.5)
        self.mag_spin.setValue(15.0)
        self.mag_spin.setStyleSheet(spin_style)
        control_panel1.addWidget(self.mag_spin)

        self.btn_load = QPushButton("Načíst data")
        self.btn_load.setStyleSheet("background-color: #008844; color: white; padding: 6px 12px; font-weight: bold;")
        self.btn_load.clicked.connect(self.start_processing)
        control_panel1.addWidget(self.btn_load)

        self.kreutz_check = QCheckBox("Kreutzův koridor I/II")
        self.kreutz_check.setToolTip(
            "Model koridorů Kreutz I a II pro 5, 10, 15, 20, 30 a 45 dní před přísluním"
        )
        self.kreutz_check.setStyleSheet("color: #ff9966; font-weight: bold;")
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
        self.coord_label.setStyleSheet(
            "color: #00ffcc; font-weight: bold; "
            "font-family: monospace; font-size: 11px;"
        )
        control_panel1.addStretch()
        control_panel1.addWidget(self.coord_label)

        main_layout.addLayout(control_panel1)

        # Panel 2: Tlačítka soumraků
        control_panel2 = QHBoxLayout()
        control_panel2.addWidget(QLabel("Rychlé nastavení soumraku:"))

        btn_style = "QPushButton { background-color: #2a2a3a; color: #ffcc00; border: 1px solid #444; padding: 4px 8px; font-size: 11px; } QPushButton:hover { background-color: #3a3a5a; }"

        self.btn_naut_eve = QPushButton("Večerní Nautický (-12°)")
        self.btn_naut_eve.setStyleSheet(btn_style)
        self.btn_naut_eve.clicked.connect(lambda: self.set_twilight(-12.0, evening=True))
        control_panel2.addWidget(self.btn_naut_eve)

        self.btn_astro_eve = QPushButton("Večerní Astro Noc (-18°)")
        self.btn_astro_eve.setStyleSheet(btn_style)
        self.btn_astro_eve.clicked.connect(lambda: self.set_twilight(-18.0, evening=True))
        control_panel2.addWidget(self.btn_astro_eve)

        self.btn_astro_morn = QPushButton("Ranní Astro Noc (-18°)")
        self.btn_astro_morn.setStyleSheet(btn_style)
        self.btn_astro_morn.clicked.connect(lambda: self.set_twilight(-18.0, evening=False))
        control_panel2.addWidget(self.btn_astro_morn)

        self.btn_naut_morn = QPushButton("Ranní Nautický (-12°)")
        self.btn_naut_morn.setStyleSheet(btn_style)
        self.btn_naut_morn.clicked.connect(lambda: self.set_twilight(-12.0, evening=False))
        control_panel2.addWidget(self.btn_naut_morn)

        for label, minutes in [("+30 min", 30), ("+1 h", 60), ("-1 h", -60), ("-30 min", -30)]:
            button = QPushButton(label)
            button.setStyleSheet(btn_style)
            button.clicked.connect(lambda checked=False, delta=minutes: self.shift_time(delta))
            control_panel2.addWidget(button)

        control_panel2.addStretch()
        main_layout.addLayout(control_panel2)

        # Plot
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('#050508')
        self.plot_widget.getViewBox().invertX(True)
        self.plot_widget.scene().sigMouseMoved.connect(self.on_mouse_move)
        main_layout.addWidget(self.plot_widget)

        self.status_label = QLabel("Připraveno.")
        main_layout.addWidget(self.status_label)

        self.polygons = []
        self.stars = []
        self.ref_time = None
        self.local_sun_azimuth = None
        self.uncovered_cache = {}
        self.redraw_timer = QTimer(self)
        self.redraw_timer.setSingleShot(True)
        self.redraw_timer.timeout.connect(self.redraw)
        self.reload_timer = QTimer(self)
        self.reload_timer.setSingleShot(True)
        self.reload_timer.timeout.connect(self.start_processing)
        self.axis_tick_timer = QTimer(self)
        self.axis_tick_timer.setSingleShot(True)
        self.axis_tick_timer.timeout.connect(self.update_axis_ticks)
        self.plot_widget.getViewBox().sigRangeChanged.connect(
            self.request_axis_tick_update
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

    def load_constellation_lines(self):
        cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "astro_cache")
        os.makedirs(cache_dir, exist_ok=True)
        json_path = os.path.join(cache_dir, "constellations.lines.json")

        if not os.path.exists(json_path):
            url = "https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json"
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=10) as resp:
                    with open(json_path, 'wb') as f:
                        f.write(resp.read())
            except Exception:
                return []

        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
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
        except Exception:
            return []

    def set_twilight(self, target_alt_deg, evening=True):
        """Nastaví dnešní večerní nebo následující ranní soumrak v UTC."""
        obs_location = OBSERVATORIES[self.obs_combo.currentText()]

        base_date = datetime.now(timezone.utc).replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        if evening:
            search_date = base_date
            start_hour, end_hour = 12, 24
        else:
            search_date = base_date + timedelta(days=1)
            start_hour, end_hour = 0, 12
        times = [search_date + timedelta(minutes=m) for m in range(start_hour * 60, end_hour * 60 + 1)]

        astropy_times = Time(times)
        frame = AltAz(obstime=astropy_times, location=obs_location)
        suns = get_sun(astropy_times).transform_to(frame)
        alts = suns.alt.deg

        crossings = np.where((alts[:-1] - target_alt_deg) * (alts[1:] - target_alt_deg) <= 0)[0]
        if len(crossings):
            best_idx = crossings[np.argmin(np.abs(alts[crossings] - target_alt_deg))]
        else:
            best_idx = int(np.argmin(np.abs(alts - target_alt_deg)))
        best_dt = times[best_idx].replace(tzinfo=timezone.utc)

        self.time_edit.blockSignals(True)
        self.time_edit.setDateTime(QDateTime(best_dt))
        self.time_edit.blockSignals(False)
        self.redraw()

    def shift_time(self, minutes):
        current = self.time_edit.dateTime().toPyDateTime()
        self.time_edit.setDateTime(QDateTime(current + timedelta(minutes=minutes)))

    def on_view_change(self):
        is_local = self.view_mode_combo.currentIndex() == 1
        local_lower_altitude = -30.0
        self.plot_widget.clear()
        if is_local:
            self.plot_widget.getViewBox().invertX(False)
            if self.projection_combo.currentIndex() == 1:
                self.plot_widget.setXRange(-3, 3)
                self.plot_widget.setYRange(-1.5, 1.5)
                self.plot_widget.setLabel('bottom', 'Relativní azimut (střed = Slunce)')
                self.plot_widget.setLabel('left', 'Výška Alt (stupně)')
            else:
                self.plot_widget.setXRange(0, 360)
                self.plot_widget.setYRange(0, 90)
                self.plot_widget.setLabel('bottom', 'Azimut vůči Slunci (střed = směr Slunce)')
                self.plot_widget.setLabel('left', 'Výška nad obzorem Alt (deg)')
        else:
            self.plot_widget.getViewBox().invertX(True)
            if self.projection_combo.currentIndex() == 1:
                self.plot_widget.setXRange(-3, 3)
                self.plot_widget.setYRange(-1.5, 1.5)
                self.plot_widget.setLabel('bottom', 'Rektascenze RA (hodiny)')
                self.plot_widget.setLabel('left', 'Deklinace Dec (stupně)')
            else:
                self.plot_widget.setXRange(0, 360)
                self.plot_widget.setYRange(-90, 90)
                self.plot_widget.setLabel('bottom', 'Rektascenze RA (hodiny)')
                self.plot_widget.setLabel('left', 'Deklinace Dec (stupně)')
        self.redraw()
        self.update_axis_ticks()

    def on_projection_change(self):
        self.on_view_change()

    def update_axis_ticks(self):
        is_local = self.view_mode_combo.currentIndex() == 1
        hammer = self.projection_combo.currentIndex() == 1
        bottom = self.plot_widget.getAxis('bottom')
        left = self.plot_widget.getAxis('left')
        if is_local:
            bottom.setLabel('Relativní azimut (střed = Slunce)' if not hammer else 'Relativní azimut (Hammer-Aitoff)')
            left.setLabel('Výška Alt (stupně)')
            x_range, y_range = self.plot_widget.getViewBox().viewRange()
            x_width = abs(x_range[1] - x_range[0])
            y_height = abs(y_range[1] - y_range[0])
            if hammer:
                tick_step = 10 if x_width < 4.5 or y_height < 2.25 else 30
                center_x = (x_range[0] + x_range[1]) / 2.0
                center_y = (y_range[0] + y_range[1]) / 2.0
                if center_x ** 2 / 8.0 + center_y ** 2 / 2.0 <= 1.0:
                    center_longitude, center_altitude = self.inverse_local_hammer(
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
                tick_step = 10 if x_width < 240.0 or y_height < 80.0 else 30
                az_values = np.arange(0, 360, tick_step)
                az_positions = self.local_azimuth(az_values)
                bottom.setTicks([[
                    (float(position), f'{value}°')
                    for position, value in zip(az_positions, az_values)
                ]])
                alt_values = np.arange(-90, 91, tick_step)
                left.setTicks([[
                    (float(value), f'{value}°') for value in alt_values
                ]])
            return

        bottom.setLabel('Rektascenze RA (hodiny)')
        left.setLabel('Deklinace Dec (stupně)')
        ra_values = np.arange(0, 360, 15)
        ra_minor = np.arange(0, 360, 7.5)
        if hammer:
            ra_positions = self.project_global(ra_values, np.zeros_like(ra_values))[0]
            ra_minor_positions = self.project_global(ra_minor, np.zeros_like(ra_minor))[0]
            dec_values = np.arange(-80, 81, 20)
            dec_minor = np.arange(-90, 91, 5)
            dec_positions = self.project_global(np.full_like(dec_values, 180), dec_values)[1]
            dec_minor_positions = self.project_global(np.full_like(dec_minor, 180), dec_minor)[1]
        else:
            ra_positions = ra_values
            ra_minor_positions = ra_minor
            dec_values = np.arange(-80, 81, 20)
            dec_minor = np.arange(-90, 91, 5)
            dec_positions = dec_values
            dec_minor_positions = dec_minor
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
            if self.projection_combo.currentIndex() == 1:
                if x ** 2 / 8.0 + y ** 2 / 2.0 > 1.0:
                    self.show_invalid_coordinates()
                    return
                local_x, altitude = self.inverse_local_hammer(x, y)
                azimuth = (local_x + (self.local_sun_azimuth or 0)) % 360
            else:
                azimuth = (x + (self.local_sun_azimuth or 0) - 180) % 360
                altitude = y
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
        if self.projection_combo.currentIndex() == 0:
            return ra % 360.0, dec

        longitude = np.radians((ra - 180.0 + 180.0) % 360.0 - 180.0)
        latitude = np.radians(dec)
        denominator = np.sqrt(1.0 + np.cos(latitude) * np.cos(longitude / 2.0))
        return (
            2.0 * np.sqrt(2.0) * np.cos(latitude) * np.sin(longitude / 2.0) / denominator,
            np.sqrt(2.0) * np.sin(latitude) / denominator,
        )

    def inverse_global(self, x, y):
        if self.projection_combo.currentIndex() == 0:
            return np.asarray(x), np.asarray(y)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.sqrt(np.maximum(0.0, 1.0 - (x ** 2 / 16.0) - (y ** 2 / 4.0)))
        longitude = 2.0 * np.arctan2(z * x, 2.0 * (2.0 * z ** 2 - 1.0))
        latitude = np.arcsin(np.clip(z * y, -1.0, 1.0))
        return (np.degrees(longitude) + 180.0) % 360.0, np.degrees(latitude)

    def project_local(self, azimuth, altitude):
        azimuth = np.asarray(azimuth, dtype=float)
        altitude = np.asarray(altitude, dtype=float)
        if self.projection_combo.currentIndex() == 0:
            return azimuth, altitude
        longitude = np.radians((azimuth - 180.0 + 180.0) % 360.0 - 180.0)
        latitude = np.radians(altitude)
        denominator = np.sqrt(1.0 + np.cos(latitude) * np.cos(longitude / 2.0))
        return (
            2.0 * np.sqrt(2.0) * np.cos(latitude) * np.sin(longitude / 2.0) / denominator,
            np.sqrt(2.0) * np.sin(latitude) / denominator,
        )

    def inverse_local_hammer(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
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
        days_values = [5, 10, 15, 20, 30, 45]
        for population, (color, line_style) in styles.items():
            corridor_x, corridor_y = [], []
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
                for index in range(len(points_x) - 1):
                    if is_local and self.projection_combo.currentIndex() == 1:
                        append_projected_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1], 0.5)
                    elif is_local:
                        append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
                    else:
                        append_projected_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1], 0.5) if self.projection_combo.currentIndex() == 1 else append_sky_segment(corridor_x, corridor_y, points_x[index], points_y[index], points_x[index + 1], points_y[index + 1])
            self.add_curve(corridor_x, corridor_y, pg.mkPen(color=color, width=1.5, style=line_style), 6)









    def add_curve(self, x_values, y_values, pen, z_value=0):
        if len(x_values) == 0:
            return
        item = pg.PlotCurveItem(x=np.array(x_values), y=np.array(y_values), pen=pen, connect="finite")
        item.setZValue(z_value)
        self.plot_widget.addItem(item)

    def add_reference_grid(self, frame_altaz, is_local):
        hammer = self.projection_combo.currentIndex() == 1
        if is_local:
            major_x, major_y, minor_x, minor_y = [], [], [], []
            for azimuth in range(0, 360, 10):
                plot_azimuth = self.local_azimuth(azimuth)
                if hammer:
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
                    color, width, z_value = (170, 105, 45, 230), 2.5, 12
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
                    if hammer:
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
                if hammer:
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
                    if self.projection_combo.currentIndex() == 1:
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
                if self.projection_combo.currentIndex() == 1:
                    append_projected_segment(ecliptic_x, ecliptic_y, line_x[0], line_y[0], line_x[1], line_y[1])
                else:
                    append_sky_segment(ecliptic_x, ecliptic_y, line_x[0], line_y[0], line_x[1], line_y[1])
        else:
            projected_x, projected_y = self.project_global(ecliptic_ra, ecliptic_dec)
            ecliptic_x, ecliptic_y = [], []
            for index in range(len(projected_x) - 1):
                if self.projection_combo.currentIndex() == 1:
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
                    append_projected_segment(galactic_x, galactic_y, line_x[0], line_y[0], line_x[1], line_y[1])
        else:
            galactic_x, galactic_y = self.project_global(galactic.ra.deg, galactic.dec.deg)
            if self.projection_combo.currentIndex() == 1:
                split_x, split_y = [], []
                for i in range(len(galactic_x) - 1):
                    if (abs(galactic_x[i + 1] - galactic_x[i]) < 0.3 and
                            abs(galactic_y[i + 1] - galactic_y[i]) < 0.3):
                        split_x.extend([galactic_x[i], galactic_x[i + 1], np.nan])
                        split_y.extend([galactic_y[i], galactic_y[i + 1], np.nan])
                galactic_x, galactic_y = split_x, split_y
            else:
                galactic_x, galactic_y = galactic_x.tolist(), galactic_y.tolist()
            if self.projection_combo.currentIndex() == 0:
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
            else:
                self.plot_widget.setYRange(lower_altitude, 90, padding=0)
        else:
            self.local_sun_azimuth = None
            sun_altitude = float(sun.dec.deg)
            moon_altitude = float(moon.dec.deg)
            sun_x, sun_y = self.project_global([sun.ra.deg], [sun.dec.deg])
            moon_x, moon_y = self.project_global([moon.ra.deg], [moon.dec.deg])

        self.add_reference_grid(frame_altaz, is_local)
        self.add_reference_circles(frame_altaz, is_local)

        if self.stars:
            star_ra = np.array([star[0] for star in self.stars])
            star_dec = np.array([star[1] for star in self.stars])
            star_mag = np.array([star[2] for star in self.stars])
            if is_local:
                star_az, star_alt = self.transform_coords(star_ra, star_dec, frame_altaz)
                visible = star_alt >= local_lower_altitude
                star_x, star_y = self.project_local(self.local_azimuth(star_az[visible]), star_alt[visible])
            else:
                visible = np.ones(len(star_ra), dtype=bool)
                star_x, star_y = self.project_global(star_ra, star_dec)
            star_size = np.clip(7.0 - star_mag[visible], 1.5, 6.0)
            self.plot_widget.addItem(pg.ScatterPlotItem(
                x=star_x, y=star_y, size=star_size, brush=(220, 230, 255, 190), pen=None))

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
                    if self.projection_combo.currentIndex() == 1:
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
                    if self.projection_combo.currentIndex() == 1:
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
                    if self.projection_combo.currentIndex() == 1:
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
                        if self.projection_combo.currentIndex() == 1:
                            append_projected_segment(c_x, c_y, line_x[0], line_y[0], line_x[1], line_y[1])
                        else:
                            append_sky_segment(c_x, c_y, line_x[0], line_y[0], line_x[1], line_y[1])
                    else:
                        x1, y1 = self.project_global([ra1], [dec1])
                        x2, y2 = self.project_global([ra2], [dec2])
                        if self.projection_combo.currentIndex() == 0:
                            append_sky_segment(c_x, c_y, x1[0], y1[0], x2[0], y2[0])
                        else:
                            append_projected_segment(c_x, c_y, x1[0], y1[0], x2[0], y2[0])

            if c_x:
                const_item = pg.PlotCurveItem(x=np.array(c_x), y=np.array(c_y), pen=pg.mkPen(color=(255, 255, 255, 120), width=0.8), connect="finite")
                const_item.setZValue(15)
                self.plot_widget.addItem(const_item)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
