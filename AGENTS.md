# SkyCoverage – pokyny pro úpravy

## Projekt

- Hlavní desktopová aplikace je `coverage.py` (PyQt6 + pyqtgraph).
- Čisté mapové výpočty patří do `sky_geometry.py`.
- Výběr astronomických časů patří do `astronomy_time.py`.
- Čisté výpočty drah patří do `orbital_mechanics.py`.
- Plánování, profily sestav a export plánů patří do `observation_planner.py`.
- Uživatelské rozhraní je v češtině; nové texty udržujte ve stejném jazyce.

## Neměnné vlastnosti

- Bez výslovného požadavku neměňte vzhled, barvy, pořadí ovládacích prvků ani výchozí nastavení.
- Všechny vrstvy musí fungovat v globálním i horizontálním zobrazení a ve všech podporovaných projekcích.
- Čáry nesmějí vytvářet falešné spojnice přes šev 0°/360°.
- Zachovejte formát `planner_data.json`, archivních plánů a textových exportů.
- Soubory v `skycov/`, `astro_cache/`, `observation_plans/` a `planner_data.json` jsou lokální data uživatele. Testy je nesmějí měnit ani mazat.
- Síťové katalogy se při běžném startu mohou obnovovat, ale automatické testy musí být nezávislé na síti.
- Výřez mapy se nastavuje automaticky jen po přepnutí pohledu nebo projekce a po resetu. Změna času, observatoře ani vrstev nesmí zahodit ruční posun a přiblížení.
- Čas mapy je vždy skutečný okamžik v UTC bez ohledu na zobrazené pásmo (`UTC`/`SEČ`/`SELČ`). Pro výpočty používejte `map_time()` nebo `qdatetime_to_time()`, nikdy `Time(...dateTime().toPyDateTime())`, protože `toPyDateTime()` zahodí časové pásmo.
- Úklid v `skycov/` nesmí mít pevný horní rok. Platnost roku z názvu souboru určuje `mpc_file_year_is_plausible()`.

## Organizace kódu

- Numerické funkce bez vazby na Qt umisťujte do samostatného modulu a testujte je přímo.
- `MainWindow` má koordinovat stav a vykreslení; nepřidávejte do ní další obecné parsery nebo matematické algoritmy.
- Drahé hodnoty závislé jen na načteném katalogu ukládejte do cache při načtení katalogu, ne při každém překreslení.
- Zachovejte veřejné importy z `coverage.py`, které používají stávající testy nebo externí skripty.
- Zápisy uživatelských dat provádějte atomicky přes dočasný soubor a `os.replace`.

## Výkon vykreslování

Mapa musí zůstat plynulá při posunu a přiblížení. Platí tato pravidla:

- **Žádné astropy v horkých cestách.** Obsluha myši, náhled mozaiky, `redraw()` ani `refresh_view_layers()` nesmějí volat `SkyCoord(...).transform_to(...)`, `directional_offset_by` ani `get_sun`/`get_body` po jednotlivých bodech nebo ve smyčkách. Jedno volání astropy stojí 5–10 ms.
- **ICRS ↔ AltAz** převádějte přes `transform_coords()`, `horizontal_to_icrs()` nebo `project_sky()`. Používají rotační matici z `horizontal_rotation()`, která je v cache pro každý čas a observatoř (chyba pod 0,01°). Astropy zůstává jen pro Slunce, Měsíc a planety (topocentrická paralaxa Měsíce) a pro jednorázové výpočty. Jejich výsledky se ukládají do cache podle času a observatoře.
- **Posuny a obrysy polí** počítejte vektorově přes `offset_coordinates()` a `field_outline_arrays()`. `offset_coordinate()` a `field_outline()` zůstávají jen jako kompatibilní obálky.
- **Čáry** stavějte vektorově přes `MainWindow.segments()`, `pair_segments()` nebo `sky_geometry.polyline_segments()`, `sky_segments()` a `projected_segments()`, ne smyčkou s `append_sky_segment()`/`append_projected_segment()`. Každá `PlotCurveItem` s NaN oddělovači musí mít `connect="finite"`, jinak pyqtgraph spojí body přes mezery a vzniknou falešné spojnice.
- **Katalogy** převádějte na numpy pole jednou po načtení (`build_star_arrays`, `build_polygon_arrays`, `build_constellation_segments`, `build_boundary_arrays`). `redraw()` je pouze promítá.
- **Vrstvy závislé na výřezu** patří do `refresh_view_layers()`, ne do `redraw()`. Jde o rastry (stínování pod obzorem, výplň Mléčné dráhy), značky ořezané na výřez (DSO, komety, planetky) a Bayerovy popisky. Spouští se se zpožděním po ustálení posunu nebo zoomu. Rastry pokrývají výřez s okrajem (`raster_extent()`, `raster_grid()`), aby se při posunu neodhalil nevykreslený okraj. Značky se ořezávají přes `culling_range()`.
- **Popisky** (`TextItem`) při obnově výřezu znovu používejte (`setText`/`setPos`), nevytvářejte je pokaždé znovu.
- **`on_mouse_move()`** musí trvat pod 1 ms. Je omezený přes `pg.SignalProxy` a během tažení se stisknutým tlačítkem nic nepočítá. Poloha Slunce pro elongaci se bere z `redraw()` (`sun_radec`).
- **`DeviceCoordinateCache`** používejte jen tam, kde nemění vzhled. U čar s velmi nízkým krytím (zorná pole MPC) mění odstín kvantizací, proto je zapnutá jen v MER. U velkých položek přes celou oblohu naopak zpomaluje.
- **Orientační limity** pro offscreen měření v okně 1600×1000: `redraw()` pod 300 ms, `refresh_view_layers()` pod 80 ms, `on_mouse_move()` včetně náhledu mozaiky 3×3 pod 5 ms. Když změna limit překročí, najděte příčinu profilováním (`cProfile`) dřív, než ji zařadíte.

## Ověření změn

Po každé změně spusťte:

```powershell
python -m py_compile coverage.py observation_planner.py sky_geometry.py orbital_mechanics.py astronomy_time.py
pytest -q
git diff --check
```

U změn vykreslování navíc proveďte offscreen náhled alespoň v rovníkovém a horizontálním režimu (`QT_QPA_PLATFORM=offscreen`). Porovnejte ho s předchozí verzí, aby se nezměnil vzhled. Projděte všechny kombinace pohledu a projekce se zapnutými vrstvami, včetně posunu, zoomu a náhledu mozaiky, a zkontrolujte, že nevznikají výjimky. U změn výkonu změřte časy podle limitů v části *Výkon vykreslování*.

Offscreen testy, které vytvářejí `MainWindow`, spouštějte jen nad kopií adresáře aplikace. Načítání dat totiž uklízí soubory ve `skycov/` a mění `planner_data.json`.

Do repozitáře nepřidávejte diagnostické obrázky ani stažené katalogy.
