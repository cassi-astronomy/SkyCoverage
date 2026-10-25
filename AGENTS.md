# SkyCoverage – pokyny pro úpravy

## Projekt

- Hlavní desktopová aplikace je `coverage.py` (PyQt6 + pyqtgraph).
- Čisté mapové výpočty patří do `sky_geometry.py`.
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

## Organizace kódu

- Numerické funkce bez vazby na Qt umisťujte do samostatného modulu a testujte je přímo.
- `MainWindow` má koordinovat stav a vykreslení; nepřidávejte do ní další obecné parsery nebo matematické algoritmy.
- Drahé hodnoty závislé jen na načteném katalogu ukládejte do cache při načtení katalogu, ne při každém překreslení.
- Zachovejte veřejné importy z `coverage.py`, které používají stávající testy nebo externí skripty.
- Zápisy uživatelských dat provádějte atomicky přes dočasný soubor a `os.replace`.

## Ověření změn

Po každé změně spusťte:

```powershell
python -m py_compile coverage.py observation_planner.py sky_geometry.py orbital_mechanics.py
pytest -q
git diff --check
```

U změn vykreslování navíc proveďte offscreen náhled alespoň v rovníkovém a horizontálním režimu. Do repozitáře nepřidávejte diagnostické obrázky ani stažené katalogy.
