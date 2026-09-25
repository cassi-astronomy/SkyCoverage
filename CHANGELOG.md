# Přehled změn

V tomto souboru jsou zaznamenány významné změny projektu SkyCoverage.

## 2026-09-25

### Přidáno

- Verzování projektu pomocí Gitu.
- Propojení s veřejným GitHub repozitářem `cassi-astronomy/SkyCoverage`.
- Automatické odeslání větve na GitHub po každém místním commitu.
- Kontrola cest a typů položek při rozbalování archivu MPC.

### Opraveno

- Odsazení v `update_axis_ticks`, které způsobovalo `IndentationError` a bránilo spuštění programu.
- Nastavení os globálního pohledu, které bylo kvůli předčasnému `return` nedosažitelné.
- Rozbalování `skycov.tgz`, které mohlo vytvářet chybnou strukturu `skycov/skycov`.
- Falešná vodorovná čára vznikající propojením konců kružnice sluneční elongace přes přechod RA 0°/360°.
- Načtení dat MPC už nepřepisuje aktuální čas datem nejnovějšího archivního souboru.

### Změněno

- Nebeský rovník je světle modrý a plný, zatímco ekliptika je výrazně oranžová a čárkovaná.
- Souřadnicový panel má pevnou velikost a samostatné řádky pro Az/Alt, RA/Dec a elongaci.
- Tlačítka soumraku počítají dnešní večer a následující ráno podle aktuálního data UTC.
