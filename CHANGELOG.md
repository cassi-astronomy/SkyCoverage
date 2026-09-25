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
- Data Sky Coverage se načítají pouze z ručně rozbalené místní složky; program už nestahuje `skycov.tgz` z MPC.
- Při načtení se uchovává posledních 62 dní od nejnovějšího `.DAT` souboru a starší soubory se mažou.
- Po vyčištění dat se odstraní také prázdné podadresáře observatoří; kořenová složka `skycov` vždy zůstane zachovaná.
- `.DAT` soubory bez podporovaného data v názvu se odstraní, protože je program nedokáže zařadit ani zobrazit.
- Ovladač „Dny zpět“ je omezený na uchovávaných 62 dní a tlačítko se jmenuje „Načíst data“.
