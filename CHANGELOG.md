# Přehled změn

V tomto souboru jsou zaznamenány významné změny projektu SkyCoverage.

## 2026-09-25

### Přidáno

- Verzování projektu pomocí Gitu.
- Propojení s veřejným GitHub repozitářem `cassi-astronomy/SkyCoverage`.
- Automatické odeslání větve na GitHub po každém místním commitu.
- Kontrola cest a typů položek při rozbalování archivu MPC.
- Elongační kružnice 70° kolem Slunce.

### Opraveno

- Odsazení v `update_axis_ticks`, které způsobovalo `IndentationError` a bránilo spuštění programu.
- Nastavení os globálního pohledu, které bylo kvůli předčasnému `return` nedosažitelné.
- Rozbalování `skycov.tgz`, které mohlo vytvářet chybnou strukturu `skycov/skycov`.
- Falešná vodorovná čára vznikající propojením konců kružnice sluneční elongace přes přechod RA 0°/360°.
- Načtení dat MPC už nepřepisuje aktuální čas datem nejnovějšího archivního souboru.
- Pohyb myši mimo platnou oblast lokální nebo globální projekce už neposílá Astropy souřadnice s výškou mimo rozsah −90° až +90°.
- Inverzní Hammerův převod používá správné měřítko, takže Az/Alt a RA/Dec pod kurzorem souhlasí s vykreslenou mřížkou.

### Změněno

- Nebeský rovník je světle modrý a plný, zatímco ekliptika je výrazně oranžová a čárkovaná.
- Souřadnicový panel má pevnou velikost a samostatné řádky pro Az/Alt, RA/Dec a elongaci.
- Tlačítka soumraku počítají dnešní večer a následující ráno podle aktuálního data UTC.
- Data Sky Coverage se načítají pouze z ručně rozbalené místní složky; program už nestahuje `skycov.tgz` z MPC.
- Při načtení se uchovává posledních 62 dní od nejnovějšího `.DAT` souboru a starší soubory se mažou.
- Po vyčištění dat se odstraní také prázdné podadresáře observatoří; kořenová složka `skycov` vždy zůstane zachovaná.
- `.DAT` soubory bez podporovaného data v názvu se odstraní, protože je program nedokáže zařadit ani zobrazit.
- Lokální mřížka obsahuje jemné čáry po 10° a výraznější čáry po 30°; popisky os se při přiblížení automaticky přepnou z 30° na 10°.
- Výškové kružnice po 10° a 20° jsou téměř stejně jasné jako kružnice po 30°, aby zůstaly čitelné nad vykresleným pokrytím.
- Kreutzovy koridory sahají až k bodu 5 dní před přísluním.
- Ovladač „Dny zpět“ je omezený na uchovávaných 62 dní a tlačítko se jmenuje „Načíst data“.
