# Přehled změn

## 2026-10-30

### Změněno

- Při zavírání okna s neuloženými poli plánu se aplikace zeptá, zda plán uložit, zavřít bez uložení, nebo zavírání zrušit. Za uložený se plán považuje po uložení do archivu nebo po načtení z archivu, dokud se pole nezmění.
- Komety se na mapě kreslí znakem komety (hlava s ohonem) místo modrého kosočtverce. Popisek při přiblížení už obsahuje jen jméno a jasnost, bez znaku navíc.
- Planetky se zobrazují do 14 mag (dříve 12 mag). Dráhy se berou z MPCORB, který se po stažení ořízne jen na planetky, které 14 mag mohou dosáhnout. Bez něj se použije původní roční seznam jasných planetek MPC (asi do 12 mag).
- Komety a planetky se už nestahují automaticky při spuštění. Na kartě `Data MPC` jsou tlačítka `Stáhnout komety` a `Stáhnout planetky` s průběhem stahování v procentech a s údajem, kdy byla data stažena. Mapa se po stažení sama obnoví.
- Planetky mají sytě žlutou barvu, aby se nepletly s oranžovými archivními poli. Při přiblížení pod 55° nesou číslo, pod 25° také jméno a jasnost (nejvýš 40 nejjasnějších popisků ve výřezu). Nečíslované planetky mají provizorní označení.
- Polohy planetek se počítají vektorově v numpy pro celý katalog najednou.

- Mapa je při posunu a přiblížení výrazně plynulejší. Překreslení mapy se zkrátilo z 1–2,5 s na zhruba 0,2 s. Hlavní brzdou byly obrysy archivních a plánovaných polí, které se počítaly přes astropy bod po bodu (asi 30 ms na pole). Teď se počítají vektorově v numpy.
- Převod ICRS ↔ Alt/Az používá matici v cache pro zvolený čas a observatoř místo opakovaných volání astropy (rozdíl pod 0,01°). Měsíc a planety se dál počítají topocentricky přes astropy, jen jednou pro daný čas.
- Údaj souřadnic pod kurzorem se počítá v řádu desetin milisekundy, nejvýš 60× za sekundu, a během tažení mapy vůbec. Náhled mozaiky pod kurzorem se zrychlil zhruba stokrát (3×3 pole z ~450 ms na ~3 ms).
- Výplň Mléčné dráhy, stínování pod obzorem, DSO, komety, planetky a Bayerovy popisky se po posunu nebo zoomu obnoví samostatně, bez celého překreslení. Rastry pokrývají i okolí výřezu, takže se při posunu neodhalí nevykreslený okraj. DSO a planetky se nově objeví i po přiblížení tlačítky.
- Katalogy hvězd, zorných polí MPC a souhvězdí se převádějí na numpy pole jednou po načtení. Čáry se skládají vektorově.

### Opraveno

- V globálním zobrazení se komety skrývaly podle toho, zda jsou nad obzorem zvolené observatoře, takže se při změně místa objevovaly a mizely. Obzor se teď uplatní jen v horizontálním zobrazení, stejně jako u planetek.
- Dithering prvního pole v expozičním skriptu byl jen kladným směrem (`rand(0, d)`), takže pole bylo v průměru posunuté o polovinu ditheru. Nyní je u všech polí rozsah symetrický kolem středu (`±d/2`) podle ditheru z profilu sestavy.
- Při zobrazení v `SEČ`/`SELČ` se celá obloha počítala pro čas posunutý o 1–2 hodiny, protože se zobrazený místní čas bral jako UTC. Export plánu nyní v řádku „Čas mapy (UTC)“ uvádí skutečný UTC čas.
- Úklid dat MPC by od roku 2031 považoval všechny soubory za neplatné a smazal je (pevný horní limit roku 2030). Limit se teď odvíjí od aktuálního roku.
- Změna času v lokálním pohledu už neresetuje ruční posun a přiblížení mapy.
- Elongační kružnice se kreslily s falešnými spojnicemi přes přerušení (chyběl `connect="finite"`), což bylo vidět hlavně ve stereografické a MER projekci.
- Kreutzovy koridory se počítaly přes chybnou dvojitou změnu počátku souřadnic (barycentrum/Slunce). Teď se počítají přímo a jsou výrazně rychlejší.
- Signály změny času a observatoře už nejsou na překreslení napojené dvakrát.

## 2026-10-25

### Přidáno

- Vnější hranice Mléčné dráhy podle katalogu José R. Vieiry s jemnou výplní (90 % průhlednost) a vlastním přepínačem vrstvy; vnitřní jasové struktury se nevykreslují.
- Měsíc ve vrstvě planet a údaj o jeho výšce a osvětlené fázi v horním souřadnicovém panelu.
- Zřetelně světlejší modré pozadí mapy, když je Měsíc nad obzorem a nachází se více než tři dny od novu.
- Jasné planetky do 12 mag z ročního katalogu MPC; při přiblížení pod 55° se zobrazují s decentním číselným popiskem.
- Archiv plánů má tlačítko `Ukázat zapnuté`, které přizpůsobí výřez mapy všem označeným plánům, i když leží ve vzdálených částech oblohy.
- Ovládání času obsahuje také jemné kroky `+10 min` a `-10 min`.
- Čas lze zobrazit v `UTC`, `SELČ` (UTC+2), nebo `SEČ` (UTC+1); zvolené pásmo se obnoví při dalším spuštění.

### Změněno

- Vodorovná mozaika s EQ kamerou je první a výchozí volbou orientace.
- Galaktický rovník se už nevykresluje; nahrazují jej katalogové obrysy Mléčné dráhy.
- Aktivní i archivní plánovaná pole používají jednu společnou oranžovou barvu; archivní pole mají 50% krytí a pole načteného archivního plánu se při současně zapnutém archivu nevykreslují podruhé.
- Mapové projekce, geometrie vrstev a výpočty drah jsou oddělené od GUI; obálka Mléčné dráhy se připravuje pouze jednou po načtení katalogu.
- Hranice souhvězdí, kružnice výšky nad obzorem a zorná pole MPC mají odlišnější barvy, styly a hierarchii; hustá vrstva zorných polí je jemnější a výškové kružnice po 10° jsou zřetelnější než azimutální mřížka.

### Opraveno

- Ranní soumrak vybírá nejbližší budoucí ranní průchod zvolené výšky Slunce. Po půlnoci proto zůstává u právě nadcházejícího rána a nepřeskakuje automaticky na další den.
- Archiv plánů vyhrazuje hlavní prostor názvu plánu, aby bylo čitelné datum pozorování na jeho začátku. Datum uložení zobrazuje jen den; celý obsah buněk zůstává dostupný v tooltipu.

## 2026-10-23

### Opraveno

- Horní ovládací karta na široké obrazovce využívá volné místo ve dvou řádcích a v užším okně se automaticky přeskupí, takže zbytečně neposouvá mapu dolů.
- Podbarvení pod obzorem ve stereografické projekci používá hladkou vektorovou kružnici bez viditelných rastrových zubů.
- Texty u číselných polí plánovače a editoru sestav jsou oddělené od hodnot, takže se při užším panelu neschovávají pod ovládací šipky ani neořezávají zleva.

### Přidáno

- Souřadnicový panel ukazuje výšku Slunce a zda je pro zvolený čas a observatoř nad, pod, nebo právě na obzoru.
- Profil sestavy obsahuje dobu vyčtení snímku a průměrnou dobu přejezdu montáže mezi poli; plánovač z nich, počtu polí, cyklů a délky expozice průběžně odhaduje celkovou délku sekvence.

### Změněno

- Vrstva planet zobrazuje Merkur až Neptun krátkými popisky `Me`, `Ve`, `Ma`, `Ju`, `Sa`, `Ur`, `Ne`; Pluto se nevykresluje.

## 2026-10-22

### Přidáno

- Plánovač pozorování s umisťováním jednotlivých polí i mozaik kliknutím do mapy.
- Profily vybavení s vlastní šířkou a výškou zorného pole, rotací kamery a ditheringem.
- Generování expozičního skriptu kompatibilního s původním `off_CO.py`, včetně filtru, délky expozice a počtu cyklů.
- Lokální archiv plánů, jehož jednotlivé vrstvy lze v mapě zapínat a vypínat.
- Indikace data a stáří nejnovějších MPC dat a tlačítko pro bezpečnou aktualizaci z oficiálního balíčku.
- Samostatné pracovní záložky pro velkou mapu pokrytí a plánování pozorování.
- Tlačítka pro přiblížení, oddálení a obnovení celého pohledu bez kolečka myši.
- Širokoúhlá MER Mercatorova projekce s výchozím zorným polem 116° a volitelným středem.
- Živý čárkovaný náhled zorného pole nebo celé mozaiky pod kurzorem před umístěním kliknutím.
- MER mapa se při tažení posouvá bez nákladného přepočtu projekce a všech MPC polí.
- Vykreslení více než 100 000 MPC polí je v MER režimu vektorizované a oříznuté na aktuální výřez.
- V lokálním MER pohledu je zakřivený obzor a oblast pod ním je světle hnědá s 90% průhledností; hvězdy, Slunce, Měsíc i plánovací vrstvy zůstávají viditelné.
- MER je výchozí projekce; lokální pohled se podle zeměpisné šířky automaticky orientuje k jihu nebo severu, aby byl obzor dole.
- Hvězdný katalog se při přiblížení rozšiřuje z 5 až na 7 mag.
- Při větším zoomu se zobrazují DSO z OpenNGC do 11 mag jako jednoduché symboly bez názvů.
- Komety do odhadované 17 mag se počítají pro zvolený čas z denně obnovovaných elementů MPC a v lokální mapě se filtrují podle obzoru.

### Změněno

- Místní archiv Sky Coverage uchovává posledních 92 dní vůči nejnovějším datům v balíčku.
- Aktualizace MPC nejprve ověří a rozbalí nový archiv do dočasného umístění a až poté atomicky nahradí místní data.
- Počet dní zobrazeného pokrytí se počítá od zvolené noci; chybějící novější MPC data se již nenahrazují starším obdobím.
- Spojnice a hranice souhvězdí i kružnice výšek nad obzorem používají v MER správnou stupňovou škálu a nejsou už chybně odfiltrovány.
- MER po dokončení posunu nebo zoomu znovu vystředí projekci na aktuální výřez; souhvězdí a ostatní vrstvy tak nejsou deformované podle původního pevného středu mapy.
- Aktualizace Sky Coverage stahuje archiv z aktuální oficiální adresy MPC `/iau/data/skycov.tgz`.
- Plánování používá pravý boční panel, zatímco mapa zůstává vysoká; hlavní ovládání využívá volné místo v hlavičce.
- Archiv vlastních polí má výrazné barvy odlišné od tyrkysových polí MPC Sky Coverage.
- Tlačítka soumraku přepnou mapu do lokálního MER pohledu a vystředí ji nad azimutem Slunce.
- Stínování pod obzorem se určuje inverzní projekcí celého výřezu, takže se ani po posunu nepřevrátí; obzor má tmavě hnědou barvu.
- Aktualizace MPC průběžně ukazuje stažená procenta a po dokončení výslovně oznámí nejnovější datum v balíčku.
- Číselná pole času, počtu dnů a mezní magnitudy mají čitelnější kontrast, dostatek místa pro hodnotu a jednoznačné šipky nahoru/dolů; mění čas po minutě, období po dni a magnitudu po 0,5 mag.
- Plánovač nabízí mozaiku vodorovnou s místním obzorem; u EQ sestavy zachová natočení kamery a automaticky zmenší rozteč podle paralaktického náklonu, aby zůstal alespoň zadaný minimální překryv.
- Ukládání plánu navrhne název ve tvaru `datum_noci_sestava_Evening/Morning` a současně vytvoří samostatný textový soubor se souřadnicemi polí a pozorovacím skriptem.

V tomto souboru jsou zaznamenány významné změny projektu SkyCoverage.

## 2026-10-04

### Přidáno

- Volitelná vrstva aktuálních poloh Merkuru, Venuše, Marsu, Jupiteru, Saturnu, Uranu a Neptunu s barevnými značkami a krátkými popisky.
- V lokálním pohledu se planety zobrazují jen tehdy, když jsou nad obzorem.

## 2026-10-03

### Přidáno

- Stereografická projekce pro úhlově věrnější tvary souhvězdí; v lokálním pohledu má zenit uprostřed a obzor jako kružnici.
- Volitelná Bayerova označení jasných hvězd řeckými písmeny. Při přiblížení se zobrazují hvězdy až do 4 mag.
- Barevné odstíny hvězd podle indexu B−V a jemná záře nejjasnějších hvězd.

### Změněno

- Grafické rozhraní používá moderní tmavé panely, sjednocené ovládací prvky, výraznější stavové informace a samostatný souřadnicový panel.
- Mapa zachovává stejné měřítko obou os, takže se geometrie projekcí neroztahuje podle poměru stran okna.
- Spojnice souhvězdí jsou zelené a jejich hranice tlumeně jantarové pro snazší rozlišení.

## 2026-09-26

### Přidáno

- Volitelná vrstva hranic 88 souhvězdí v globálním i lokálním zobrazení.
- Přepínač „Hranice souhvězdí“, kterým lze novou vrstvu kdykoliv skrýt.
- Lambertova válcová plochojevná projekce jako třetí možnost zobrazení mapy.

### Opraveno

- Ranní soumrak vybírá vzestup Slunce přes zadanou výšku, takže u amerických observatoří nezamění doznívající večerní soumrak krátce po 00:00 UTC za ranní.
- Časové pole i tlačítka pro posun času zachovávají časové pásmo UTC a nepřebírají místní pásmo Windows.

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
- Čára obzoru je silnější, sytě hnědá a vykresluje se nad ostatními čarami oblohy.
- Kreutzovy koridory sahají až k bodu 5 dní před přísluním.
- Střední dráhy obou Kreutzových populací mají značky pro 5, 10, 20, 30 a 40 dní před přísluním.
- Ovladač „Dny zpět“ je omezený na uchovávaných 62 dní a tlačítko se jmenuje „Načíst data“.
