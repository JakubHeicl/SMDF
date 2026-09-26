# Simple MD Framework

Malý modulární framework pro klasickou molekulární dynamiku v Pythonu. Projekt
odděluje stav systému, interakce, integrátory, termostaty, konfiguraci a výstup,
aby bylo možné jednotlivé části měnit bez zásahu do hlavní simulační smyčky.

Framework je určen pro výukové a validační úlohy.
Referenční scénáře v adresáři `tests/` obsahují stejné počáteční stavy pro
framework a LAMMPS.

## Funkce

- konfigurace simulace pomocí TOML,
- načtení souřadnic z XYZ nebo vytvoření kubické mřížky,
- ortorombický periodický box a minimum-image convention,
- Lennardův–Jonesův potenciál s volitelným cutoffem a posunem energie,
- harmonický potenciál,
- skládání více potenciálů pomocí `SumInteraction`,
- integrátory velocity Verlet a BAOAB,
- Andersenův a Langevinův termostat,
- Numba kompilace párové Lennardovy–Jonesovy smyčky,
- příkazy CLI pro kontrolu a spuštění konfigurace,
- pět referenčních scénářů s výstupy frameworku a LAMMPS.

## Požadavky a instalace

Projekt vyžaduje Python 3.11 nebo novější. Běhové závislosti jsou uvedeny
v `requirements.txt`:

- NumPy,
- Numba.

Příklad instalace ve Windows PowerShell:

```powershell
git clone <URL>
cd simple_md_framework
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

V Linuxu nebo macOS se prostředí aktivuje pomocí:

```bash
source .venv/bin/activate
```

První výpočet Lennardovy–Jonesovy interakce může být pomalejší, protože Numba
při prvním volání zkompiluje párovou smyčku.

## Rychlé spuštění

Konfiguraci lze nejprve zkontrolovat bez vytvoření výstupních souborů:

```powershell
python -m smdf check tests/2/framework.toml
```

Simulace se spustí příkazem:

```powershell
python -m smdf run tests/2/framework.toml --output results/demo
```

Argument `--output` přepíše výstupní cestu uvedenou v TOML. Bez tohoto
argumentu se cesta vyhodnocuje relativně ke konfiguračnímu souboru.

Program záměrně nepřepisuje existující výstupní adresář. Pro opakovaný běh je
proto nutné zvolit nový adresář nebo starý výstup odstranit.

## TOML konfigurace

Následující konfigurace vytvoří 125 částic na kubické mřížce a provede krátkou
NVE simulaci:

```toml
[simulation]
units = "reduced"
dt = 0.002
steps = 100
seed = 42

[initial]
source = "lattice"
temperature = 1.0
n_side = 5
density = 0.5
label = "Ar"

[system]
mass = 1.0

[interaction]
type = "lennard_jones"
epsilon = 1.0
sigma = 1.0
cutoff = 2.5
shift_energy = true

[integrator]
type = "velocity_verlet"

[thermostat]
type = "none"

[output]
directory = "results/framework"
thermo_stride = 1
trajectory_stride = 10
```

### Sekce konfigurace

| Sekce | Význam |
|---|---|
| `simulation` | Jednotky, časový krok, počet kroků a seed generátoru náhodných čísel. |
| `initial` | Zdroj souřadnic a počáteční teplota. |
| `system` | Hmotnost částic a případně rozměry boxu. |
| `interaction` | Typ a parametry potenciálu. |
| `integrator` | Algoritmus propagace. |
| `thermostat` | Regulace teploty. |
| `output` | Výstupní adresář a intervaly zápisu. |

`seed` řídí počáteční náhodné hybnosti i stochastické termostaty. Stejný
vstup a seed proto vytvoří stejný běh frameworku.

### Vstup z mřížky

Pro `source = "lattice"` vznikne `n_side³` částic. Délka kubického boxu se
spočítá z počtu částic a hustoty. V této variantě se `system.box` nezadává.

```toml
[initial]
source = "lattice"
temperature = 1.0
n_side = 5
density = 0.5
label = "Ar"
```

### Vstup z XYZ

```toml
[initial]
source = "xyz"
path = "particles.xyz"
temperature = 0.0

[system]
mass = 1.0
box = [10.0, 10.0, 10.0]
```

Cesta k XYZ je relativní k TOML souboru. Framework načítá první snímek běžného
XYZ formátu:

```text
2
optional comment
Ar 1.0 2.0 3.0
Ar 4.0 5.0 6.0
```

XYZ neobsahuje parametry potenciálu ani standardizovaný periodický box, a proto
tyto údaje zůstávají v TOML. Pokud se `system.box` vynechá, systém je
neperiodický.

### Interakce

Dostupné typy:

```toml
[interaction]
type = "lennard_jones"
epsilon = 1.0
sigma = 1.0
cutoff = 2.5
shift_energy = true
```

```toml
[interaction]
type = "harmonic"
k = 1.0
```

Harmonický potenciál je definován vzhledem k počátku souřadnic a vyžaduje
neperiodický XYZ vstup. Pro periodický LJ systém musí být cutoff menší než
polovina nejkratší strany boxu.

Více příspěvků lze sečíst:

```toml
[interaction]
type = "sum"

[[interaction.terms]]
type = "lennard_jones"
epsilon = 1.0
sigma = 1.0
cutoff = 2.5
shift_energy = true

[[interaction.terms]]
type = "harmonic"
k = 0.1
```

Tato konkrétní kombinace je v současné implementaci použitelná pouze pro
neperiodický XYZ vstup, protože `HarmonicPotential` periodický box odmítá.

### Integrátory a termostaty

| Konfigurační typ | Implementace |
|---|---|
| `velocity_verlet` | Velocity Verlet, termostat se případně aplikuje na konci kroku. |
| `baoab` | B/2–A/2–O–A/2–B/2 splitting, převážně určený pro Langevinův termostat|
| `none` | Bez termostatu. |
| `andersen` | Náhodné obnovení hybností s pravděpodobností `dt / tau`. |
| `langevin` | Langevinův termostat s `gamma = 2 / tau`. |

Příklad Langevinovy konfigurace:

```toml
[integrator]
type = "baoab"

[thermostat]
type = "langevin"
temperature = 1.0
tau = 1.0
```

## Jednotky

Framework podporuje pouze redukované Lennardovy–Jonesovy jednotky:

```toml
[simulation]
units = "reduced"
```

Interně platí

```text
sigma = epsilon = m_reference = k_B = 1.
```

Délka je tedy vyjádřena v násobcích `sigma`, energie v násobcích
`epsilon`, hmotnost v násobcích referenční hmotnosti a teplota jako
`k_B T / epsilon`. Vstupní hodnoty nejsou interpretovány jako metry, jouly
ani kelviny.

## Výstup simulace

Každý úspěšný běh vytvoří:

| Soubor | Obsah |
|---|---|
| `thermo.csv` | Krok, čas, kinetická, potenciální a celková energie a teplota. |
| `trajectory.xyz` | Vícesnímková trajektorie souřadnic. |
| `initial.xyz` | Počáteční souřadnice. |
| `final.xyz` | Konečné souřadnice. |
| `input.toml` | Kopie původní konfigurace. |
| `metadata.json` | Vyhodnocené nastavení, počet částic, box a počet stupňů volnosti. |

`thermo_stride` a `trajectory_stride` nastavují interval zápisu. Počáteční
a konečný stav se uloží vždy.

## Architektura

```text
TOML / XYZ
    │
    ▼
config.py ── vytvoření SystemState, Interaction, Thermostat a Integrator
    │
    ▼
Simulation.run()
    │
    ├── Integrator.step()
    │      ├── Interaction.compute()
    │      └── Thermostat.apply()
    │
    └── Output.write()
```

| Modul | Odpovědnost |
|---|---|
| `state.py` | Datová reprezentace částic: polohy, hybnosti, hmotnosti, značky a box. |
| `interactions.py` | Rozhraní interakcí, LJ, harmonický potenciál a součet interakcí. |
| `thermostat.py` | Rozhraní a implementace termostatů. |
| `integrator.py` | Velocity Verlet a BAOAB. |
| `simulation.py` | Čas, počítadlo kroků a hlavní simulační smyčka. |
| `io.py` | TOML/XYZ vstup a zápis výsledků. |
| `config.py` | Validace konfigurace a sestavení simulace. |
| `__main__.py` | CLI příkazy `check` a `run`. |

Interakce vrací dvojici `(forces, potential_energy)`. Integrátor vlastní
interakci a termostat a mění `SystemState` na místě. Díky společným
abstraktním rozhraním lze přidat novou implementaci bez změny
`Simulation.run()`.

## Referenční scénáře

Adresář `tests/` obsahuje vstupy i již vytvořené výstupy frameworku a LAMMPS.

| Test | Účel |
|---|---|
| `tests/1` | Dvě částice v minimu LJ potenciálu; analyticky `U = -1` a `F = 0`. |
| `tests/2` | Deterministická NVE simulace 125 částic porovnaná krok po kroku s LAMMPS. |
| `tests/3` | Interakce přes hranici periodického boxu a minimum-image convention. |
| `tests/4` | Langevinův a Andersenův termostat, LAMMPS poskytuje Langevinovu referenci. |
| `tests/5` | Harmonický oscilátor porovnaný s LAMMPS a řešením `x(t) = cos(t)`. |

V každé složce jsou:

- frameworkové vstupy `framework*.toml` a případné XYZ,
- frameworkové výstupy v `results/`,
- LAMMPS vstup `*.in` a `initial.data`,
- LAMMPS termodynamický log a atomový dump.

Při opakovaném běhu frameworku použijte nový výstupní adresář, protože
výsledky uložené v `tests/*/results` se automaticky nepřepisují.

## Známá omezení

- LJ výpočet používá přímou párovou smyčku s časovou složitostí
  `O(N²)`; není implementován seznam sousedů.
- Podporován je pouze ortorombický box.
- LJ parametry jsou společné pro všechny částice, nejsou implementovány
  parametry podle dvojic atomových typů ani směšovací pravidla.
- Nejsou implementovány vazby, úhly, constraints
- Framework neprovádí paralelní doménový rozklad.
- Periodické vzdálenosti používají minimum-image convention, ale souřadnice
  trajektorie se po propagaci explicitně nebalí zpět do primární buňky.
- Posunutý LJ potenciál je spojitý v cut-offu, jeho síla však spojitá není.

## Přidání nové komponenty

Nová interakce implementuje:

```python
class MyInteraction(Interaction):
    def compute(self, state):
        # return forces with shape (N, 3), scalar potential energy
        return forces, potential_energy
```

Nový termostat implementuje `apply(state)` a nový integrátor
`step(state)`. Aby byla komponenta dostupná z TOML, je následně potřeba
doplnit její validaci a konstrukci do `config.py`.
