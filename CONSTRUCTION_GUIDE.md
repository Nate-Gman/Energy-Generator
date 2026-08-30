# CryoLavaTunnel — Construction & Build Guide

<p align="center">
  <strong>Complete construction instructions for the CryoLavaTunnel energy harvester</strong><br>
  <em>From site selection through commissioning — every subsystem, every component, every step</em>
</p>

---

## Table of Contents

1. [Safety Warning & Prerequisites](#safety-warning--prerequisites)
2. [Site Selection & Geological Survey](#site-selection--geological-survey)
3. [Phase 1 — Site Preparation & Infrastructure](#phase-1--site-preparation--infrastructure)
4. [Phase 2 — Access Tunnel Construction](#phase-2--access-tunnel-construction)
5. [Phase 3 — Cold-Air Cavern Excavation](#phase-3--cold-air-cavern-excavation)
6. [Phase 4 — Cavern Lining & Ultra Thermal Insulation](#phase-4--cavern-lining--ultra-thermal-insulation)
7. [Phase 5 — Tunnel Bore Construction](#phase-5--tunnel-bore-construction)
8. [Phase 6 — Lava Heat Exchanger Installation](#phase-6--lava-heat-exchanger-installation)
9. [Phase 7 — Turbine & Generator Installation](#phase-7--turbine--generator-installation)
10. [Phase 8 — Bottoming Cycle Construction](#phase-8--bottoming-cycle-construction)
11. [Phase 9 — Exit Fans, Stack & Nozzle](#phase-9--exit-fans-stack--nozzle)
12. [Phase 10 — Electrical & Switchyard](#phase-10--electrical--switchyard)
13. [Phase 11 — Control Systems & SCADA](#phase-11--control-systems--scada)
14. [Phase 12 — Auxiliary Systems](#phase-12--auxiliary-systems)
15. [Phase 13 — Commissioning & Testing](#phase-13--commissioning--testing)
16. [Phase 14 — Operation & Maintenance](#phase-14--operation--maintenance)
17. [BOM Cross-Reference Index](#bom-cross-reference-index)
18. [Critical Engineering Challenges](#critical-engineering-challenges)
19. [Permitting & Regulatory](#permitting--regulatory)

---

## Safety Warning & Prerequisites

### WARNING

**This is a conceptual design study, NOT a construction blueprint.** No such
system has ever been built. The 3000°C lava temperatures, 300 bar cryogenic
storage, and 110 TW power levels are far beyond any existing engineering
practice. Building this system would require:

- **Materials that do not exist** (turbine blades rated for 3000°C)
- **Construction methods never attempted** (tunneling within 200m of active lava)
- **Scale never achieved** (6 km³ underground cavern — 1/6 the volume of Lake Mead)
- **Thermal management never demonstrated** (keeping -150°C next to 3000°C)

**Do NOT attempt to build any part of this system without consulting:
- Geologists (volcanology, geothermal, rock mechanics)
- Structural engineers (underground, high-pressure, cryogenic)
- Thermodynamicists (heat transfer, refrigeration, turbine design)
- Electrical engineers (high-voltage, grid connection, protection)
- Safety engineers (HAZOP, SIL, emergency systems)
- Environmental scientists (lava interaction, gas emissions, seismic)**

### Prerequisites

Before any construction:

1. **Geological feasibility study** — confirm lava body location, temperature,
   depth, composition, and stability
2. **Geotechnical investigation** — rock quality, stress regime, groundwater,
   thermal profile
3. **Environmental impact assessment** — emissions, seismic risk, water use,
   ecological impact
4. **Permitting** — mining, environmental, electrical, nuclear-equivalent safety
5. **Financing** — estimated CAPEX is in the trillions of USD
6. **Technology development** — materials, turbomachinery, insulation systems
   that do not currently exist at the required scale and temperature

---

## Site Selection & Geological Survey

### Requirements

The site must have:

| Requirement | Specification | Why |
|-------------|---------------|-----|
| Active lava body | 2000-3000°C, stable for 50+ years | Heat source |
| Lava depth | 500-2000 m below surface | Accessible by tunneling |
| Lava volume | >1 km³ | Sustained heat supply |
| Rock quality | Hard rock (granite, basalt, gneiss) | Cavern structural integrity |
| Groundwater | Minimal or controllable | Cavern sealing |
| Surface access | Road, rail, power grid nearby | Logistics |
| Seismic stability | Low to moderate seismic activity | Structural safety |
| Climate | Cool climate preferred | Aids cavern cooling |

### Survey Steps

1. **Satellite thermal imaging** — Identify lava bodies via IR hotspot detection
2. **Seismic tomography** — Map lava body geometry, depth, and temperature
3. **Exploratory drilling** — Core samples at proposed cavern and tunnel locations
4. **Geothermal gradient measurement** — Confirm rock temperature vs depth
5. **Rock mechanics testing** — Uniaxial compressive strength, Young's modulus,
   Poisson's ratio, thermal conductivity, thermal expansion coefficient
6. **Stress field measurement** — In-situ stress via hydraulic fracturing tests
7. **Groundwater assessment** — Piezometric levels, flow rates, chemistry
8. **Seismic hazard analysis** — Historical seismicity, fault mapping,
   ground motion prediction
9. **Lava composition analysis** — Viscosity, temperature, gas content,
   crystallization behavior
10. **Topographic survey** — Surface mapping for site layout, roads,
    switchyard, buildings

### Equipment Needed

- Drilling rig (diamond core, 1000m+ depth capability)
- Seismograph network (temporary, 20+ stations)
- Thermal gradient probes
- Borehole camera and acoustic televiewer
- In-situ stress testing equipment
- Laboratory rock testing equipment
- UAV/drone for topographic survey
- Ground-penetrating radar

---

## Phase 1 — Site Preparation & Infrastructure

### 1.1 Access Roads

**BOM ref:** [71] Site Civil Works

```
Specification:
  - Width: 8 m (two-lane, heavy haul)
  - Length: 5-20 km from nearest public road
  - Load capacity: 100+ tonnes (for transformer, turbine transport)
  - Surface: Gravel with geotextile base, paved at switchyard
  - Drainage: Culverts every 200 m, ditches both sides
  - Grade: Maximum 8% for heavy haul routes
```

**Construction steps:**
1. Clear vegetation along surveyed alignment
2. Cut and fill to design grade
3. Install geotextile and compact subgrade
4. Lay gravel base course (300 mm)
5. Install drainage culverts and ditches
6. Compact wearing surface
7. Pave switchyard access and transformer pad areas

### 1.2 Site Fencing & Security

**BOM ref:** [71] Site Civil Works

```
Specification:
  - Perimeter fence: 3 m chain-link with barbed wire top
  - Gates: 2 (main vehicle, pedestrian)
  - Security cameras: 16 (coverage of all access points)
  - Lighting: LED every 50 m along perimeter
```

### 1.3 Temporary Construction Power

```
Specification:
  - Diesel generators: 2 x 2 MW (redundant)
  - Distribution: 480V, 415V, 240V
  - Construction water: 100 m³/day from wells or tanker
  - Construction air: 2 x 100 m³/min compressors
```

### 1.4 Construction Camp

```
Specification:
  - Offices: 20-person modular building
  - Locker rooms: 100-person capacity
  - Workshop: 500 m², crane-equipped
  - Warehouse: 1000 m², climate-controlled for sensitive equipment
  - First aid station: with trauma kit and oxygen
  - Helipad: for emergency evacuation
```

### 1.5 Meteorological Station

**BOM ref:** [73] Meteorological Station

```
Specification:
  - Anemometer (wind speed/direction)
  - Thermometer (dry bulb, wet bulb)
  - Barometer
  - Rain gauge
  - Solar pyranometer
  - Data logger with satellite uplink
```

---

## Phase 2 — Access Tunnel Construction

### 2.1 Access Tunnel

**BOM ref:** [18] Access Tunnel

```
Specification:
  - Diameter: 6 m
  - Length: 420 m (from surface to cavern)
  - Grade: 1:10 (descending toward cavern)
  - Lining: 350 mm precast concrete segments
  - Ventilation: 800 mm duct, forced ventilation
  - Lighting: LED every 20 m with battery backup
  - Drainage: 150 mm pipe, gravity flow to sump
```

**Construction steps:**
1. **Portal construction** — Excavate portal area, install ground support
   (shotcrete, rock bolts, steel arches), pour portal collar
2. **TBM launch** — Install TBM in portal, begin boring
3. **Segment installation** — Erect precast concrete segments behind TBM
   (350 mm thick, EPDM gaskets between segments [BOM 82])
4. **Segment bolting** — Torque bolts [BOM 81] to design specification
5. **Grouting** — Backfill annular gap with pea gravel and cement grout
6. **Ventilation duct installation** — Suspend 800 mm duct from crown
7. **Drainage installation** — Install 150 mm drainage pipe in invert
8. **Lighting installation** — LED fixtures every 20 m [BOM 86]
9. **Communication cable** — Fiber optic [BOM 87] for construction comms
10. **Breakthrough** — TBM breaks into cavern excavation area

**Equipment:**
- TBM (Tunnel Boring Machine) — 6 m diameter, hard rock
- Segment erector
- Grout pump
- Ventilation fan (forced, 100 m³/min)
- Locomotive for muck haulage
- Concrete segment casting yard (surface)

### 2.2 Access Tunnel Safety

**BOM ref:** [24] Emergency Ventilation & Refuges, [100] Escape Respirators

```
Specification:
  - Escape refuges: 1 per 200 m (2 total), pressurized, 24-hour capacity
  - Emergency ventilation: reversible fan, 200 m³/min
  - Escape respirators: 50 units at portal, 50 at cavern end
  - Fire suppression: water mist at portal and cavern end
  - Communication: hardwired phone + radio + fiber
```

---

## Phase 3 — Cold-Air Cavern Excavation

### 3.1 Cavern Excavation

**BOM ref:** [1] Cold-Air Cavern, [137] Roadheader, [139] Shotcrete Robot,
[140] Rock Bolter, [104] Rock Bolts

```
Specification (per system):
  - Volume: 6.0 x 10^9 m³ (6 km³)
  - Depth: 30 m below surface
  - Dimensions: ~1820 m x 1820 m x 1820 m (cubic equivalent)
  - Flat span: 45 m (largest clear span)
  - Crown height: 30 m
  - Floor area: 167,000 m²
  - Excavation method: drill-and-blast + roadheader
  - Lava proximity: 200 m from 3000°C lava body
  - Rock temperature at cavern: ~671°C (lava thermal halo)
```

**CRITICAL:** The cavern is 200 m from a 3000°C lava body. The surrounding
rock is ~671°C. Excavation will encounter extremely hot rock. All personnel
must work in cooled environments with heat-protective equipment.

**Construction steps:**

1. **Pilot drift** — Drive a 4m x 4m exploration drift around the cavern
   perimeter to confirm rock quality and temperature
2. **Rock bolt installation** — Install rock bolts [BOM 104] on 2m x 2m pattern
   as excavation progresses (C30/37 grade, stainless steel mesh [BOM 105])
3. **Top heading excavation** — Excavate the crown (top 8m) in 12m advances
   using roadheader [BOM 138] for soft rock, drill-and-blast for hard rock
4. **Initial support** — Apply 50mm shotcrete layer immediately after
   excavation using shotcrete robot [BOM 139]
5. **Rock bolt installation** — Install full-length rock bolts through
   shotcrete
6. **Bench excavation** — Remove bench material in 4m lifts, continuing
   rock bolting and shotcrete
7. **Final shotcrete** — Apply 600mm steel-fibre reinforced shotcrete
   (C30/37 + stainless steel mesh) [BOM 12, 105] in multiple layers
8. **Floor slab** — Pour reinforced concrete floor slab [BOM 79]
9. **Roof arches** — Install roof support arches [BOM 78] at 6m spacing
10. **Drainage channels** — Cut drainage channels [BOM 156] in floor
11. **Sump excavation** — Excavate drainage sump (2000 m³) [BOM 14, 157]
12. **Survey** — Laser scan the completed excavation to verify dimensions

**Equipment:**
- Roadheader (2 x 200kW cutting heads) [BOM 138]
- Drill jumbo (3-boom, for blast holes and rock bolts)
- Shotcrete robot (remote-controlled, 20 m³/hr) [BOM 139]
- Rock bolter (mechanized, 50 bolts/hr) [BOM 140]
- LHD (Load-Haul-Dump) loaders for muck removal
- Dump trucks (50-tonne, underground)
- Ventilation fan (high-capacity, 500 m³/min for hot rock)
- Cooling fans (spot cooling for personnel in 671°C rock)

### 3.2 Cavern Monitoring Installation (Pre-Lining)

**BOM ref:** [3] Cavern Monitoring, [88] Pressure Transmitters,
[89] Temperature Transmitters, [76] DTS Mapping, [77] DAS Monitor

```
Specification:
  - Pressure sensors: 12 (piezoresistive, 0-10 bar) [BOM 88]
  - Temperature sensors: 24 (RTD Pt100, -50 to 200°C) [BOM 89]
  - DTS fiber: 8.0 km (distributed temperature sensing) [BOM 76]
  - DAS fiber: 8.0 km (distributed acoustic sensing) [BOM 77]
  - Geophones: 16 (microseismic monitoring)
```

**Installation steps:**
1. Mount pressure sensors at 12 locations on cavern walls (4 per level:
   crown, mid, invert, sump)
2. Mount temperature sensors at 24 locations (8 per level: 3 levels)
3. Install DTS fiber optic cable in a serpentine pattern across all walls
   and crown, clipped to rock bolts every 2m
4. Install DAS fiber optic cable parallel to DTS, 1m offset
5. Install 16 geophones in boreholes (3m deep) around cavern perimeter
6. Route all cables to a junction box at the access tunnel portal
7. Connect to surface monitoring station via fiber optic [BOM 72]

---

## Phase 4 — Cavern Lining & Ultra Thermal Insulation

### 4.1 HDPE Waterproof Membrane

**BOM ref:** [12] Cavern Lining System, [13] Cavern Hydraulic Isolation Door

```
Specification:
  - Material: 8 mm welded HDPE membrane
  - Gas-tight: < 0.01% leak rate per day at 300 bar
  - Welding: hot-wedge fusion, 100% spark tested
  - Coverage: 100% of cavern walls, crown, and floor
```

**Installation steps:**
1. Clean shotcrete surface (remove loose material, dust)
2. Roll out HDPE sheets (3m wide rolls) starting from crown
3. Weld sheets using hot-wedge fusion welder
4. Spark test all welds (10kV, no pinholes accepted)
5. Install HDPE anchors through membrane into shotcrete
6. Seal around all penetrations (sensors, pipes, door)
7. Pressure test: pressurize cavern to 2 bar, monitor for 24 hours

### 4.2 Ultra Thermal Insulation (CRITICAL)

**BOM ref:** [2] Ultra Thermal Insulation System

```
Specification:
  - Total thickness: 500 mm
  - Layers (outside to inside):
    Layer 1: 50 mm aerogel blanket (k = 0.014 W/m·K)
    Layer 2: 100 mm vacuum insulated panels (VIP, k = 0.004 W/m·K)
    Layer 3: 30 layers multi-layer insulation (MLI, k = 0.00005 W/m·K)
    Layer 4: Reflective foil barriers between each layer
  - Combined R-value: 30 m²·K/W
  - Effective U_ground: 0.008 W/(m²·K)
  - Coverage: 100% of cavern exterior surface (1,000,000 m²)
```

**CRITICAL:** Without this insulation, the heat leak from the 671°C
surrounding rock would be ~246 MW — far exceeding the chiller capacity.
With ultra insulation, the leak drops to ~6.4 MW.

**Installation steps:**

1. **Surface preparation** — Clean HDPE membrane surface, repair any
   defects, ensure dry conditions

2. **Layer 1 — Aerogel blanket (50 mm)**
   - Cut aerogel blankets to size (2m x 1m sheets)
   - Adhere to HDPE using high-temperature silicone adhesive
   - Overlap joints by 100mm
   - Seal all joints with aerogel tape
   - Wear PPE: aerogel dust is a respiratory hazard

3. **Layer 2 — Vacuum insulated panels (100 mm)**
   - Cut VIP panels to size on surface (panels cannot be cut underground)
   - Transport pre-cut panels to cavern
   - Adhere to aerogel layer using thermal break tape
   - Butt joints with 10mm gap filled with aerogel tape
   - **DO NOT puncture VIP panels** — loss of vacuum destroys insulation
   - Install protective board over VIPs to prevent damage

4. **Layer 3 — Multi-layer insulation (30 layers)**
   - Install MLI blankets (aluminized Mylar + spacer material)
   - Layer 30 sheets in alternating orientation
   - Secure with non-conductive standoffs every 500mm
   - Seal perimeter with reflective tape
   - Install in 4m x 2m modules for handling

5. **Layer 4 — Reflective foil barriers**
   - Install reflective foil between each major layer
   - Overlap joints by 150mm
   - Tape all seams with aluminum foil tape

6. **Quality control**
   - Thermal imaging: scan entire surface for cold spots
   - Thickness verification: measure at 100 locations per 1000 m²
   - R-value verification: hot-plate test on sample panels
   - Vacuum integrity: acoustic test on all VIP panels

### 4.3 PU Foam Insulation (Inner Layer)

**BOM ref:** [12] Cavern Lining System

```
Specification:
  - Thickness: 200 mm closed-cell polyurethane foam
  - k-value: 0.024 W/m·K
  - Density: 40 kg/m³
  - Application: spray-applied on warm side (interior)
```

**Installation:**
1. Spray PU foam onto interior surface of ultra insulation
2. Apply in 50mm lifts (4 passes for 200mm total)
3. Allow curing between lifts (2 hours minimum)
4. Trim and smooth surface
5. Apply UV-protective coating (prevents degradation)

### 4.4 Hydraulic Isolation Door

**BOM ref:** [13] Cavern Hydraulic Isolation Door, [45] Pressure Relief Valve

```
Specification:
  - Diameter: 4000 mm
  - Actuation: hydraulic, fail-close
  - Pressure rating: 300 bar (matching cavern charge pressure)
  - Seal: double O-ring with monitoring port
  - Closure time: < 30 seconds
```

**Installation:**
1. Cast door frame into access tunnel lining at cavern entrance
2. Install hydraulic door in frame
3. Connect hydraulic power unit (redundant pumps)
4. Install pressure relief valve [BOM 45] bypassing door
5. Test: 10 open/close cycles, leak test at 300 bar

### 4.5 Active Refrigeration Plant

**BOM ref:** [4] Active Refrigeration Plant, [15] Cascade Refrigeration,
[16] Lava-Heated Absorption Chiller

```
Specification:
  - Cooling capacity: 1,000,000 kW_thermal (1 GW_thermal)
  - COP: 1.0 (base chiller), 0.3 (cascade cryogenic)
  - Power draw: 1,000,000 kW_elec (base), supplemented by absorption
  - Cascade cooling: multi-stage N2/He for -150°C
  - Absorption chiller: LiBr/H2O, COP 0.3, lava-heated, 85% of cooling load
```

**Installation steps:**

1. **Absorption chiller [BOM 16]**
   - Pour foundation pad at surface (near lava heat pipe exit)
   - Install LiBr/H2O absorption chiller unit
   - Connect lava heat pipe supply (hot side)
   - Connect cooling water return (cold side)
   - Connect chilled water supply to cavern cooling loop
   - Test: verify COP and capacity at design conditions

2. **Cascade refrigeration [BOM 15]**
   - Install multi-stage cascade refrigeration unit (N2 pre-cooling,
     He final stage)
   - Connect to cavern cooling loop
   - Install refrigerant storage and recovery system
   - Test: achieve -150°C at cavern design heat load

3. **Cooling loop piping**
   - Install insulated chilled water/brine piping from surface plant
     to cavern (through access tunnel)
   - Install circulation pumps (redundant, 3 x 50% capacity)
   - Install expansion tank and make-up system
   - Test: circulation at design flow, verify temperature drop

---

## Phase 5 — Tunnel Bore Construction

### 5.1 Main Tunnel Bores

**BOM ref:** [5] Tunnel Bore, [19] Parallel Tunnel Bores, [20] Tunnel Casing,
[137] TBM

```
Specification (per system):
  - Number of parallel bores: 48
  - Length: 7000 m (4.35 miles) per bore
  - Diameter: 20 m
  - Cross-section: 314.16 m² per bore
  - Total tunnel length: 336,000 m (48 x 7000m)
  - Lava contact length: 6000 m per bore
  - Bore method: TBM (hard rock) + drill-and-blast (lava contact zone)
  - Lining: 350 mm precast concrete segments
  - Casing: 4500 mm OD steel casing (non-contact zones)
  - Refractory: 120 mm alumina-silica firebrick (lava zone)
```

**Construction steps:**

1. **TBM assembly** — Assemble 20m-diameter TBM at portal [BOM 137]
   (This is a custom machine — no 20m TBM exists commercially)
2. **Launch chamber** — Excavate 25m x 25m x 50m launch chamber
3. **TBM launch** — Begin boring first bore
4. **Segment installation** — Erect 350mm precast segments [BOM 80]
   behind TBM, bolt segments [BOM 81], seal with EPDM gaskets [BOM 82]
5. **Backfill grouting** — Fill annular gap with cement grout
6. **Steel casing** — Install 4500mm OD steel casing [BOM 20] in
   non-contact zones (API L80 grade, Inconel 625 cladding in lava zone)
7. **Refractory lining** — Install 120mm alumina-silica firebrick [BOM 21]
   in lava contact zone using refractory anchors [BOM 83]
8. **Expansion joints** — Install expansion joints [BOM 22, 205-207]
   every 50m (140 total per bore) with 450mm stroke capacity
9. **Drainage** — Install 150mm drainage pipe [BOM 23, 157] in invert
10. **Lighting** — Install LED lighting every 20m [BOM 86]
11. **Repeat** — Repeat for all 48 parallel bores

**Equipment:**
- 20m diameter TBM (custom design, does not exist commercially)
- Segment transport cars
- Grouting equipment
- Refractory installation robot
- Crane for casing installation
- Ventilation fans (high-capacity for hot zone)

### 5.2 Tunnel Monitoring & Safety

**BOM ref:** [10] Monitoring & Safety System, [24] Emergency Ventilation & Refuges,
[38] Tunnel Monitoring, [39] Site Monitoring & Safety, [40] Dual System Interconnect,
[85] Access Platforms, [90] Flow Transmitters, [91] Vibration Transmitters,
[92] Thermocouple Wells

```
Specification:
  - Temperature sensors: 90 (every 20m along tunnel, RTD)
  - Pressure taps: 45 (every 40m)
  - Flow transmitters: 48 (one per bore) [BOM 90]
  - Vibration transmitters: 12 per turbine module [BOM 91]
  - Lava thermocouple wells: 6 (into lava contact zone) [BOM 92]
  - Access platforms: at each turbine module [BOM 85]
  - Dual system interconnect: cross-tie valves between systems [BOM 40]
  - Site monitoring: seismic, gas, weather [BOM 39]
```

**Installation:**
1. Mount temperature sensors on tunnel wall at 20m intervals
2. Install pressure taps with tubing to transmitter panels
3. Install flow transmitters at bore entrances
4. Drill 6 thermocouple wells 5m into lava contact zone
5. Install access platforms [BOM 85] at each turbine module
6. Install dual system interconnect valves [BOM 40] between parallel systems
7. Route all cables to junction boxes at 500m intervals
8. Connect to SCADA via fiber optic
9. Install site monitoring network [BOM 39]: seismic, gas, weather stations

### 5.3 Tunnel Expansion Joints & Drainage

**BOM ref:** [22] Thermal Expansion Joints, [23] Tunnel Condensate Drainage,
[84] Expansion Joint Bellows, [205] EJ Tie Rods, [206] EJ Internal Sleeve,
[207] EJ External Shroud

```
Specification:
  - Expansion joints: 140 per bore @ 50m spacing [BOM 22]
  - Stroke: 450mm per joint (handles 11m thermal growth over 1800m)
  - Bellows: Inconel 625, multi-ply [BOM 84]
  - Tie rods: limit axial movement [BOM 205]
  - Internal sleeve: protects bellows from flow [BOM 206]
  - External shroud: protects bellows from environment [BOM 207]
  - Drainage pipe: 150mm, full tunnel length [BOM 23]
```

**Installation:**
1. Install expansion joint assemblies at 50m intervals
2. Install internal sleeve [BOM 206] in flow direction
3. Install external shroud [BOM 207] over bellows
4. Adjust tie rods [BOM 205] to design cold-set position
5. Install 150mm drainage pipe [BOM 23] in tunnel invert
6. Connect drainage to condensate trays [BOM 158] and sumps

### 5.4 Tunnel Safety Systems

**BOM ref:** [24] Emergency Ventilation & Refuges, [101] Tunnel Fire Suppression

```
Specification:
  - Escape refuges: 6 per bore (pressurized, 24-hour capacity)
  - Emergency ventilation: reversible, 500 m³/min
  - Fire suppression: water mist deluge at turbine locations
  - Communication: fiber + radio + hardwired phone
  - Escape respirators: 20 per refuge [BOM 100]
```

---

## Phase 6 — Lava Heat Exchanger Installation

### 6.1 Shell-and-Tube Heat Exchanger

**BOM ref:** [6] Lava Heat-Exchange Contact, [25] Shell-and-Tube Lava HX,
[26] Heat Pipes, [179-181] Heat Pipe Sections

```
Specification (per system):
  - HX tubes: 200,000 (Inconel 617 or Haynes 230)
  - Tube OD: 25 mm
  - Tube length: 1500 m
  - U-value: 3500 W/(m²·K)
  - Fin enhancement: 30x effective area
  - Total UA: 1.71 TW/K
  - Heat pipes: sodium-charged, 50mm OD
```

**Installation steps:**

1. **Tube fabrication** — Manufacture 200,000 Inconel 617 tubes
   (25mm OD, 1500m length each — this is an unprecedented manufacturing run)
2. **Tube installation** — Feed tubes into pre-drilled boreholes from
   the tunnel wall into the lava body
3. **Tube sheet installation** — Install tube sheets at tunnel wall
   penetration points
4. **Header installation** — Install inlet and outlet headers connecting
   tube bundles
5. **Fin installation** — Install external fins on tubes (30x area enhancement)
6. **Pressure test** — Test each tube bundle at 1.5x design pressure
7. **Leak test** — Helium leak test all tube-to-tube-sheet welds

### 6.2 Heat Pipes

**BOM ref:** [26] Heat Pipes, [179] Evaporator, [180] Adiabatic,
[181] Condenser

```
Specification:
  - Type: sodium-charged, gravity-assisted thermosyphon
  - OD: 50 mm
  - Material: Inconel 617
  - Evaporator section: embedded in lava [BOM 179]
  - Adiabatic section: through refractory [BOM 180]
  - Condenser section: in tunnel air stream [BOM 181]
```

**Installation steps:**
1. Drill 50mm diameter boreholes from tunnel wall into lava body
   (depth: 10-50m into lava)
2. Insert heat pipe evaporator section into borehole
3. Pack borehole with thermally conductive grout
4. Route adiabatic section through refractory lining
5. Install condenser section in tunnel air flow path
6. Install heat pipe support brackets
7. Test: verify heat pipe thermal performance

---

## Phase 7 — Turbine & Generator Installation

### 7.1 Turbine Modules

**BOM ref:** [7] Turbine Array, [27] Stator Vanes, [28] Shaft & Bearings,
[106-115] Turbine Components, [102] Turbine Foundations

```
Specification (per system):
  - Stages: 28 axial-flow turbine stages
  - Rotor diameter: 3200 mm
  - Blades: 18 per stage, Inconel 718 single-crystal [BOM 106]
  - Stator vanes: 18 per stage [BOM 107]
  - Stage spacing: 4.5 m
  - RPM: 3600 (60 Hz)
  - Inlet T max: 650°C (material limit; model assumes 3000°C with
    ceramic coatings that do not exist)
  - Total turbine length: 126 m (28 x 4.5m)
```

**Installation steps:**

1. **Foundation preparation** [BOM 102]
   - Pour reinforced concrete pedestals at each turbine stage location
   - Install anchor bolts (M48, embedded 1m into concrete)
   - Grout base plates to design tolerance (±0.1mm)
   - Install vibration isolation pads

2. **Turbine casing installation** [BOM 108]
   - Lower split casing halves onto foundation
   - Bolt lower casing to base plate
   - Install diaphragms [BOM 109] in lower casing

3. **Rotor assembly** [BOM 110, 111]
   - Assemble rotor discs [BOM 110] on main shaft [BOM 111]
   - Install rotor blades [BOM 106] (18 per stage, single-crystal Inconel 718)
   - Dynamic balance the rotor assembly
   - Install labyrinth seals [BOM 115]

4. **Stator vane installation** [BOM 107]
   - Install stator vanes in casing diaphragms
   - Verify vane angles and clearances

5. **Bearing installation** [BOM 112, 113]
   - Install tilting-pad journal bearings [BOM 113] (2 per stage group)
   - Install thrust bearing [BOM 112]
   - Connect lube oil supply [BOM 57, 173, 174]

6. **Rotor installation**
   - Lower rotor assembly into lower casing
   - Install upper casing half
   - Torque casing bolts to specification
   - Check clearances (tip clearance, labyrinth clearance)

7. **Turning gear** [BOM 114]
   - Install turning gear for slow rotation during warm-up
   - Connect to turbine shaft

8. **Seal installation**
   - Install labyrinth + buffer air seals [BOM 115]
   - Connect seal air supply

### 7.2 Generators

**BOM ref:** [29] Generator, [116-122] Generator Components,
[58] H2 Cooling, [59] Seal Oil

```
Specification:
  - Rating: 45 MVA, 13.8 kV, 0.85 PF
  - Cooling: hydrogen-cooled, 75°C hot-spot limit
  - Type: 2-pole, direct-coupled, 3600 RPM
  - Efficiency: 96%
```

**Installation steps:**

1. **Stator installation** [BOM 116, 117]
   - Lower stator core [BOM 116] onto foundation
   - Install stator windings [BOM 117]
   - Terminal bushings [BOM 121]

2. **Rotor installation** [BOM 118]
   - Install generator rotor
   - Connect to turbine shaft (rigid coupling)
   - Align to within 0.05mm

3. **Exciter & AVR** [BOM 119, 120]
   - Install exciter [BOM 119] on outboard bearing
   - Install AVR panel [BOM 120]
   - Connect field wiring

4. **Isophase busduct** [BOM 122]
   - Install busduct from generator terminals to switchgear
   - Phase isolation, forced air cooling

5. **Hydrogen cooling system** [BOM 58, 176, 177]
   - Install hydrogen supply cylinders [BOM 176]
   - Install hydrogen control panel [BOM 177]
   - Install seal oil system [BOM 59, 175]
   - Purge with CO2 [BOM 178] before H2 fill
   - Pressurize with H2 to design pressure

6. **Generator protection** [BOM 49, 163]
   - Install generator protection relay [BOM 163]
   - Install generator circuit breaker [BOM 160]
   - Install surge arresters [BOM 161]
   - Install PT/CT instrument transformers [BOM 162]

### 7.3 Reheat Sections

**BOM ref:** [30] Reheat Sections

```
Specification:
  - Count: 48 reheat sections (between turbine stages)
  - Heat source: lava HX (same tubes as main heating)
  - Function: reheat air to T_hot between expansion stages
```

**Installation:**
1. Install reheat tube bundles between turbine stages
2. Connect to lava HX supply/return headers
3. Install reheat temperature sensors
4. Test: verify reheat achieves T_hot at design flow

### 7.4 Valves & Piping

**BOM ref:** [42] Main Isolation Valve, [43] Turbine Bypass Valve,
[44] Anti-Surge Valve, [45] Pressure Relief Valve, [46] Condensate Drain Valve,
[60] Cavern Recharge Piping, [61] Bottoming Piping, [62] Cooling Water,
[173] Lube Oil Console, [174] Lube Oil Cooler, [175] Seal Oil System

```
Valves:
  - Main isolation valve (cavern outlet): 300 bar rated [BOM 42]
  - Turbine bypass valve: routes air around turbine for startup [BOM 43]
  - Anti-surge valve: per turbine stage, prevents flow instability [BOM 44]
  - Pressure relief valve: cavern overpressure protection [BOM 45]
  - Condensate drain valve: removes water from tunnel [BOM 46]

Piping:
  - Cavern recharge piping: from compressor to cavern [BOM 60]
  - Bottoming cycle interconnecting piping: K to sCO2 to steam to ORC [BOM 61]
  - Cooling water system: to all heat exchangers and coolers [BOM 62]
```

**Installation:**
1. Install main isolation valve [BOM 42] at cavern outlet
2. Install turbine bypass valve [BOM 43] around each turbine module
3. Install anti-surge valves [BOM 44] at each turbine stage
4. Install pressure relief valve [BOM 45] on cavern
5. Install condensate drain valves [BOM 46] at low points in tunnel
6. Install cavern recharge piping [BOM 60] from compressor to cavern
7. Install bottoming cycle piping [BOM 61] between K, sCO2, steam, ORC
8. Install cooling water system [BOM 62] to all heat exchangers
9. Install lube oil console [BOM 173] and cooler [BOM 174]
10. Install seal oil system [BOM 175]
11. Pressure test all piping at 1.5x design pressure
12. Leak test with helium

### 7.5 Turbine Auxiliary Systems

**BOM ref:** [57] Lube Oil, [58] H2 Cooling, [59] Seal Oil,
[64] Instrument Air, [170-172] Air System

```
Lube oil system [BOM 57, 173, 174]:
  - Lube oil console with pumps (2x100%), cooler, filter
  - Oil piping to all bearings
  - Oil mist elimination

Seal oil system [BOM 59, 175]:
  - Seal oil console with pumps (2x100%)
  - Seal oil piping to generator seals
  - Oil-water separator [BOM 159]

Instrument air [BOM 64, 170-172]:
  - Compressor (2x100%) [BOM 170]
  - Air dryer [BOM 171]
  - Air receiver [BOM 172]
  - Distribution piping to all pneumatic valves
```

---

## Phase 8 — Bottoming Cycle Construction

### 8.1 Potassium Vapor Rankine Cycle

**BOM ref:** [31] Potassium Cycle, [123] Potassium Turbine,
[124] Potassium Condenser, [96] Potassium Pump

```
Specification:
  - Inlet T: 2000+°C
  - Condenser T: 800°C
  - Efficiency: 50%
  - Working fluid: potassium vapor
```

**Installation:**
1. Pour foundation for potassium turbine [BOM 123]
2. Install potassium turbine (specialty, high-temperature)
3. Install potassium condenser [BOM 124] (rejects heat to sCO2)
4. Install potassium feed pump [BOM 96]
5. Install potassium evaporator (receives heat from exhaust)
6. Connect piping [BOM 61] with high-temperature joints
7. Install potassium inventory and storage
8. Pressure test at 1.5x design pressure
9. Leak test with helium

### 8.2 Supercritical CO2 Brayton Cycle

**BOM ref:** [32] sCO2 Cycle, [125] sCO2 Turbine, [126] Recuperator,
[127] PHX, [128] Compressor, [93] Feed Pump

```
Specification:
  - Inlet T: 1000+°C
  - Compressor inlet T: 35°C (near-critical)
  - Efficiency: 48%
  - Working fluid: supercritical CO2 (73.8 bar)
```

**Installation:**
1. Pour foundation for sCO2 turbomachinery
2. Install sCO2 turbine [BOM 125]
3. Install sCO2 compressor [BOM 128]
4. Install recuperator [BOM 126]
5. Install primary heat exchanger [BOM 127] (receives heat from K condenser)
6. Install sCO2 feed pump [BOM 93]
7. Install sCO2 precooler (rejects heat to steam cycle)
8. Connect piping with supercritical-rated fittings
9. Pressure test at 1.5x design pressure (110 bar)
10. Leak test
11. Charge with CO2

### 8.3 Steam Rankine Cycle

**BOM ref:** [33] Steam Cycle, [129-131] Steam Turbines,
[132] Condenser, [133] Boiler, [94] Condensate Pump

```
Specification:
  - Inlet T: 500+°C
  - Condenser T: 30°C
  - Efficiency: 40%
  - Working fluid: water/steam
```

**Installation:**
1. Install HP steam turbine [BOM 129]
2. Install IP steam turbine [BOM 130]
3. Install LP steam turbine [BOM 131]
4. Install steam condenser [BOM 132] (rejects heat to ORC)
5. Install steam boiler/evaporator [BOM 133] (receives heat from sCO2 precooler)
6. Install condensate pump [BOM 94]
7. Install feedwater heaters and deaerator
8. Install steam piping with hangers and supports
9. Hydrostatic test at 1.5x design pressure
10. Boil-out and steam blow

### 8.4 Organic Rankine Cycle

**BOM ref:** [34] ORC Cycle, [134] Evaporator, [135] Turbine,
[136] Condenser, [95] Pump

```
Specification:
  - Working fluid: R245fa or silicone oil
  - Evaporator T: 120°C
  - Condenser T: 35°C
  - Efficiency: 12%
```

**Installation:**
1. Install ORC evaporator [BOM 134] (receives heat from steam condenser)
2. Install ORC turbine [BOM 135]
3. Install ORC condenser [BOM 136] (rejects heat to cooling tower)
4. Install ORC working fluid pump [BOM 95]
5. Charge with working fluid (R245fa)
6. Leak test

### 8.5 Cooling Tower

**BOM ref:** [70] Cooling Tower, [182-187] Tower Components,
[141] Demin Water, [142] CW Treatment, [187] CW Pumps

```
Specification:
  - Type: mechanical draft, counterflow
  - Fill media: PVC cross-flute [BOM 182]
  - Drift eliminators: [BOM 183]
  - Fans: [BOM 184]
  - Water distribution: [BOM 185]
  - Basin: [BOM 186]
  - Pumps: [BOM 187]
```

**Installation:**
1. Pour cooling tower foundation
2. Erect structural frame
3. Install basin [BOM 186] with waterproofing
4. Install fill media [BOM 182]
5. Install drift eliminators [BOM 183]
6. Install water distribution headers [BOM 185]
7. Install fan assemblies [BOM 184]
8. Install cooling water pumps [BOM 187]
9. Install water treatment system [BOM 142]
10. Install demineralized water plant [BOM 141]
11. Connect to ORC condenser and other cooling loads
12. Test: water circulation, fan operation, water quality

---

## Phase 9 — Exit Fans, Stack & Nozzle

### 9.1 Stack / Chimney

**BOM ref:** [9] Stack, [37] Stack Detail, [103] Stack Structure

```
Specification:
  - Height: 1200 m above tunnel exit
  - Diameter: 20 m
  - Function: buoyancy draft + exit jet acceleration
  - Structural support: reinforced concrete with steel liner
```

**Installation:**
1. Pour stack foundation (massive, 50m diameter x 10m deep)
2. Erect stack structure using slip-form concrete [BOM 103]
3. Install steel liner (Inconel, for high-temperature exhaust)
4. Install stack lighting (aviation obstruction lights)
5. Install lightning protection [BOM 66]
6. Install stack monitoring (temperature, flow, emissions)

### 9.2 Exit Nozzles

**BOM ref:** [35] Exit Nozzle & Jet

```
Specification:
  - Type: converging, fixed-geometry
  - Area: 2.0 m² per nozzle (48 nozzles per system)
  - Material: Inconel 617
  - Exit velocity: 959 m/s (Mach 0.85)
```

**Installation:**
1. Install nozzle assemblies at stack exit (48 per system)
2. Bolt to stack liner flanges
3. Install nozzle temperature monitoring
4. Test: verify flow distribution across nozzles

### 9.3 Exit Fans

**BOM ref:** [8] Exit Fans, [36] Fan Generators

```
Specification:
  - Type: ducted axial, generator-coupled
  - Diameter: 2800 mm
  - Blades: 8, carbon-fibre composite, 150°C rated
  - RPM: 1800
  - Generator: 850 kW PM direct-drive
  - Efficiency: 75-90%
  - Count: 48 per system
```

**Installation:**
1. Install fan housings in nozzle exit ducts
2. Install fan rotors (carbon-fibre blades)
3. Install PM direct-drive generators [BOM 36]
4. Connect generator output to switchgear
5. Install fan vibration monitoring [BOM 91]
6. Test: individual fan rotation, generator output

---

## Phase 10 — Electrical & Switchyard

### 10.1 Generator Switchgear (13.8 kV)

**BOM ref:** [11] Switchyard & Grid Connection, [41] Step-Up Transformer & Switchyard,
[51] 13.8kV Switchgear, [160] Generator Breaker,
[161] Surge Arresters, [162] PT/CT, [163] Gen Protection,
[164] Bus Differential

```
Specification:
  - Voltage: 13.8 kV
  - Configuration: metal-clad switchgear, indoor
  - Breakers: vacuum SF6, per turbine generator
  - Protection: differential, overcurrent, ground, reverse power
```

**Installation:**
1. Pour switchgear room floor (control room building [BOM 67])
2. Install switchgear cubicles
3. Install generator circuit breakers [BOM 160]
4. Install surge arresters [BOM 161]
5. Install PT/CT instrument transformers [BOM 162]
6. Install protection relays [BOM 163, 164]
7. Connect isophase busduct from generators
8. Wire CT/PT circuits to protection panels
9. Test: relay calibration, breaker timing, insulation

### 10.2 Step-Up Transformer

**BOM ref:** [200] Step-Up Transformer, [201] OLTC,
[202] Bushings, [203] Conservator, [204] Radiators

```
Specification:
  - Rating: matches total turbine MVA (28 x 45 = 1260 MVA per system)
  - Voltage: 13.8 kV / 132 kV
  - Cooling: ONAN/ONAF/ODAF
  - OLTC: on-load tap changer [BOM 201]
```

**Installation:**
1. Pour transformer foundation (oil containment basin, 110% oil volume)
2. Install fire wall between transformer and building
3. Transport transformer to site (special heavy haul, 200+ tonnes)
4. Position on foundation
5. Install bushings [BOM 202] (HV and LV)
6. Install conservator + Buchholz relay [BOM 203]
7. Install cooling radiators [BOM 204] and fans
8. Install OLTC [BOM 201] and tap changer control
9. Connect 13.8kV busduct from switchgear
10. Connect 132kV cable to switchyard
11. Oil filling and filtration
12. Test: insulation resistance, turns ratio, impedance, sweep frequency

### 10.3 132 kV Switchyard

**BOM ref:** [52] 132kV Switchyard, [188-194] Switchyard Components

```
Specification:
  - Voltage: 132 kV
  - Configuration: GIS (gas-insulated) or AIS (air-insulated)
  - Bays: 4 (2 transformer incomers, 2 line feeders)
  - Busbar: 2 main buses with transfer bus
```

**Installation:**
1. Pour switchyard foundation and cable trenches
2. Install SF6 circuit breakers [BOM 188]
3. Install disconnectors [BOM 189]
4. Install current transformers [BOM 190]
5. Install voltage transformers [BOM 191]
6. Install surge arresters [BOM 192]
7. Install post insulators [BOM 193]
8. Install busbar system [BOM 194]
9. Connect to step-up transformer and transmission lines
10. Install protection and control panels
11. Test: breaker timing, CT polarity, insulation, SF6 density

### 10.4 Cables

**BOM ref:** [195] MV Cable, [196] Control Cable, [197] Fiber Cable,
[198] Station Cable, [199] Cable Tray, [65] Cable Tray & Conduit

```
Specification:
  - MV cable: 13.8kV, XLPE insulated, per generator [BOM 195]
  - Control cable: multicore, 0.6kV rated [BOM 196]
  - Fiber cable: single-mode, SCADA network [BOM 197]
  - Station cable: 480V power, station service [BOM 198]
  - Cable tray: galvanized steel, ladder type [BOM 199]
```

**Installation:**
1. Install cable tray system [BOM 199] on walls and ceilings
2. Pull MV cables from generators to switchgear
3. Pull control cables from field instruments to control room
4. Pull fiber optic cables for SCADA network
5. Pull station service power cables
6. Terminate all cables
7. Test: insulation, continuity, fiber attenuation

### 10.5 Station Service & UPS

**BOM ref:** [53] Station Service Transformer, [54] UPS,
[55] DC Battery, [56] Diesel Backup

```
Specification:
  - Station service: 132kV/480V transformer [BOM 53]
  - UPS: 100 kVA, online double-conversion [BOM 54]
  - DC battery: 125V DC, 8-hour capacity [BOM 55]
  - Diesel backup: 2 MW, auto-start [BOM 56]
```

**Installation:**
1. Install station service transformer [BOM 53]
2. Install UPS system [BOM 54]
3. Install DC battery bank [BOM 55] and charger
4. Install diesel generator [BOM 56] with auto-transfer switch
5. Connect to all critical loads (SCADA, lighting, ventilation, seals)
6. Test: battery discharge, diesel start, UPS transfer

### 10.6 Grounding & Lightning Protection

**BOM ref:** [66] Grounding & Lightning

```
Specification:
  - Ground grid: copper conductor, 70 mm², buried 0.5m deep
  - Ground rods: copper-clad steel, 3m long, every 10m
  - Lightning mast: at switchyard, 30m height
  - Ground resistance: < 1 ohm
```

---

## Phase 11 — Control Systems & SCADA

### 11.1 SCADA System

**BOM ref:** [47] SCADA, [166] HMI, [168] Historian, [169] Network

```
Specification:
  - SCADA points: 2400 total
  - HMI workstations: 2 (operator + engineering) [BOM 166]
  - Historian server: 1, 10-year data retention [BOM 168]
  - Network: redundant industrial Ethernet [BOM 169]
  - Architecture: distributed, hot-standby
```

**Installation:**
1. Install SCADA server hardware in control room [BOM 67]
2. Install HMI workstations [BOM 166]
3. Install historian server [BOM 168]
4. Install network switches [BOM 169] (redundant ring topology)
5. Install fiber optic media converters
6. Configure SCADA database (all 2400 points)
7. Develop HMI screens (overview, cavern, tunnel, turbine, electrical)
8. Configure alarms and trends
9. Configure historian logging
10. Test: point-to-point verification, alarm response, trend display

### 11.2 PLC Controllers

**BOM ref:** [48] Governor, [167] PLC Controllers

```
Specification:
  - Turbine governor PLC: per turbine module [BOM 167]
  - Cavern monitoring PLC: for chiller, valves, door
  - Tunnel monitoring PLC: for ventilation, drainage, safety
  - I/O: remote I/O racks at each equipment location
```

**Installation:**
1. Install PLC cabinets at each equipment location
2. Wire all field instruments to PLC I/O
3. Program PLC logic (ladder/structured text)
4. Configure communication to SCADA
5. Test: I/O checkout, logic verification, communication

### 11.3 Protection & Control

**BOM ref:** [49] Generator Protection, [50] Sync Panel,
[165] Transformer Differential

```
Specification:
  - Generator protection: differential, overcurrent, loss of excitation,
    reverse power, out-of-step, stator earth fault [BOM 49]
  - Synchronizing panel: auto-sync to grid [BOM 50]
  - Transformer protection: differential, Buchholz, OLTC [BOM 165]
  - Bus protection: differential [BOM 164]
```

### 11.4 Control Room

**BOM ref:** [67] Control Room Building, [151] Precision AC

```
Specification:
  - Building: 200 m², climate-controlled [BOM 67]
  - HVAC: precision cooling [BOM 151], 22°C ± 2°C, 50% RH
  - Operator console: 2 positions, each with 3 monitors
  - Wall display: large-screen overview
  - Communication: phone, radio, satellite backup
  - Security: card access, CCTV
```

---

## Phase 12 — Auxiliary Systems

### 12.1 Buildings

**BOM ref:** [67] Control Room, [68] Turbine Hall, [69] Bottoming Building

```
Turbine hall [BOM 68]:
  - Size: 200m x 50m x 25m (for 28 turbine modules)
  - Crane: 50-tonne overhead [BOM 145-148]
  - HVAC: supply ducts [BOM 149], exhaust fans [BOM 150]
  - Lighting: high-bay LED

Bottoming cycle building [BOM 69]:
  - Size: 100m x 40m x 20m
  - Houses K, sCO2, steam, ORC equipment
  - Crane: 20-tonne overhead
```

### 12.2 Fire Protection

**BOM ref:** [63] Fire Protection, [152] Fire Pump, [153] Sprinklers,
[154] FACP, [155] Clean Agent

```
Specification:
  - Fire water pump skid: 2x100%, diesel + electric [BOM 152]
  - Sprinkler heads: throughout buildings [BOM 153]
  - Fire alarm control panel: addressable [BOM 154]
  - Clean agent: in control room (no water on electronics) [BOM 155]
  - Tunnel fire: water mist deluge [BOM 101]
```

### 12.3 Drainage & Condensate

**BOM ref:** [14] Drainage Sump, [97] Condensate Pump,
[156-159] Drainage Components

```
Specification:
  - Cavern drainage channels: [BOM 156]
  - Tunnel drainage pipe: 150mm [BOM 157]
  - Condensate collection trays: [BOM 158]
  - Oil-water separator: [BOM 159]
  - Condensate return pump: [BOM 97]
  - Sump capacity: 2000 m³ [BOM 14]
```

### 12.4 HVAC

**BOM ref:** [143] Turbine Hall HVAC, [144] Control Room HVAC,
[149-151] Ducts and AC

```
Specification:
  - Turbine hall: supply ducts [BOM 149], exhaust fans [BOM 150]
  - Control room: precision AC [BOM 151]
  - Bottoming building: general ventilation
  - Cavern: cooling via refrigeration plant (Phase 4)
```

### 12.5 Crane

**BOM ref:** [145-148] Crane Components

```
Specification:
  - Capacity: 50 tonnes
  - Type: overhead, double-girder
  - Span: 45m (turbine hall)
  - Components: bridge girder [BOM 145], hoist [BOM 146],
    end carriages [BOM 147], runway beam [BOM 148]
```

### 12.6 Gas Detection & Safety

**BOM ref:** [75] Gas Analysis, [99] ESD, [100] Respirators,
[101] Tunnel Fire

```
Specification:
  - Gas sensors: SO2, H2S, CO2, CO, O2, CH4 at exit [BOM 75]
  - ESD (Emergency Shutdown System): SIL-3 rated [BOM 98]
  - Escape respirators: at all refuges [BOM 100]
  - Gas detection: throughout tunnel and cavern
```

### 12.7 Vibration Monitoring

**BOM ref:** [74] Vibration Monitor, [91] Vibration Transmitters

```
Specification:
  - Accelerometers: 12 per turbine module [BOM 91]
  - Vibration monitoring system: [BOM 74]
  - Trip levels: per API 670
  - Analysis: FFT, trend, alarm
```

---

## Phase 13 — Commissioning & Testing

### 13.1 Pre-Commissioning (Cold)

1. **Loop checks** — Verify every instrument, every valve, every control
   signal from field to SCADA
2. **Insulation testing** — Megger all MV cables, transformers, generators
3. **Protection testing** — Inject test current, verify relay operation
4. **Valve testing** — Stroke all valves, verify limit switches
5. **Pump testing** — Run all pumps on test fluid
6. **Communication testing** — Verify all SCADA points, PLC comms
7. **Cavern pressure test** — Pressurize to 50% design, hold 24h, check leaks
8. **Cavern insulation verification** — Thermal imaging, verify no hot spots

### 13.2 Hot Commissioning

1. **Lava HX activation** — Slowly open HX to lava heat, monitor tube temps
2. **Refrigeration startup** — Start absorption chiller, then cascade
3. **Cavern cooldown** — Cool cavern to -150°C (this will take weeks)
4. **Cavern pressurization** — Compress air to 300 bar in stages
5. **Turbine first rotation** — Turn on turning gear, then slow roll
6. **Turbine speed ramp** — Ramp to rated speed (3600 RPM)
7. **Generator excitation** — Apply field, build voltage
8. **Synchronization** — Auto-sync to grid [BOM 50]
9. **Load ramp** — Ramp load at 2%/min [MONITOR_HW ramp limit]
10. **Full load test** — Operate at 100% valve for 24 hours
11. **Conservation audit** — Verify E_in = E_out + dU/dt (residual < 5%)
12. **Carnot audit** — Verify P_net ≤ η_Carnot × Q_lava
13. **Performance test** — Measure power, efficiency, EROI

### 13.3 Acceptance Tests

| Test | Method | Acceptance |
|------|--------|------------|
| Cavern leak rate | Pressure decay over 24h | < 0.1%/day |
| Insulation performance | Heat leak measurement | < 10 MW |
| Cavern temperature | RTD + DTS | -150°C ± 2°C |
| Turbine vibration | API 670 | < 2.8 mm/s RMS |
| Generator temperature | RTD + H2 cooler | < 75°C hot spot |
| Conservation residual | First Law audit | < 5% |
| Carnot audit | P_net vs η_Carnot × Q_lava | P_net < ceiling |
| EROI | E_out / E_recharge | > 1.0 |
| SCADA availability | Uptime during test | > 99.5% |
| Protection trip time | Injection test | < 80 ms |

---

## Phase 14 — Operation & Maintenance

### 14.1 Normal Operation

```
Startup sequence:
  1. Verify cavern T = -150°C, P = 300 bar
  2. Start lube oil pumps, seal oil pumps
  3. Start instrument air
  4. Open cavern isolation valve [BOM 42]
  5. Start turning gear, then roll turbine
  6. Ramp to rated speed (3600 RPM)
  7. Excite generator, build voltage
  8. Synchronize to grid [BOM 50]
  9. Ramp load at 2%/min
  10. Monitor all systems via SCADA

Shutdown sequence:
  1. Ramp load down at 2%/min
  2. Open generator breaker [BOM 160]
  3. Cooldown turbine (turning gear for 4h)
  4. Close cavern isolation valve [BOM 42]
  5. Stop lube oil after bearing temps < 50°C
  6. Maintain cavern cooling
```

### 14.2 Maintenance Schedule

| Component | Interval | Task |
|-----------|----------|------|
| Turbine bearings | 4000h | Oil sample, vibration check |
| Turbine blades | 8000h | Boroscope inspection |
| Generator | 4000h | H2 purity check, oil check |
| HX tubes | 6 months | Leak test, thickness check |
| Heat pipes | 12 months | Performance verification |
| Cavern insulation | 12 months | Thermal imaging |
| Cavern seal | 6 months | Leak test |
| Valves | 6 months | Stroke test, leak test |
| SCADA | 1 month | Backup, point verification |
| Protection relays | 6 months | Injection test |
| Battery bank | 1 month | Discharge test |
| Diesel generator | 1 month | Start test, load test |
| Fire system | 6 months | Full functional test |
| Gas detectors | 1 month | Calibration |
| Vibration monitors | 6 months | Calibration |

### 14.3 Recharge Cycle

```
Recharge sequence (after discharge):
  1. Close discharge valve, open recharge valve [BOM 42]
  2. Start 20-stage intercooled compressor [BOM 17]
  3. Compress atmospheric air to 300 bar (542 kJ/kg)
  4. Cool compressed air from 20°C to -150°C:
     a. Absorption chiller (85% of load, lava-powered) [BOM 16]
     b. Cascade refrigeration (15% of load, electric) [BOM 15]
  5. Fill cavern until P = 300 bar, T = -150°C
  6. Close recharge valve
  7. Cavern ready for next discharge
```

---

## BOM Cross-Reference Index

All 207 BOM assemblies mapped to construction phases:

| Phase | BOM Items | Description |
|-------|-----------|-------------|
| 1. Site Prep | 71, 73 | Civil works, met station |
| 2. Access Tunnel | 18, 24, 80-82, 86-87, 100 | Tunnel, safety, lining |
| 3. Cavern Excavation | 1, 3, 76-80, 88-92, 104, 137-140 | Excavation, monitoring |
| 4. Cavern Lining | 2, 12-16, 45, 78-79, 83, 105, 156 | Lining, insulation, cooling |
| 5. Tunnel Bore | 5, 19-24, 38, 83-87, 101, 137 | Bores, casing, refractory |
| 6. Lava HX | 6, 25-26, 92, 179-181 | HX tubes, heat pipes |
| 7. Turbines | 7, 27-30, 48, 57-59, 64, 102-115, 170-172 | Turbines, generators, aux |
| 8. Bottoming | 31-34, 61-62, 70, 93-96, 123-136, 141-142, 182-187 | K, sCO2, steam, ORC, CT |
| 9. Exit & Stack | 8-9, 35-37, 103 | Fans, nozzle, stack |
| 10. Electrical | 11, 29, 41-42, 49-56, 66, 160-199 | Switchgear, transformer, cables |
| 11. Control | 47, 67, 74-75, 98-99, 166-169 | SCADA, PLC, control room |
| 12. Auxiliary | 63, 68-69, 143-155, 157-159 | Buildings, fire, HVAC, crane |
| Valves | 42-46 | Isolation, bypass, relief, drain |
| Construction | 137-140 | TBM, roadheader, shotcrete robot |

---

## Critical Engineering Challenges

### 1. Materials (UNSOLVED)

| Challenge | Current State | Required |
|-----------|---------------|----------|
| Turbine blades at 3000°C | Inconel 718 rated to 700°C | 3000°C ceramic — does not exist |
| Lava containment | No material survives 3000°C for decades | Unknown |
| Cavern seal at -150°C + 300 bar | HDPE becomes brittle below -50°C | New seal material needed |
| HX tubes in lava | Inconel 617 rated to 1200°C | 3000°C service — does not exist |

### 2. Thermal Management (PARTIALLY SOLVED IN MODEL)

| Challenge | Model Solution | Real Status |
|-----------|----------------|-------------|
| Cavern at -150°C near 671°C rock | Ultra insulation (R=30) | Aerogel+VIP+MLI exists but not at this scale |
| Tunnel expansion (11m over 1800m) | 140 expansion joints per bore | Standard engineering |
| Lava HX at 3000°C | Inconel 617 + heat pipes | Material does not exist for 3000°C |

### 3. Scale (UNPRECEDENTED)

| Component | Required | Largest Existing |
|-----------|----------|------------------|
| Cavern volume | 6 km³ | ~0.005 km³ (existing CAES) |
| Tunnel diameter | 20 m | 15 m (Gotthard Base Tunnel) |
| Parallel bores | 48 | 1-2 (typical) |
| HX tubes | 200,000 | ~10,000 (refinery) |
| Turbine stages | 28 | 10-15 (gas turbine) |

### 4. Safety (EXTREME RISK)

- Tunneling within 200m of active lava — never attempted
- 300 bar cryogenic storage — never attempted at this scale
- 3000°C heat exchanger — never attempted
- 110 TW power handling — no grid can accept this

---

## Permitting & Regulatory

### Required Permits

1. **Mining permit** — for underground excavation
2. **Environmental impact assessment** — NEPA or equivalent
3. **Geothermal permit** — for lava heat extraction
4. **Air quality permit** — for exhaust emissions
5. **Water use permit** — for cooling tower make-up
6. **Electrical generation license** — for grid connection
7. **High-pressure vessel permit** — for 300 bar cavern
8. **Cryogenic facility permit** — for -150°C storage
9. **Seismic safety review** — for underground construction
10. **Volcanic risk assessment** — for lava proximity
11. **Emergency response plan** — for all hazards
12. **Decommissioning plan** — for end-of-life

### Regulatory Standards

- ASME B31.1 — Power Piping
- ASME B31.3 — Process Piping
- ASME Section VIII — Pressure Vessels
- API 670 — Machinery Protection
- API 617 — Centrifugal Compressors
- API 618 — Reciprocating Compressors
- IEEE C37 — Switchgear Standards
- IEEE 80 — Grounding
- NFPA 70 — National Electrical Code
- NFPA 850 — Fire Protection for Power Plants
- IEC 61511 — Functional Safety (SIL)
- ISO 13702 — Fire and Explosion Control

---

## Summary

This construction guide covers all 14 phases of building the CryoLavaTunnel
energy harvester, from site selection through commissioning and operation.
Every subsystem is addressed with specifications, construction steps,
equipment lists, and BOM cross-references.

**The physics is real. The engineering is conceptual. Many of the required
materials and construction methods do not exist. This guide is intended for
research and educational purposes — it is NOT a construction blueprint.**

**Do not build a tunnel over lava without consulting actual geologists,
structural engineers, and thermodynamicists.**
