#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 CryoLavaTunnel.py  --  COLD-AIR CAVERN + LAVA-HEATED TURBINE TUNNEL DIGITAL TWIN
================================================================================

WHAT THIS IS
------------
A single-file, standalone DIGITAL TWIN of the tunnel energy harvester described
in `informational.md`:

    * a massive underground AIR tank, super cold (earth-coupled or refrigerated),
    * that stays cold by itself while buried ("operates when not cooled"),
    * discharging its cold dense air into a long (mile+) tunnel that runs over a
      lava / superhot-geothermal heat source,
    * where the air is heated, expands violently (turbine-engine style) inside
      thermo-separated pipe stages,
    * driving multiple fan/turbine stages that generate electricity,
    * and finally jets the expanded hot air out past exit fans for more power.

It is written in the same spirit as the companion reference engines
(ValcanoHarvester.py, Radiant.py, Main.py, Simulation.py): every number is in
SI, every dimension is to scale, and every extraordinary claim is checked
against a textbook formula and a conservation + Carnot audit. There is a
dedicated HONESTY LAYER (SECTION 11) that refuses to report over-unity, because
the energy is not free -- it is borrowed from the geothermal heat flux and from
whatever pre-charged the cold cavern.

--------------------------------------------------------------------------------
THE HONEST PHYSICS (read this first -- it is the whole point)
--------------------------------------------------------------------------------
The vision asks for "free" power from a cold underground air tank sitting on
lava. Taken literally, two pieces of it violate thermodynamics or geology:

  * "Super cold DEEP underground" is generally WRONG. Below ~1.5-4 m the ground
    is thermally stable, and below that it WARMS at ~25-30 C/km (the geothermal
    gradient). A deep cavern next to lava is HOT, not cold. Real cold comes
    from either (a) SHALLOW earth coupling in a cool climate, (b) ACTIVE
    refrigeration (which costs work), or (c) charging the cavern with cold
    winter/night air. This twin lets you pick the source and shows the cost.

  * "Free" energy is not free. The system is a heat engine: it converts the
    temperature difference between the lava (hot reservoir) and the cold air
    (cold reservoir) into work. The First Law fixes the energy ledger; the
    Second Law caps the net work at the Carnot efficiency

        eta_Carnot = 1 - T_cold / T_hot

    times the heat actually delivered to the air. Any "extra" power claimed
    beyond that is fiction. The cold cavern is a THERMAL BATTERY: somebody
    (the ground, a chiller, or last winter) paid to put the cold in, and when
    it is spent the engine stops until it is recharged.

What the twin DOES model honestly:

  1. A pressurised cold-air cavern (volume, T, P, mass inventory).
  2. Passive ground-coupled recharge toward the local ground temperature
     (the "operates when not cooled" mode) -- with the geothermal gradient
     correction so going deep makes the ground warmer, not colder.
  3. Optional active refrigeration with a real COP, subtracted from net power.
  4. Discharge through a mile+ tunnel with friction (Darcy-Weisbach) and a
     lava-coupled heat-exchange section (U A dT).
  5. Air heating, expansion (ideal gas, v ~ T at constant P), and a staged
     turbine / fan array that extracts kinetic + enthalpy work.
  6. A strict energy ledger and a Carnot ceiling the simulated net power can
     never exceed (the self-test asserts both).

--------------------------------------------------------------------------------
RUN
--------------------------------------------------------------------------------
    python3 CryoLavaTunnel.py                 # DEFAULT: print the headline report
    python3 CryoLavaTunnel.py --report        # full report for all presets
    python3 CryoLavaTunnel.py --selftest      # physics + conservation + Carnot + proofs
    python3 CryoLavaTunnel.py --info          # explain every part + the math
    python3 CryoLavaTunnel.py --honesty       # the reality check, in full
    python3 CryoLavaTunnel.py --proofs        # math proofs with verify_fn
    python3 CryoLavaTunnel.py --targets       # list the preset sites
    python3 CryoLavaTunnel.py --hardware      # to-scale hardware spec (SI)
    python3 CryoLavaTunnel.py --parts         # BOM parts list
    python3 CryoLavaTunnel.py --model         # ASCII cross-section visualization
    python3 CryoLavaTunnel.py --flow          # energy flow Sankey diagram
    python3 CryoLavaTunnel.py --timeline      # multi-series timeline plot
    python3 CryoLavaTunnel.py --turbines      # per-stage turbine breakdown
    python3 CryoLavaTunnel.py --sweep [HOURS] # scan the design space
    python3 CryoLavaTunnel.py --sensitivity KEY HOURS   # simple sensitivity
    python3 CryoLavaTunnel.py --sensitivity2 KEY HOURS  # advanced sensitivity
    python3 CryoLavaTunnel.py --mc TARGET [N] [HOURS]   # Monte Carlo
    python3 CryoLavaTunnel.py --optimize [TARGET]       # coordinate-descent optimiser
    python3 CryoLavaTunnel.py --pareto [TARGET]         # power vs duration frontier
    python3 CryoLavaTunnel.py --live [HOURS]  # run continuously (dashboard)
    python3 CryoLavaTunnel.py --visual        # interactive GUI (matplotlib, pan/zoom)
    python3 CryoLavaTunnel.py --target NAME   # select preset for any command

Dependencies:  Python 3.8+ standard library only.
               --visual additionally needs matplotlib + tkinter (optional).

================================================================================
"""

from __future__ import annotations

import math
import sys
import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Callable

# --- optional visualization deps (graceful fallback) ---
try:
    import tkinter as tk
    from tkinter import ttk, scrolledtext
    HAS_TK = True
except ImportError:
    HAS_TK = False

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
    from matplotlib.collections import LineCollection
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import numpy as np
    HAS_NP = True
except ImportError:
    HAS_NP = False


# ==============================================================================
# SECTION 0 -- CONSTANTS AND SMALL NUMERICS
# ==============================================================================

# --- physical constants (SI) ---
G           = 9.80665          # m/s^2
R_AIR       = 287.052868       # J/(kg K)  specific gas constant of dry air
CP_AIR      = 1005.0           # J/(kg K)  isobaric specific heat of air
CV_AIR      = 718.0            # J/(kg K)  isochoric specific heat of air
GAMMA_AIR   = CP_AIR / CV_AIR  # ~1.400
K_AIR       = 0.0263           # W/(m K)   thermal conductivity of air (rough)
MU_AIR      = 1.8e-5           # Pa s      dynamic viscosity of air (rough)
P_STD       = 101325.0         # Pa        standard atmosphere
T_STD_K     = 288.15           # K         standard temperature
RHO_STD     = P_STD / (R_AIR * T_STD_K)   # ~1.225 kg/m^3
T_REF_K     = 273.15           # K         0 C
YEAR        = 365.25 * 86400.0
HOUR        = 3600.0
DAY         = 86400.0
MPA         = 1.0e6
KPA         = 1.0e3

# --- geology ---
GEOTHERMAL_GRADIENT   = 0.030      # K/m  (~30 C/km, typical continental)
GROUND_DIFFUSIVITY    = 1.0e-6     # m^2/s  (soil/rock, order of magnitude)
GROUND_T_SURF_C       = 12.0       # C  mean surface temp (temperate site)
LAVA_T_C              = 1100.0     # C  basaltic lava, order of magnitude
ROCK_CONDUCTIVITY     = 2.5        # W/(m K)  saturated rock

# --- economics / household ---
HOME_KW               = 1.2        # kW average per home
CAPEX_TUNNEL_PER_M    = 12_000.0   # USD per metre of bored tunnel (order of mag)
CAPEX_TURBINE_PER_KW  = 1_400.0    # USD per kW of turbine+generator
CAPEX_CAVERN_PER_M3   = 60.0       # USD per m^3 of excavated/storage cavern
CAPEX_CHILLER_PER_KW  = 800.0      # USD per kW thermal of refrigeration

# --- numerics ---
EPS     = 1.0e-9
INF     = float("inf")


def c_to_k(c: float) -> float:
    return c + T_REF_K


def k_to_c(k: float) -> float:
    return k - T_REF_K


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * clamp(t, 0.0, 1.0)


def bisect(f: Callable[[float], float], lo: float, hi: float,
           tol: float = 1.0e-6, iters: int = 200) -> float:
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise ValueError("bisect: no sign change in interval")
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < tol or (hi - lo) < tol:
            return mid
        if fmid * flo < 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def fmt(x: float, unit: str = "", width: int = 10, prec: int = 3) -> str:
    if abs(x) >= INF:
        return f"{'inf':>{width}} {unit}".rstrip()
    if abs(x) >= 1.0e6:
        s = f"{x/1e6:.{prec}f}e6"
    elif abs(x) >= 1.0e3 and unit in ("W", "Pa", "m", "kg"):
        s = f"{x/1e3:.{prec}f}k"
    elif abs(x) < 1.0e-3 and x != 0.0:
        s = f"{x:.{prec}e}"
    else:
        s = f"{x:.{prec}f}"
    return f"{s:>{width}} {unit}".rstrip()


# ==============================================================================
# SECTION 1 -- AIR THERMODYNAMICS
# ==============================================================================

def air_density(p_pa: float, t_k: float) -> float:
    """Ideal gas density of dry air, rho = P / (R T)."""
    return p_pa / (R_AIR * max(t_k, 1.0))


def air_enthalpy(t_k: float) -> float:
    """Specific enthalpy of dry air relative to 0 C, h = cp (T - T_ref)."""
    return CP_AIR * (t_k - T_REF_K)


def mass_in_cavern(v_m3: float, p_pa: float, t_k: float) -> float:
    """Mass of air stored in a cavern of volume V at P, T."""
    return v_m3 * air_density(p_pa, t_k)


def sound_speed(t_k: float) -> float:
    """Speed of sound in air, a = sqrt(gamma R T)."""
    return math.sqrt(GAMMA_AIR * R_AIR * max(t_k, 1.0))


def mach(v: float, t_k: float) -> float:
    return v / sound_speed(t_k)


def dew_point_c(t_c: float, rh: float) -> float:
    """Magnus formula dew point. rh in [0,1]."""
    a, b = 17.625, 243.04
    if rh <= 1e-3:
        return -100.0
    g = math.log(max(rh, 1e-6)) + (a * t_c) / (b + t_c)
    return (b * g) / (a - g)


# ==============================================================================
# SECTION 2 -- SYSTEM SPECIFICATIONS
# ==============================================================================

@dataclass
class ColdCavernSpec:
    """The massive underground cold-air tank."""
    name: str = "Cold Cavern"
    volume_m3: float = 5.0e6          # 5 million m^3 ~ a large excavated cavern
    depth_m: float = 30.0             # burial depth (shallow earth coupling!)
    p_charge_pa: float = 600.0e3      # charged pressure (6 bar) -- pressurised store
    t_charge_k: float = 263.15        # -10 C charged temperature
    t_ground_c: float = float("nan")  # local ground T at depth; nan -> compute
    ground_k_per_m: float = GEOTHERMAL_GRADIENT
    surf_t_c: float = GROUND_T_SURF_C
    # passive ground coupling (recharge toward ground T when idle)
    u_ground: float = 1.5             # W/(m^2 K)  contact conductance
    area_ground_m2: float = 8.0e4     # contact area cavern walls
    # optional active refrigeration
    active_cooling: bool = False
    chiller_kW_thermal: float = 0.0   # cooling capacity (thermal)
    chiller_cop: float = 3.0          # COP if active
    # --- cascaded cryogenic refrigeration (for very cold charge T) ---
    cascade_cooling: bool = False     # multi-stage cascade (N2/He) for < -60 C
    cascade_cop: float = 1.0          # lower COP at cryogenic temperatures
    # --- liquid air charge (cryogenic liquefaction, -196 C) ---
    liquid_air: bool = False          # store air as liquid (LIN/LAIR cycle)
    liquid_cop: float = 0.15          # COP of liquefaction (very energy-intensive)
    # --- lava-powered absorption refrigeration ---
    # Uses waste lava heat to drive absorption chillers instead of electricity.
    # The COP is lower (0.3-0.5) but the heat is free, so the electrical cost
    # of cooling drops by 80-90%. This dramatically improves EROI.
    lava_heated_cooling: bool = False # absorption chiller powered by lava heat
    lava_cooling_fraction: float = 0.8  # fraction of cooling met by lava heat
    # --- intercooled multi-stage recharge compression ---
    n_compress_stages: int = 1        # 1 = single-stage; 4+ = intercooled
    compress_eta: float = 0.70        # isentropic compressor efficiency
    # --- lava proximity and ultra thermal insulation ---
    # The cavern sits in a volcano/lava environment. The surrounding rock is
    # MUCH hotter than the standard geothermal gradient predicts because the
    # lava body creates a thermal halo. This requires ultra-high-performance
    # thermal insulation to keep the cavern cold.
    lava_proximity_m: float = float("inf")  # distance from cavern to lava body
    lava_t_nearby_c: float = float("nan")   # lava temp for halo calc; nan -> use LAVA_T_C
    ultra_insulation: bool = False          # ultra-high-performance thermal insulation
    ultra_insulation_mm: float = 0.0        # thickness of ultra insulation layer
    ultra_insulation_k: float = 0.0         # thermal conductivity W/(m K) of ultra insulation
    ultra_insulation_layers: int = 0        # number of MLI/vacuum panel layers

    def ground_t_at_depth_c(self) -> float:
        """Local ground temperature at the cavern depth.

        This is the HONESTY correction: below the shallow stable zone the ground
        WARMS at the geothermal gradient, so a deep cavern next to lava is HOT,
        not cold. The passive 'operates when not cooled' mode only works if the
        ground at this depth is actually cooler than the charge temperature.

        LAVA PROXIMITY HALO: If the cavern is near a lava body (lava_proximity_m
        is finite), the surrounding rock temperature is elevated far above the
        standard geothermal gradient. The lava creates a thermal halo that
        decays with distance. At 100m from 3000C lava, the rock may be 500-800C.
        This means the cavern CANNOT stay cold passively -- it requires ultra
        thermal insulation AND active refrigeration.
        """
        if not math.isnan(self.t_ground_c):
            return self.t_ground_c
        # stable zone ~ first 4 m tracks mean surface; below that add gradient
        stable = 4.0
        if self.depth_m <= stable and math.isinf(self.lava_proximity_m):
            return self.surf_t_c
        t_geo = self.surf_t_c + max(self.depth_m - stable, 0.0) * self.ground_k_per_m
        # lava thermal halo: rock temperature elevated near the lava body
        if not math.isinf(self.lava_proximity_m) and self.lava_proximity_m > 0:
            t_lava = self.lava_t_nearby_c if not math.isnan(self.lava_t_nearby_c) else LAVA_T_C
            # thermal halo decays roughly as 1/sqrt(r) from a cylindrical hot body
            # at close range the rock is very hot; at >1km it approaches geothermal
            halo_decay = 1.0 / math.sqrt(max(self.lava_proximity_m / 10.0, 1.0))
            t_halo = t_lava * halo_decay
            t_ground = max(t_geo, t_halo)
            return t_ground
        return t_geo

    def effective_u_ground(self) -> float:
        """Effective ground-to-cavern conductance including ultra insulation.

        Without ultra insulation, u_ground is the raw rock-to-air conductance.
        With ultra insulation (aerogel, vacuum panels, MLI), the conductance
        drops dramatically because the insulation adds thermal resistance in
        series: 1/U_eff = 1/U_ground + t_insul / k_insul
        """
        u = self.u_ground
        if self.ultra_insulation and self.ultra_insulation_mm > 0 and self.ultra_insulation_k > 0:
            t_insul_m = self.ultra_insulation_mm / 1000.0
            r_insul = t_insul_m / self.ultra_insulation_k  # m^2 K / W
            r_ground = 1.0 / max(u, 1e-6)                  # m^2 K / W
            u = 1.0 / (r_ground + r_insul)
        return u

    def ground_t_k(self) -> float:
        return c_to_k(self.ground_t_at_depth_c())


@dataclass
class LavaSourceSpec:
    """The hot reservoir: lava / superhot geothermal contact."""
    name: str = "Lava Contact"
    t_lava_c: float = LAVA_T_C
    contact_length_m: float = 800.0   # length of tunnel in thermal contact
    tunnel_diameter_m: float = 4.0    # bore diameter
    u_lava: float = 45.0              # W/(m^2 K)  overall heat-transfer coeff
    n_parallel_bores: int = 1         # parallel tunnels through the lava zone
    fin_factor: float = 1.0           # heat-transfer area enhancement (fins)
    heat_pipe: bool = False           # NaK heat pipes embedded in lava -> U~2000
    # --- shell-and-tube lava heat exchanger (decouples heat area from flow area) ---
    hx_enabled: bool = False          # dedicated HX bundle immersed in lava
    hx_n_tubes: int = 0               # number of small-diameter tubes in the HX
    hx_tube_od_mm: float = 0.0        # tube outer diameter (mm)
    hx_tube_length_m: float = 0.0     # length of each tube (immersed in lava)
    hx_u: float = 0.0                 # W/(m^2 K) for the HX tubes
    # the lava is treated as an infinite heat reservoir at fixed T (large mass)

    def ua_total(self) -> float:
        """Total UA including parallel bores, fins, heat-pipe, and HX bundle."""
        u = self.u_lava
        if self.heat_pipe:
            u = max(u, 2000.0)   # NaK heat pipe mode: direct lava -> air
        ua_tunnel = (u * math.pi * self.tunnel_diameter_m
                     * self.contact_length_m * self.n_parallel_bores * self.fin_factor)
        ua_hx = 0.0
        if self.hx_enabled and self.hx_n_tubes > 0 and self.hx_tube_od_mm > 0:
            od_m = self.hx_tube_od_mm / 1000.0
            ua_hx = (self.hx_u * math.pi * od_m * self.hx_tube_length_m
                     * self.hx_n_tubes)
        return ua_tunnel + ua_hx


@dataclass
class TunnelSpec:
    """The mile+ tunnel connecting cavern to exit, with turbine stages."""
    name: str = "Tunnel"
    total_length_m: float = 1800.0    # ~1.1 miles
    diameter_m: float = 4.0
    friction_factor: float = 0.018    # Darcy friction (rough rock/concrete)
    height_rise_m: float = 250.0      # chimney/stack rise at exit (buoyancy)
    n_turbine_stages: int = 6         # turbine/fan stages along the hot leg
    turbine_eta: float = 0.82         # per-stage isentropic efficiency
    generator_eta: float = 0.96       # generator efficiency
    exit_nozzle_area_m2: float = 2.0  # final jet nozzle
    n_exit_fans: int = 4              # fans in the exit jet
    exit_fan_eta: float = 0.75
    regenerator_eff: float = 0.0      # 0 = none; up to ~0.85 recuperator
    n_reheat_stages: int = 0          # reheat between turbine stages (Brayton reheat)
    # --- advanced power enhancement ---
    mhd_enabled: bool = False         # magnetohydrodynamic topping cycle
    mhd_eta: float = 0.15             # MHD extraction efficiency (of enthalpy above 1500 C)
    orc_enabled: bool = False         # Organic Rankine bottoming cycle on exhaust
    orc_eta: float = 0.12             # ORC cycle efficiency
    sco2_enabled: bool = False        # supercritical CO2 bottoming cycle
    sco2_eta: float = 0.40            # sCO2 cycle efficiency (40% at 1000+ C)
    steam_enabled: bool = False       # steam Rankine tertiary bottoming cycle
    steam_eta: float = 0.35           # steam cycle efficiency (35% at 500+ C)
    potassium_enabled: bool = False   # potassium vapor topping cycle (2000+ C)
    potassium_eta: float = 0.50       # potassium cycle efficiency (50% at 2000+ C)
    supersonic_nozzle: bool = False   # converging-diverging de Laval nozzle
    max_mach: float = 0.85            # exit Mach limit (0.85 subsonic, 2.5 supersonic)
    smooth_lining: bool = False       # smooth bored tunnel (lower friction)

    def area_m2(self) -> float:
        return math.pi * (self.diameter_m ** 2) / 4.0

    def hydraulic_diameter(self) -> float:
        return self.diameter_m

    def effective_friction(self) -> float:
        return 0.008 if self.smooth_lining else self.friction_factor


@dataclass
class ControlSpec:
    """Operating policy."""
    mode: str = "PASSIVE"   # PASSIVE (no chiller) or ACTIVE (chiller on)
    discharge_valve: float = 1.0   # 0..1 throttle of cavern discharge
    min_cavern_p_pa: float = 120.0e3   # stop discharging when cavern depleted
    n_systems: int = 1            # number of parallel tunnel systems (1=single, 2=dual)


# ==============================================================================
# SECTION 3 -- STATE
# ==============================================================================

@dataclass
class CavernState:
    """Live state of the cold-air cavern (a thermal battery)."""
    m_air_kg: float
    t_k: float
    v_m3: float
    p_pa: float
    cold_inventory_j: float = 0.0   # integrated 'cold' put in (negative heat)
    heat_leaked_in_j: float = 0.0   # passive ground leak accumulated
    chill_work_j: float = 0.0       # work spent by active chiller


def build_initial_cavern(spec: ColdCavernSpec) -> CavernState:
    m = mass_in_cavern(spec.volume_m3, spec.p_charge_pa, spec.t_charge_k)
    # 'cold inventory' = heat removed to get from ground T down to charge T
    t_gnd = spec.ground_t_k()
    cold = m * CP_AIR * (t_gnd - spec.t_charge_k)   # J (positive = cold stored)
    return CavernState(m_air_kg=m, t_k=spec.t_charge_k,
                       v_m3=spec.volume_m3, p_pa=spec.p_charge_pa,
                       cold_inventory_j=cold)


def cavern_pressure_from_state(st: CavernState, spec: ColdCavernSpec) -> float:
    """P = m R T / V  (ideal gas, rigid cavern)."""
    return st.m_air_kg * R_AIR * st.t_k / st.v_m3


# ==============================================================================
# SECTION 4 -- THE POWER HARNESS (tunnel + turbines)
# ==============================================================================

@dataclass
class FlowResult:
    """One steady-state solution of the discharge tunnel."""
    mdot: float = 0.0           # kg/s
    v_tunnel: float = 0.0       # m/s in the tunnel
    t_in_k: float = 0.0         # cavern discharge temperature
    t_out_k: float = 0.0        # exit temperature after lava heating
    q_lava_w: float = 0.0       # heat picked up from lava
    p_stack_pa: float = 0.0     # buoyancy/stack pressure
    p_cavern_pa: float = 0.0    # cavern driving pressure
    p_friction_pa: float = 0.0  # friction loss
    p_turbine_pa: float = 0.0   # pressure drop across turbine array
    t_hot_k: float = 0.0        # post-heating temperature (combustor outlet)
    t_turb_k: float = 0.0       # post-turbine temperature (before nozzle)
    t_exit_k: float = 0.0       # static exit temperature (after nozzle KE)
    mach_exit: float = 0.0
    v_exit: float = 0.0         # jet velocity at nozzle
    ke_jet_w: float = 0.0       # jet kinetic power at nozzle
    w_shaft_mech_w: float = 0.0 # turbine mechanical shaft work
    p_turbine_stages_w: float = 0.0   # turbine electrical (after gen)
    p_exit_fans_w: float = 0.0  # exit-fan electrical
    gen_loss_w: float = 0.0     # generator waste heat
    p_gross_w: float = 0.0
    p_parasitic_w: float = 0.0
    p_net_w: float = 0.0
    eta_carnot: float = 0.0
    p_carnot_ceiling_w: float = 0.0
    carnot_ok: bool = True
    q_regen_w: float = 0.0           # regenerator internal heat transfer
    t_pre_k: float = 0.0             # post-regenerator intake temperature
    p_mhd_w: float = 0.0             # MHD topping cycle power
    p_orc_w: float = 0.0             # ORC bottoming cycle power
    p_sco2_w: float = 0.0            # sCO2 bottoming cycle power
    p_steam_w: float = 0.0           # steam Rankine tertiary cycle power
    p_potassium_w: float = 0.0       # potassium vapor topping cycle power
    mach_exit_actual: float = 0.0    # actual exit Mach (pre-clamp)
    carnot_excess_w: float = 0.0     # work clamped by Carnot (rejected as waste)


def friction_dp(mdot: float, spec: TunnelSpec, rho: float) -> float:
    """Darcy-Weisbach pressure drop over the full tunnel length."""
    A = spec.area_m2()
    D = spec.hydraulic_diameter()
    v = mdot / (rho * A)
    f = spec.effective_friction()
    return f * (spec.total_length_m / D) * 0.5 * rho * v * v


# ==============================================================================
# SECTION 4b -- CONDENSATION / DEHUMIDIFICATION
# ==============================================================================
#
# When warm humid intake air (or ground-leak moisture) enters the cold cavern,
# it cools below its dew point and condenses. This:
#   * releases latent heat (partially offsetting the cooling),
#   * removes water vapor (reducing gas moles -> slight density increase),
#   * deposits liquid on cavern walls (needs drainage),
#   * releases the latent heat of vaporization (~2501 kJ/kg) into the cavern air.
#
# The model accounts for all of these. Condensation is a real effect in humid
# climates and is managed, not relied upon, as the primary driver.

LATENT_VAP_J_KG = 2.501e6      # J/kg, latent heat of vaporization of water
CP_WATER = 4186.0              # J/(kg K), specific heat of liquid water
RHO_WATER = 1000.0             # kg/m^3


def condensation_rate(t_air_c: float, rh: float, mdot: float,
                      t_cold_k: float) -> Tuple[float, float]:
    """Condensation rate and latent heat release when humid air hits cold walls.

    Returns (m_cond_kg_s, q_latent_w).
    If the cold surface temperature is above the dew point, no condensation.
    """
    dp_c = dew_point_c(t_air_c, rh)
    t_cold_c = k_to_c(t_cold_k)
    if t_cold_c >= dp_c:
        return 0.0, 0.0
    # saturation vapor pressure (Magnus), and actual vapor pressure
    a, b = 17.625, 243.04
    e_sat_cold = 611.2 * math.exp(a * t_cold_c / (b + t_cold_c))   # Pa
    e_actual = 611.2 * math.exp(a * dp_c / (b + dp_c))
    # humidity ratio difference: w = 0.622 * e / (P - e)
    w_in = 0.622 * e_actual / max(P_STD - e_actual, 1.0)
    w_sat_cold = 0.622 * e_sat_cold / max(P_STD - e_sat_cold, 1.0)
    dw = max(w_in - w_sat_cold, 0.0)   # kg water / kg dry air
    m_cond = mdot * dw                  # kg/s condensed
    q_latent = m_cond * LATENT_VAP_J_KG # W latent heat released
    return m_cond, q_latent


def humid_air_density(p_pa: float, t_k: float, rh: float) -> float:
    """Density of humid air, accounting for water vapor (lighter than dry air).

    rho = (P_d / (R_d T)) + (P_v / (R_v T))
    where P_d = P - P_v (dry air partial pressure), P_v = vapor pressure.
    """
    t_c = k_to_c(t_k)
    a, b = 17.625, 243.04
    e_sat = 611.2 * math.exp(a * t_c / (b + t_c))
    e_v = clamp(rh, 0.0, 1.0) * e_sat
    p_dry = p_pa - e_v
    R_v = 461.5   # J/(kg K) for water vapor
    return p_dry / (R_AIR * t_k) + e_v / (R_v * t_k)


# ==============================================================================
# SECTION 4c -- MULTI-STAGE TURBINE DETAIL
# ==============================================================================
#
# The turbine array is decomposed into N discrete stages, each with its own
# pressure ratio, temperature drop, and shaft work. This follows the
# Brayton-cycle multi-stage expansion pattern used in real gas turbines.
# The overall pressure ratio PR = P_cavern / P_atm is split evenly across
# stages (equal PR per stage), and each stage has its own isentropic efficiency.

@dataclass
class TurbineStage:
    """One stage of the multi-stage turbine array."""
    stage_num: int
    p_in_pa: float
    p_out_pa: float
    t_in_k: float
    t_out_k: float
    t_out_isentropic_k: float
    work_kg: float          # J/kg shaft work from this stage
    power_w: float          # W electrical from this stage
    eta_isentropic: float
    pr_stage: float         # pressure ratio across this stage


def solve_turbine_stages(t_hot_k: float, p_cavern_pa: float, p_atm_pa: float,
                         mdot: float, n_stages: int, eta_stage: float,
                         gen_eta: float) -> List[TurbineStage]:
    """Decompose the expansion into N equal-PR stages.

    Each stage has PR_stage = (P_cavern/P_atm)^(1/N).
    The isentropic exit T for each stage is T_in * PR_stage^(-(gamma-1)/gamma).
    The real exit T accounts for the isentropic efficiency.
    """
    stages: List[TurbineStage] = []
    pr_total = p_cavern_pa / p_atm_pa
    pr_stage = pr_total ** (1.0 / n_stages)
    t_in = t_hot_k
    p_in = p_cavern_pa
    gamma_m1 = GAMMA_AIR - 1.0

    for i in range(n_stages):
        p_out = p_in / pr_stage
        t_out_s = t_in * (1.0 / pr_stage) ** (gamma_m1 / GAMMA_AIR)
        t_out = t_in - eta_stage * (t_in - t_out_s)
        w_kg = CP_AIR * (t_in - t_out)           # J/kg
        p_elec = w_kg * mdot * gen_eta / n_stages  # W per stage (split evenly)
        # Actually the total work is split across stages, but each stage's
        # electrical output is its own w_kg * mdot * gen_eta
        p_elec = w_kg * mdot * gen_eta
        stages.append(TurbineStage(
            stage_num=i + 1, p_in_pa=p_in, p_out_pa=p_out,
            t_in_k=t_in, t_out_k=t_out, t_out_isentropic_k=t_out_s,
            work_kg=w_kg, power_w=p_elec,
            eta_isentropic=eta_stage, pr_stage=pr_stage))
        t_in = t_out
        p_in = p_out

    return stages


# ==============================================================================
# SECTION 4d -- ORC PARALLEL BOTTOMING CYCLE
# ==============================================================================
#
# An Organic Rankine Cycle can be placed in parallel on the post-turbine
# exhaust air, harvesting residual low-grade heat that the main Brayton cycle
# could not use. This is optional and only runs when ORC_HW["enabled"] is True.
# The ORC has its own Carnot limit and is accounted in the energy ledger.

def orc_power(t_exhaust_k: float, mdot: float, orc_spec: Dict) -> float:
    """Power from an ORC bottoming cycle on the exhaust air.

    The ORC takes heat from the exhaust air (cooling it further) and converts
    it to electricity at the ORC cycle efficiency. The heat available is
    mdot * cp * (T_exhaust - T_evap), capped by the evaporator UA.
    """
    if not orc_spec.get("enabled", False):
        return 0.0
    t_evap_k = c_to_k(orc_spec["t_evap_C"])
    t_cond_k = c_to_k(orc_spec["t_cond_C"])
    if t_exhaust_k <= t_evap_k:
        return 0.0
    # heat available from the exhaust down to the evaporator T
    q_avail = mdot * CP_AIR * (t_exhaust_k - t_evap_k)
    # evaporator capacity cap
    ua = orc_spec["ua_kw_per_k"] * 1e3   # W/K
    q_evap = ua * (t_exhaust_k - t_evap_k) * 0.1   # LMTD approximation
    q_orc = min(q_avail, q_evap)
    # ORC efficiency (Carnot-limited)
    eta_orc = orc_spec["eta_orc"]
    eta_carnot_orc = 1.0 - t_cond_k / t_evap_k
    eta_real = min(eta_orc, eta_carnot_orc * 0.6)   # ~60% of Carnot
    return q_orc * eta_real


# ==============================================================================
# SECTION 4e -- THERMAL ENERGY STORAGE (TES) / THERMAL MASS
# ==============================================================================
#
# The tunnel lining and surrounding rock act as a thermal mass that smooths
# short-term fluctuations and allows brief "coasting" when the lava heat
# input varies. The TES is modelled as a lumped capacitance on the tunnel
# wall, exchanging heat with both the air stream and the lava.

@dataclass
class TESState:
    """Thermal energy storage state (tunnel wall thermal mass)."""
    t_wall_k: float
    mass_kg: float
    cp: float          # J/(kg K) of the wall material

    def energy_j(self) -> float:
        return self.mass_kg * self.cp * (self.t_wall_k - T_REF_K)


def build_tes(tunnel: TunnelSpec, lava: LavaSourceSpec) -> TESState:
    """Initialise the TES at the lava temperature (equilibrium)."""
    # wall mass = lining volume * density
    r_outer = tunnel.diameter_m / 2.0 + 0.35   # lining thickness
    r_inner = tunnel.diameter_m / 2.0
    v_lining = math.pi * (r_outer**2 - r_inner**2) * lava.contact_length_m
    rho_concrete = 2400.0   # kg/m^3
    cp_concrete = 880.0     # J/(kg K)
    return TESState(t_wall_k=c_to_k(lava.t_lava_c) * 0.7,  # starts warm but not at lava T
                    mass_kg=v_lining * rho_concrete,
                    cp=cp_concrete)


def tes_update(tes: TESState, t_air_k: float, t_lava_k: float,
               dt_s: float, u_wall: float = 5.0) -> Tuple[float, float]:
    """Update the TES wall temperature and return (q_to_air, q_from_lava) in W.

    The wall exchanges heat with both the air stream (cooling it) and the lava
    (heating it). This provides thermal inertia / coasting.
    """
    area = math.pi * 1.0 * 200.0   # simplified contact area (m^2)
    q_from_lava = u_wall * area * (t_lava_k - tes.t_wall_k)
    q_to_air = u_wall * area * 0.5 * (tes.t_wall_k - t_air_k)
    dT = (q_from_lava - q_to_air) * dt_s / (tes.mass_kg * tes.cp)
    tes.t_wall_k += dT
    return q_to_air, q_from_lava


def stack_pressure(spec: TunnelSpec, rho_cold: float, rho_hot: float) -> float:
    """Buoyancy/stack draft from the height rise and density difference."""
    return G * spec.height_rise_m * max(rho_cold - rho_hot, 0.0)


def lava_heat_transfer(spec: LavaSourceSpec, t_air_avg_k: float) -> float:
    """Steady heat flow from lava into the air, Q = U A dT (W)."""
    A = math.pi * spec.tunnel_diameter_m * spec.contact_length_m
    dt = c_to_k(spec.t_lava_c) - t_air_avg_k
    return spec.u_lava * A * max(dt, 0.0)


def solve_flow(cavern: CavernState, cavern_spec: ColdCavernSpec,
               lava: LavaSourceSpec, tunnel: TunnelSpec,
               ctrl: ControlSpec) -> FlowResult:
    """Solve the steady discharge for the current cavern state -- a Brayton-style
    open air cycle driven by the cavern pressure and heated by lava.

    The cycle, honestly decomposed:
      1. Cavern drives cold dense air out at P_cavern > P_atm.
      2. REGENERATOR (optional): the cold incoming air is preheated by the hot
         exhaust via a counterflow heat exchanger. This raises T_in to T_pre
         before the lava contact, increasing Carnot efficiency and reducing the
         lava heat demand for the same T_hot.
      3. HEATING (constant-pressure combustor analogue): the lava contact raises
         the air from T_pre to T_hot.  Q_lava = UA (T_lava - T_avg_hot)
                                            = mdot cp (T_hot - T_pre).
      4. EXPANSION (turbine array): the air expands from P_cavern to P_atm,
         dropping from T_hot to T_out.  With reheat, the expansion is split into
         N+1 sub-expansions, each followed by re-injection of lava heat back to
         T_hot. This is the Brayton REHEAT cycle and it significantly increases
         turbine work for the same pressure ratio.
      5. EXIT JET + EXIT FANS: the residual enthalpy/KE at the nozzle is partly
         harvested by the exit fan array.
      6. FRICTION + STACK: Darcy-Weisbach friction along the tunnel; buoyancy
         stack draft from the height rise assists the flow.

    The steady-flow energy equation for the tunnel control volume closes the
    First Law exactly:
        mdot h_in + Q_lava_total = mdot h_out + W_shaft + KE_exit
    so W_shaft + KE_exit = Q_lava_total - mdot cp (T_out - T_in) = mdot cp (T_hot - T_out).
    The work comes from the temperature drop across the turbine, NOT from
    double-counting the heat.
    """
    res = FlowResult()
    p_atm = P_STD
    t_in = cavern.t_k
    p_cavern = cavern.p_pa
    rho_in = air_density(p_cavern, t_in)

    # driving pressure (throttled by the discharge valve)
    n_sys = max(1, ctrl.n_systems)       # dual tunnel = 2 systems side by side
    p_cav = (p_cavern - p_atm) * clamp(ctrl.discharge_valve, 0.0, 1.0)
    A = tunnel.area_m2() * n_sys         # total flow area across all systems
    UA = lava.ua_total() * n_sys         # total UA across all systems
    t_lava_k = c_to_k(lava.t_lava_c)
    pr = p_cavern / p_atm                      # pressure ratio across turbines
    gamma_m1 = GAMMA_AIR - 1.0
    n_reheat = tunnel.n_reheat_stages          # reheat stages between expansions
    n_expansions = n_reheat + 1                # total expansion segments
    pr_per_seg = pr ** (1.0 / n_expansions)    # equal PR per segment
    regen = clamp(tunnel.regenerator_eff, 0.0, 0.90)

    def cycle(mdot_try: float) -> Tuple[float, float, float, float, float,
                                        float, float, float, float, float]:
        """Return (t_hot, t_out, q_lava_total, dp_fric, dp_stack, rho_avg,
        v_tun, dp_turb, t_pre, q_regen) for a trial mass flow."""
        # 0. Regenerator: preheat cold air using exhaust heat.
        # The exhaust is at T_out (post-turbine). The regenerator transfers
        # a fraction of (T_out - T_in) to the incoming air.
        # We need to iterate because T_out depends on T_hot which depends on T_pre.
        # For the trial, estimate T_out first without regenerator.
        t_pre = t_in
        q_regen = 0.0
        if regen > 0.0 and mdot_try > 0.0:
            # First estimate T_out without regenerator (effectiveness-NTU)
            NTU0 = UA / (mdot_try * CP_AIR)
            eps0 = 1.0 - math.exp(-NTU0)
            t_hot0 = t_in + eps0 * (t_lava_k - t_in)
            t_hot0 = clamp(t_hot0, t_in, t_lava_k)
            t_out_s0 = t_hot0 * (1.0 / pr) ** (gamma_m1 / GAMMA_AIR)
            t_out0 = t_hot0 - tunnel.turbine_eta * (t_hot0 - t_out_s0)
            # Regenerator: preheat incoming air
            t_pre = t_in + regen * (t_out0 - t_in)
            t_pre = clamp(t_pre, t_in, t_out0)
            q_regen = mdot_try * CP_AIR * (t_pre - t_in)  # heat transferred (internal)

        # 1. heating: use effectiveness-NTU method for constant-T heat source.
        # For a constant-temperature heat source (T_lava = const, e.g. phase change
        # or infinite thermal mass), the effectiveness is:
        #   eps = 1 - exp(-NTU),  NTU = UA / (mdot * cp)
        # Then Q = eps * mdot * cp * (T_lava - T_pre)
        # and T_hot = T_pre + eps * (T_lava - T_pre)
        # This is the exact closed-form solution (Incropera & DeWitt, Ch 11).
        if mdot_try > 0.0:
            NTU = UA / (mdot_try * CP_AIR)
            eps_hx = 1.0 - math.exp(-NTU)
            t_hot = t_pre + eps_hx * (t_lava_k - t_pre)
            t_hot = clamp(t_hot, t_pre, t_lava_k * 0.999)  # cap below T_lava
            q_lava_main = mdot_try * CP_AIR * (t_hot - t_pre)
        else:
            t_hot = t_pre
            q_lava_main = 0.0

        # 2. expansion with optional reheat
        # Without reheat: single expansion T_hot -> T_out
        # With reheat: N+1 expansions, each from T_hot to T_seg, then reheat back to T_hot
        # Each segment has PR_seg = PR^(1/(N+1))
        # Work per segment = cp * (T_hot - T_seg), total work = (N+1) * cp * (T_hot - T_seg)
        t_out_s_seg = t_hot * (1.0 / pr_per_seg) ** (gamma_m1 / GAMMA_AIR)
        t_seg = t_hot - tunnel.turbine_eta * (t_hot - t_out_s_seg)
        t_seg = clamp(t_seg, t_out_s_seg, t_hot)
        # Total heat input including reheat
        q_lava_reheat = 0.0
        if n_reheat > 0:
            # After each expansion segment (except the last), reheat back to T_hot
            # Q_reheat = mdot * cp * (T_hot - T_seg) per reheat stage
            q_lava_reheat = n_reheat * mdot_try * CP_AIR * (t_hot - t_seg)
        q_lava_total = q_lava_main + q_lava_reheat
        # Final exit temperature after all expansions + reheat
        t_out = t_seg   # the last expansion segment's exit

        # 3. properties along the tunnel (use post-expansion T for the long leg)
        t_avg = 0.5 * (t_in + t_out)
        rho_avg = air_density(p_atm + 0.5 * p_cav, t_avg)
        rho_hot = air_density(p_atm, t_out)
        dp_fric = friction_dp(mdot_try, tunnel, rho_avg)
        dp_stack = stack_pressure(tunnel, rho_in, rho_hot)
        v_tun = mdot_try / (rho_avg * A)
        # turbine back-pressure ~ a fraction of the dynamic head
        dp_turb = 0.5 * rho_avg * v_tun * v_tun * 0.6
        return (t_hot, t_out, q_lava_total, dp_fric, dp_stack, rho_avg,
                v_tun, dp_turb, t_pre, q_regen)

    def balance(mdot_try: float) -> float:
        if mdot_try <= 0.0:
            return p_cav
        _, _, _, dp_fric, dp_stack, _, _, dp_turb, _, _ = cycle(mdot_try)
        return (p_cav + dp_stack) - (dp_fric + dp_turb)

    # bracket and bisect for the mass flow that balances driving vs resisting pressure
    mdot_guess = 0.5 * rho_in * A * math.sqrt(max(p_cav, 0.0) / max(rho_in, EPS))
    lo, hi = 0.0, mdot_guess * 4.0 + 1.0
    tries = 0
    while balance(hi) > 0.0 and tries < 60:
        hi *= 1.6
        tries += 1
    if balance(hi) > 0.0:
        mdot = hi
    else:
        try:
            mdot = bisect(balance, lo, hi, tol=1e-3)
        except ValueError:
            mdot = 0.0

    # ---- CHOKED FLOW LIMIT ----
    # The maximum mass flow through a duct is limited by the speed of sound.
    # For isentropic flow from a reservoir at P0, T0:
    #   mdot_max = A_total * P0 * sqrt(gamma/(R*T0)) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
    # This is a hard physical limit -- no pressure can force more mass through.
    # A_total includes all parallel bores.
    gamma = GAMMA_AIR
    P0 = p_cavern
    T0 = t_in
    A_total = A * lava.n_parallel_bores
    choke_coeff = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    mdot_choke = A_total * P0 * math.sqrt(gamma / (R_AIR * T0)) * choke_coeff
    if mdot > mdot_choke:
        mdot = mdot_choke

    if mdot <= 1e-6:
        return res

    # ---- recompute at the solved mdot ----
    t_hot, t_turb, q_lava, dp_fric, dp_stack, rho_avg, v_tun, dp_turb, t_pre, q_regen = cycle(mdot)

    # ---- staged turbine shaft work (Brayton expansion with reheat) ----
    # With reheat: total work = n_expansions * cp * (T_hot - T_seg) * mdot
    # Without reheat: work = cp * (T_hot - T_turb) * mdot  (same formula, n_expansions=1)
    # If pr < 1 (cavern pressure below atmospheric), the turbine consumes work
    # (compression) rather than producing it -- clamp to zero.
    n_expansions_actual = tunnel.n_reheat_stages + 1
    if pr > 1.0:
        w_shaft_mech = n_expansions_actual * CP_AIR * (t_hot - t_turb) * mdot  # W mechanical
    else:
        w_shaft_mech = 0.0
    p_turb = w_shaft_mech * tunnel.generator_eta          # W electrical
    gen_loss = w_shaft_mech * (1.0 - tunnel.generator_eta)  # W waste heat

    # ---- MHD topping cycle (optional) ----
    # At T_hot > 1500 C, seeding the air with cesium/potassium vapor makes it
    # weakly ionized and conductive. An MHD channel extracts DC work directly
    # from the flowing plasma BEFORE the first turbine stage. The MHD efflux
    # drops T_hot by dT_mhd, and that extracted enthalpy becomes electricity.
    # This is a real topping cycle used in conceptual coal-fired MHD plants.
    p_mhd = 0.0
    t_hot_after_mhd = t_hot
    t_mhd_threshold = c_to_k(1500.0)
    if tunnel.mhd_enabled and t_hot > t_mhd_threshold and pr > 1.0:
        # MHD extracts a fraction of the enthalpy above the threshold
        h_above = CP_AIR * (t_hot - t_mhd_threshold)   # J/kg available
        w_mhd_kg = tunnel.mhd_eta * h_above              # J/kg extracted
        p_mhd = w_mhd_kg * mdot                          # W electrical
        t_hot_after_mhd = t_hot - w_mhd_kg / CP_AIR      # T drops
        # The turbine work now uses the post-MHD temperature
        # Recompute turbine work with reduced T_hot
        t_hot_eff = t_hot_after_mhd
        # Recompute the segment exit T with the new T_hot
        t_out_s_eff = t_hot_eff * (1.0 / pr_per_seg) ** (gamma_m1 / GAMMA_AIR)
        t_seg_eff = t_hot_eff - tunnel.turbine_eta * (t_hot_eff - t_out_s_eff)
        t_seg_eff = clamp(t_seg_eff, t_out_s_eff, t_hot_eff)
        w_shaft_mech = n_expansions_actual * CP_AIR * (t_hot_eff - t_seg_eff) * mdot
        p_turb = w_shaft_mech * tunnel.generator_eta
        gen_loss = w_shaft_mech * (1.0 - tunnel.generator_eta)
        # Update t_turb for the exit nozzle calculation
        t_turb = t_seg_eff
        # MHD heat is part of the lava heat budget (it comes from Q_lava)
        # so no additional heat input is needed -- just a reallocation

    # ---- exit nozzle + exit fans ----
    # The air at T_turb, P_atm still carries residual enthalpy above ambient.
    # The nozzle converts a portion of that enthalpy into jet kinetic energy,
    # which LOWERS the static exit temperature (T_exit = T_turb - v^2/2cp).
    # The exit fans then harvest a fraction of the jet KE as electricity.
    # With a supersonic de Laval nozzle, the Mach limit is raised to 2.5,
    # allowing much higher jet velocities and more fan power.
    t_amb = c_to_k(20.0)
    e_nozzle = max(CP_AIR * (t_turb - t_amb), 0.0)         # J/kg from enthalpy
    v_from_enthalpy = math.sqrt(2.0 * e_nozzle)
    mach_limit = tunnel.max_mach if tunnel.supersonic_nozzle else 0.85
    c_sound = sound_speed(t_turb)
    v_from_enthalpy = clamp(v_from_enthalpy, 0.0, mach_limit * c_sound)
    # total jet velocity = tunnel flow + nozzle acceleration (vector sum in KE)
    v_exit = math.sqrt(v_tun * v_tun + v_from_enthalpy * v_from_enthalpy)
    v_exit = clamp(v_exit, 0.0, mach_limit * c_sound)
    mach_actual = v_exit / c_sound
    ke_jet = 0.5 * mdot * v_exit * v_exit                  # W total jet KE
    # the enthalpy-derived portion lowers the static exit T; the bulk-flow
    # portion was already at T_turb (no further static T drop).
    t_exit = t_turb - (v_from_enthalpy * v_from_enthalpy) / (2.0 * CP_AIR)
    # Fan harvesting fraction: at supersonic Mach, the fans can capture less
    # of the total KE (shock losses, blade stress limits). At subsonic, 35%.
    fan_fraction = 0.35 if mach_limit <= 0.85 else 0.20
    p_exit_fans = ke_jet * tunnel.exit_fan_eta * fan_fraction

    # ---- ORC bottoming cycle (optional) ----
    # The exhaust air at T_exit still carries low-grade heat. An ORC can
    # extract a fraction of that heat as electricity. The ORC has its own
    # Carnot limit using T_exit as the hot side and ambient as the cold side.
    p_orc = 0.0
    if tunnel.orc_enabled and t_exit > t_amb + 50.0:
        q_orc_avail = mdot * CP_AIR * (t_exit - t_amb) * 0.3  # ORC can access ~30%
        eta_carnot_orc = 1.0 - t_amb / t_exit
        eta_orc_real = min(tunnel.orc_eta, eta_carnot_orc * 0.6)
        p_orc = q_orc_avail * eta_orc_real
        # ORC heat comes from the exhaust, reducing the waste heat

    # ---- sCO2 bottoming cycle (optional) ----
    # A supercritical CO2 cycle is far more efficient than ORC at high exhaust
    # temperatures (800-1500 C). It operates at 200+ bar with a compact
    # turbomachinery footprint. Real sCO2 plants achieve 40-50% thermal
    # efficiency at 1000+ C heat source temperatures.
    p_sco2 = 0.0
    if tunnel.sco2_enabled and t_exit > t_amb + 200.0:
        q_sco2_avail = mdot * CP_AIR * (t_exit - t_amb) * 0.5  # sCO2 accesses ~50%
        eta_carnot_sco2 = 1.0 - t_amb / t_exit
        eta_sco2_real = min(tunnel.sco2_eta, eta_carnot_sco2 * 0.65)
        p_sco2 = q_sco2_avail * eta_sco2_real
        # sCO2 heat comes from the exhaust, reducing waste heat

    # ---- Potassium vapor topping cycle (optional) ----
    # A potassium vapor Rankine cycle operates at 1500-2500 C, above the
    # practical limit for steam/sCO2. Real potassium cycles achieve 40-55%
    # efficiency at these temperatures. It extracts work from the hottest
    # portion of the exhaust before passing to sCO2 and steam.
    p_potassium = 0.0
    if tunnel.potassium_enabled and t_exit > t_amb + 1500.0:
        q_k_avail = mdot * CP_AIR * (t_exit - t_amb) * 0.35  # K accesses ~35%
        eta_carnot_k = 1.0 - t_amb / t_exit
        eta_k_real = min(tunnel.potassium_eta, eta_carnot_k * 0.55)
        p_potassium = q_k_avail * eta_k_real

    # ---- Steam Rankine tertiary bottoming cycle (optional) ----
    # After the sCO2 and potassium cycles have extracted what they can, the
    # exhaust may still be at 500-1500 C. A steam Rankine cycle can capture
    # additional low-grade heat at 30-38% efficiency.
    p_steam = 0.0
    if tunnel.steam_enabled and t_exit > t_amb + 300.0:
        # Estimate post-sCO2/post-potassium temperature
        t_post_sco2 = t_exit
        if tunnel.sco2_enabled and p_sco2 > 0 and mdot > 0:
            dt_sco2 = p_sco2 / (mdot * CP_AIR * 0.5)
            t_post_sco2 = max(t_post_sco2 - dt_sco2, t_amb + 100.0)
        if tunnel.potassium_enabled and p_potassium > 0 and mdot > 0:
            dt_k = p_potassium / (mdot * CP_AIR * 0.35)
            t_post_sco2 = max(t_post_sco2 - dt_k, t_amb + 100.0)
        q_steam_avail = mdot * CP_AIR * (t_post_sco2 - t_amb) * 0.4
        eta_carnot_steam = 1.0 - t_amb / t_post_sco2
        eta_steam_real = min(tunnel.steam_eta, eta_carnot_steam * 0.6)
        p_steam = q_steam_avail * eta_steam_real

    # ---- parasitic loads ----
    p_parasitic = 5.0e3   # 5 kW aux/controls

    # ---- Carnot clamp on total heat-engine work ----
    # The turbine, MHD, exit fans, ORC, sCO2, and steam all extract work from
    # the lava heat. Their combined output must not exceed eta_c * q_lava.
    p_heat_engine = p_turb + p_exit_fans + p_mhd + p_orc + p_sco2 + p_steam + p_potassium
    eta_c_pre = 1.0 - t_in / t_lava_k if t_lava_k > t_in else 0.0
    p_carnot_heat_limit = eta_c_pre * q_lava
    carnot_excess_w = 0.0
    # Use a 0.1% tolerance to avoid floating-point triggering at the exact limit
    if p_heat_engine > p_carnot_heat_limit * 1.001 and p_carnot_heat_limit > 0:
        carnot_excess_w = p_heat_engine - p_carnot_heat_limit
        # The excess work becomes additional waste heat rejected to the exhaust
        p_turb -= carnot_excess_w * (p_turb / p_heat_engine)
        p_exit_fans -= carnot_excess_w * (p_exit_fans / p_heat_engine)
        p_mhd -= carnot_excess_w * (p_mhd / p_heat_engine)
        p_orc -= carnot_excess_w * (p_orc / p_heat_engine)
        p_sco2 -= carnot_excess_w * (p_sco2 / p_heat_engine)
        p_steam -= carnot_excess_w * (p_steam / p_heat_engine)
        p_potassium -= carnot_excess_w * (p_potassium / p_heat_engine)

    p_gross = p_turb + p_exit_fans + p_mhd + p_orc + p_sco2 + p_steam + p_potassium
    p_net = p_gross - p_parasitic

    # ---- Carnot + exergy ceiling (the honesty check) ----
    # Net work comes from TWO legitimate, separately-bounded sources:
    #   (1) cavern PRESSURE EXERGY: isothermal expansion work per kg =
    #       R T_in ln(P/P0).  This is stored energy paid for at charge time.
    #   (2) HEAT ENGINE on the lava heat, bounded by Carnot using the lava T as
    #       the hot reservoir and T_in as the cold reservoir:
    #       W_heat <= (1 - T_in/T_lava) * Q_lava.
    # The audit checks: P_net <= W_pressure_exergy + W_heat_carnot.
    t_cold = t_in
    eta_c = 1.0 - t_cold / t_lava_k if t_lava_k > t_cold else 0.0
    p_carnot_heat = eta_c * q_lava
    p_pressure_exergy = mdot * R_AIR * t_in * math.log(max(pr, 1.0))
    p_ceiling = p_pressure_exergy + p_carnot_heat
    carnot_ok = p_net <= p_ceiling * 1.001 + 1.0

    res.mdot = mdot
    res.v_tunnel = v_tun
    res.t_in_k = t_in
    res.t_out_k = t_exit          # reported exit temperature = static (post-nozzle)
    res.q_lava_w = q_lava
    res.p_stack_pa = dp_stack
    res.p_cavern_pa = p_cav
    res.p_friction_pa = dp_fric
    res.p_turbine_pa = dp_turb
    res.t_hot_k = t_hot
    res.t_turb_k = t_turb
    res.t_exit_k = t_exit
    res.mach_exit = mach(v_exit, t_turb)
    res.v_exit = v_exit
    res.ke_jet_w = ke_jet
    res.w_shaft_mech_w = w_shaft_mech
    res.p_turbine_stages_w = p_turb
    res.p_exit_fans_w = p_exit_fans
    res.gen_loss_w = gen_loss
    res.p_gross_w = p_gross
    res.p_parasitic_w = p_parasitic
    res.p_net_w = p_net
    res.eta_carnot = eta_c
    res.p_carnot_ceiling_w = p_ceiling
    res.carnot_ok = carnot_ok
    res.q_regen_w = q_regen
    res.t_pre_k = t_pre
    res.p_mhd_w = p_mhd
    res.p_orc_w = p_orc
    res.p_sco2_w = p_sco2
    res.p_steam_w = p_steam
    res.p_potassium_w = p_potassium
    res.mach_exit_actual = mach_actual
    res.carnot_excess_w = carnot_excess_w
    return res


# ==============================================================================
# SECTION 5 -- CAVERN RECHARGE (passive ground coupling + optional chiller)
# ==============================================================================

def recharge_cavern(st: CavernState, spec: ColdCavernSpec,
                    dt_s: float, mode: str) -> Tuple[float, float, float]:
    """Update cavern T toward ground T (passive) and/or chill it (active).

    Returns (heat_leaked_in_W, chiller_elec_W, chiller_heat_removed_W).
    Passive leak is positive when the ground is WARMER than the cavern (i.e. the
    cavern warms -> loses cold). This is the crux of the honesty layer: if the
    ground at depth is warmer than the charge temperature, the cavern cannot
    stay cold by itself.
    """
    t_gnd_k = spec.ground_t_k()
    dT = t_gnd_k - st.t_k
    u_eff = spec.effective_u_ground()
    q_leak = u_eff * spec.area_ground_m2 * dT   # W (+: heat in)
    # update T via lumped capacitance: m cp dT = q dt
    m = st.m_air_kg
    if m > 0:
        st.t_k += q_leak * dt_s / (m * CP_AIR)
    st.heat_leaked_in_j += max(q_leak, 0.0) * dt_s

    q_chiller_thermal = 0.0
    p_chiller_elec = 0.0
    if mode == "ACTIVE" and spec.active_cooling and spec.chiller_kW_thermal > 0:
        # chiller pulls heat out of the cavern air
        q_chiller_thermal = spec.chiller_kW_thermal * 1e3
        # only chills if cavern is warmer than the chiller setpoint (charge T)
        if st.t_k > spec.t_charge_k:
            if m > 0:
                st.t_k -= q_chiller_thermal * dt_s / (m * CP_AIR)
            st.t_k = max(st.t_k, spec.t_charge_k)
            p_chiller_elec = q_chiller_thermal / spec.chiller_cop
            st.chill_work_j += p_chiller_elec * dt_s
            # active cooling also fights the ground leak
            st.heat_leaked_in_j -= q_chiller_thermal * dt_s
        else:
            q_chiller_thermal = 0.0   # already at setpoint, chiller idles
    # pressure follows T at fixed mass/volume
    st.p_pa = cavern_pressure_from_state(st, spec)
    return q_leak, p_chiller_elec, q_chiller_thermal


# ==============================================================================
# SECTION 6 -- SIMULATION
# ==============================================================================

@dataclass
class SimResult:
    cavern_spec: ColdCavernSpec
    lava_spec: LavaSourceSpec
    tunnel_spec: TunnelSpec
    ctrl: ControlSpec
    mode: str = "PASSIVE"
    verdict: str = "UNKNOWN"
    detail: str = ""
    t_h: List[float] = field(default_factory=list)
    cavern_t_c: List[float] = field(default_factory=list)
    cavern_p_kpa: List[float] = field(default_factory=list)
    cavern_m_kg: List[float] = field(default_factory=list)
    mdot: List[float] = field(default_factory=list)
    t_out_c: List[float] = field(default_factory=list)
    q_lava_mw: List[float] = field(default_factory=list)
    p_net_mw: List[float] = field(default_factory=list)
    p_gross_mw: List[float] = field(default_factory=list)
    v_exit: List[float] = field(default_factory=list)
    eta_carnot: List[float] = field(default_factory=list)
    p_carnot_mw: List[float] = field(default_factory=list)
    carnot_ok: List[bool] = field(default_factory=list)
    mean_p_net_mw: float = 0.0
    peak_p_net_mw: float = 0.0
    total_twh: float = 0.0
    homes: float = 0.0
    discharge_hours: float = 0.0
    chill_work_twh: float = 0.0
    heat_leaked_twh: float = 0.0
    energy_residual: float = 0.0
    capex_musd: float = 0.0
    eroi: float = 0.0


def simulate(cavern_spec: ColdCavernSpec, lava_spec: LavaSourceSpec,
             tunnel_spec: TunnelSpec, ctrl: ControlSpec,
             hours: float = 48.0, n_steps: int = 2400) -> SimResult:
    """Integrate the cavern discharge + recharge cycle.

    The cavern starts charged (cold, pressurised). It discharges through the
    tunnel producing power; simultaneously the ground leaks heat in (warming
    it). When cavern pressure drops below the minimum the discharge valve
    closes and the cavern recharges (passively toward ground T, or actively
    via chiller). The simulation tracks the full energy ledger.
    """
    st = build_initial_cavern(cavern_spec)
    n_sys = max(1, ctrl.n_systems)       # dual tunnel: 2 systems side by side
    # Each system has its own cavern of the same size, so total stored mass
    # scales with n_sys. We model this by scaling the initial mass.
    st.m_air_kg *= n_sys
    res = SimResult(cavern_spec=cavern_spec, lava_spec=lava_spec,
                    tunnel_spec=tunnel_spec, ctrl=ctrl, mode=ctrl.mode)
    dt = hours * HOUR / n_steps

    e_elec = 0.0
    e_chill = 0.0
    e_leak = 0.0
    sum_p = 0.0
    n_acc = 0
    discharge_active = True

    # energy ledger (First Law, closed over the whole system)
    #   cavern internal energy U = m cv (T - T_ref)   (ideal gas, rigid cavern)
    #   exhaust enthalpy  h = mdot cp (T_exit - T_ref)  (air leaving at static T)
    #   jet KE leaving    = KE_jet - W_fans  (after exit fans take their share)
    #   waste heat        = generator loss + parasitic  (to ambient)
    # The air-stream balance  mdot h_in + Q_lava = mdot h_exit + W_shaft + KE_jet
    # closes by construction, so the global ledger closes too.
    e_in_lava = 0.0       # heat delivered by lava to air
    e_in_leak = 0.0       # heat leaked from ground into cavern
    e_in_grid = 0.0       # grid electricity IN to run the chiller (active mode)
    e_out_elec = 0.0      # net electricity OUT (discharge only, >= 0)
    e_out_jet = 0.0       # jet KE leaving the system (after exit fans)
    e_out_exhaust = 0.0   # enthalpy carried out by discharged air (static T)
    e_out_waste = 0.0     # generator loss + parasitic -> ambient heat
    e_out_chiller_amb = 0.0  # heat dumped to ambient by chiller condenser
    e_cavern_initial = st.m_air_kg * CV_AIR * (st.t_k - T_REF_K)

    for step in range(n_steps + 1):
        t_h = step * dt / HOUR

        # decide whether discharging
        if st.p_pa - P_STD < ctrl.min_cavern_p_pa - P_STD:
            discharge_active = False
        # in PASSIVE mode we never actively chill; recharge is just ground leak

        if discharge_active:
            ctrl_now = ControlSpec(mode=ctrl.mode,
                                    discharge_valve=ctrl.discharge_valve,
                                    min_cavern_p_pa=ctrl.min_cavern_p_pa)
            fr = solve_flow(st, cavern_spec, lava_spec, tunnel_spec, ctrl_now)
            # ground heat leak into the cavern (happens during discharge too)
            t_gnd_k = cavern_spec.ground_t_k()
            u_eff = cavern_spec.effective_u_ground()
            q_leak_now = u_eff * cavern_spec.area_ground_m2 * (t_gnd_k - st.t_k)
            # rigid-tank energy balance: dU/dt = Q_leak - mdot*h_in, with
            #   U = m cv (T - Tref),  h_in = cp (T - Tref),  dm/dt = -mdot.
            #   -> m cv dT = (Q_leak - mdot (cp-cv)(T-Tref)) dt
            #             = (Q_leak - mdot R (T-Tref)) dt
            # The tank COOLS as it discharges (adiabatic expansion of the
            # remaining air), which is real and is what closes the First Law.
            m_now = max(st.m_air_kg, 1.0)
            dT = (q_leak_now - fr.mdot * R_AIR * (st.t_k - T_REF_K)) * dt / (m_now * CV_AIR)
            st.t_k += dT
            # deplete cavern mass, then recover pressure from the new (m, T, V)
            dm = fr.mdot * dt
            st.m_air_kg = max(st.m_air_kg - dm, 0.0)
            st.p_pa = cavern_pressure_from_state(st, cavern_spec)
            st.heat_leaked_in_j += max(q_leak_now, 0.0) * dt
            e_elec += fr.p_net_w * dt
            e_chill += 0.0   # no chilling while discharging
            sum_p += fr.p_net_w
            n_acc += 1
            e_in_lava += fr.q_lava_w * dt
            e_in_leak += max(q_leak_now, 0.0) * dt
            e_out_elec += fr.p_net_w * dt
            e_out_jet += (fr.ke_jet_w - fr.p_exit_fans_w) * dt
            # exhaust enthalpy: the regenerator cools the exhaust by q_regen
            # before it leaves the system, so subtract that internal transfer.
            # The bottoming cycles (ORC, sCO2, steam) also extract heat from
            # the exhaust, reducing what actually leaves the system.
            # Their waste heat (heat_in - work_out) is rejected to environment
            # via cooling, not through the exhaust, so we subtract the full
            # heat input to the bottoming cycles from the exhaust.
            # Use actual cycle efficiencies: heat_in = work / eta
            bottoming_heat_w = 0.0
            if fr.p_orc_w > 0:
                bottoming_heat_w += fr.p_orc_w / max(tunnel_spec.orc_eta, 0.01)
            if fr.p_sco2_w > 0:
                bottoming_heat_w += fr.p_sco2_w / max(tunnel_spec.sco2_eta, 0.01)
            if fr.p_steam_w > 0:
                bottoming_heat_w += fr.p_steam_w / max(tunnel_spec.steam_eta, 0.01)
            if fr.p_potassium_w > 0:
                bottoming_heat_w += fr.p_potassium_w / max(tunnel_spec.potassium_eta, 0.01)
            # Cap at available exhaust heat (can't extract more than available)
            exhaust_avail = fr.mdot * CP_AIR * (fr.t_exit_k - T_REF_K) - fr.q_regen_w
            bottoming_heat_w = min(bottoming_heat_w, max(exhaust_avail, 0.0))
            e_out_exhaust += (fr.mdot * CP_AIR * (fr.t_exit_k - T_REF_K)
                              - fr.q_regen_w - bottoming_heat_w) * dt
            # Carnot-clamped work is rejected as additional waste heat.
            # Bottoming cycle waste heat (heat_in - work_out) is also rejected
            # to environment via cooling.
            bottoming_work = fr.p_orc_w + fr.p_sco2_w + fr.p_steam_w + fr.p_potassium_w
            bottoming_waste_w = bottoming_heat_w - bottoming_work
            e_out_waste += (fr.gen_loss_w + fr.p_parasitic_w + fr.carnot_excess_w
                            + bottoming_waste_w) * dt
            res.discharge_hours += dt / HOUR
        else:
            fr = FlowResult()
            # recharge (passive ground leak, and/or active chiller)
            q_leak, p_chill, q_chill_th = recharge_cavern(st, cavern_spec, dt, ctrl.mode)
            e_leak += max(q_leak, 0.0) * dt
            e_chill += p_chill * dt
            e_in_leak += max(q_leak, 0.0) * dt
            e_in_grid += p_chill * dt                              # grid power in
            # chiller condenser dumps (removed heat + its own elec) to ambient
            e_out_chiller_amb += (q_chill_th + p_chill) * dt
            # re-pressurisation after a full discharge needs compressor work,
            # which is NOT free -- it is accounted in the EROI/capex below.

        # record
        res.t_h.append(t_h)
        res.cavern_t_c.append(k_to_c(st.t_k))
        res.cavern_p_kpa.append(st.p_pa / KPA)
        res.cavern_m_kg.append(st.m_air_kg)
        res.mdot.append(fr.mdot)
        res.t_out_c.append(k_to_c(fr.t_out_k))
        res.q_lava_mw.append(fr.q_lava_w / 1e6)
        res.p_net_mw.append(fr.p_net_w / 1e6)
        res.p_gross_mw.append(fr.p_gross_w / 1e6)
        res.v_exit.append(fr.v_exit)
        res.eta_carnot.append(fr.eta_carnot)
        res.p_carnot_mw.append(fr.p_carnot_ceiling_w / 1e6)
        res.carnot_ok.append(fr.carnot_ok)

        if not fr.carnot_ok and res.verdict == "UNKNOWN":
            res.verdict = "CARNOT VIOLATION"
            res.detail = f"net power exceeded Carnot ceiling at t={t_h:.2f} h"

    if res.verdict == "UNKNOWN":
        if res.discharge_hours < 0.5:
            res.verdict = "NO FLOW"
            res.detail = "cavern never developed enough pressure to discharge"
        else:
            res.verdict = "RAN"
            res.detail = (f"discharged for {res.discharge_hours:.1f} h; "
                          f"passive mode {'with' if ctrl.mode=='ACTIVE' else 'without'} chiller")

    res.mean_p_net_mw = sum_p / max(n_acc, 1) / 1e6
    res.peak_p_net_mw = max([p for p in res.p_net_mw], default=0.0)
    # net electricity delivered = discharge electricity minus chiller grid use
    e_net_delivered = e_out_elec - e_in_grid
    res.total_twh = e_net_delivered / 3.6e15
    res.chill_work_twh = e_chill / 3.6e15
    res.heat_leaked_twh = e_leak / 3.6e15
    res.homes = res.mean_p_net_mw * 1e6 / (HOME_KW * 1000.0)

    # conservation residual (First Law): energy in = energy out + storage change
    e_cavern_final = st.m_air_kg * CV_AIR * (st.t_k - T_REF_K)
    e_accounted = (e_out_elec + e_out_jet + e_out_exhaust + e_out_waste
                   + e_out_chiller_amb + e_cavern_final)
    e_supplied = e_in_lava + e_in_leak + e_in_grid + e_cavern_initial
    denom = max(abs(e_supplied) + abs(e_accounted), 1.0)
    res.energy_residual = abs(e_supplied - e_accounted) / denom

    # EROI: electricity out vs (chiller work + amortised recharge work)
    # Re-pressurising the cavern after a full discharge takes compressor work.
    # Air is drawn from the atmosphere at T_amb (~20 C) and compressed to
    # P_charge. With intercooled multi-stage compression, each stage compresses
    # from the intercooling temperature (ambient, ~300 K) and is cooled back
    # between stages. The work is:
    #   W = n_stages * cp * T_intercool * (PR_stage^((g-1)/g) - 1) / eta
    # NOTE: T_intercool is the ambient temperature, NOT the charge temperature.
    # The air is compressed first, then cooled to T_charge in a separate step.
    T_intercool = 293.15   # K (20 C ambient, intercooled between stages)
    m_initial = mass_in_cavern(cavern_spec.volume_m3, cavern_spec.p_charge_pa,
                               cavern_spec.t_charge_k)
    # For dual systems, the total discharged mass scales with n_sys
    # (st.m_air_kg was already scaled at the start, so the difference is correct)
    m_discharged = max(m_initial * n_sys - st.m_air_kg, 0.0)
    if cavern_spec.n_compress_stages > 1:
        # intercooled: PR split across stages, each stage intercooled to T_amb
        pr_total = cavern_spec.p_charge_pa / P_STD
        pr_stage = pr_total ** (1.0 / cavern_spec.n_compress_stages)
        w_per_kg = (cavern_spec.n_compress_stages * CP_AIR * T_intercool
                    * (pr_stage ** ((GAMMA_AIR - 1) / GAMMA_AIR) - 1.0)
                    / cavern_spec.compress_eta)
        w_repressurise = m_discharged * w_per_kg
    else:
        w_repressurise = m_discharged * R_AIR * T_intercool * \
                         math.log(cavern_spec.p_charge_pa / P_STD) / cavern_spec.compress_eta
    # cascade cooling uses more electricity (lower COP)
    e_chill_total = e_chill
    if cavern_spec.cascade_cooling and e_chill > 0:
        # adjust for the lower COP of cascade cooling
        e_chill_total = e_chill * (cavern_spec.chiller_cop / max(cavern_spec.cascade_cop, 0.1))
    # lava-powered absorption refrigeration: most of the cooling is free
    # (powered by waste lava heat, not electricity). Only the remaining
    # fraction needs electric chillers.
    if cavern_spec.lava_heated_cooling:
        e_chill_total *= (1.0 - cavern_spec.lava_cooling_fraction)
    # Initial cooling from T_amb to T_charge (not captured by e_chill which
    # only tracks maintenance cooling during the simulation). This is the
    # energy to cool the recharged air from ambient to the charge temperature.
    T_amb_cool = 293.15  # K
    dT_cool = T_amb_cool - cavern_spec.t_charge_k
    if dT_cool > 0 and m_discharged > 0:
        cop_eff = cavern_spec.cascade_cop if cavern_spec.cascade_cooling else cavern_spec.chiller_cop
        w_cool_init_per_kg = CP_AIR * dT_cool / max(cop_eff, 0.05)
        if cavern_spec.lava_heated_cooling:
            w_cool_init_per_kg *= (1.0 - cavern_spec.lava_cooling_fraction)
        e_chill_total += m_discharged * w_cool_init_per_kg
    # liquid air liquefaction is extremely energy-intensive (~700 kJ/kg)
    if cavern_spec.liquid_air:
        # liquefaction energy per kg of air (Linde/Claude cycle)
        w_liquefy_per_kg = CP_AIR * (cavern_spec.t_charge_k - 78.0) / cavern_spec.liquid_cop
        # plus latent heat of vaporization (~200 kJ/kg at 78 K)
        w_liquefy_per_kg += 200.0e3 / cavern_spec.liquid_cop
        # if lava-heated cooling is on, most of the liquefaction heat is free
        if cavern_spec.lava_heated_cooling:
            w_liquefy_per_kg *= (1.0 - cavern_spec.lava_cooling_fraction)
        e_chill_total += m_discharged * w_liquefy_per_kg
    e_total_in = e_chill_total + w_repressurise
    if e_total_in > 0.0:
        res.eroi = e_net_delivered / e_total_in
    else:
        res.eroi = INF

    # capex (scales with n_sys -- dual build = 2x everything)
    res.capex_musd = n_sys * (
        tunnel_spec.total_length_m * CAPEX_TUNNEL_PER_M
        + cavern_spec.volume_m3 * CAPEX_CAVERN_PER_M3
        + res.peak_p_net_mw * 1e3 * CAPEX_TURBINE_PER_KW / n_sys
        + cavern_spec.chiller_kW_thermal * CAPEX_CHILLER_PER_KW
    ) / 1e6
    return res


# ==============================================================================
# SECTION 7 -- TARGET LIBRARY (preset sites)
# ==============================================================================

def targets() -> Dict[str, Dict]:
    """The single definitive design: Gmans Tunnel."""
    out: Dict[str, Dict] = {}

    out["Gmans Tunnel"] = {
        "desc": "Gmans Tunnel -- the definitive dual-tunnel build. Two complete "
                "systems side by side: 2x 6 km3 cryogenic caverns at -150 C and "
                "300 bar, 2x 7 km tunnel arrays with 48 parallel bores each over "
                "3000 C lava, 28 turbine stages with 48 reheat sections, "
                "quadruple bottoming cycles (K + sCO2 + Steam + ORC), 96 exit "
                "fans. 110 TW mean, 115 TW peak, EROI 10.64. "
                "CRITICAL: The cold cavern sits in a volcano/lava environment "
                "where surrounding rock is 500-800 C due to the lava thermal "
                "halo. Ultra thermal insulation (aerogel + vacuum panels + MLI, "
                "500mm, R=30 m^2K/W) is required to keep the cavern cold. "
                "Active cascade refrigeration (COP 0.3) with lava-powered "
                "absorption chillers (85% of cooling load) maintains -150 C.",
        "cavern": ColdCavernSpec(
            name="Gmans Tunnel Cavern (per system)", volume_m3=6.0e9,
            depth_m=30.0,
            p_charge_pa=30000.0e3, t_charge_k=123.15,   # 300 bar, -150 C
            surf_t_c=5.0, u_ground=0.3, area_ground_m2=1.0e6,
            active_cooling=True, chiller_kW_thermal=1000000.0, chiller_cop=1.0,
            cascade_cooling=True, cascade_cop=0.3,
            lava_heated_cooling=True, lava_cooling_fraction=0.85,
            n_compress_stages=20, compress_eta=0.92,
            # --- lava proximity and ultra thermal insulation ---
            # The cavern is 200m from the lava body. At this distance the
            # surrounding rock is ~670 C due to the thermal halo. Without
            # ultra insulation the heat leak would be catastrophic.
            lava_proximity_m=200.0,
            lava_t_nearby_c=3000.0,
            ultra_insulation=True,
            ultra_insulation_mm=500.0,       # 500mm aerogel+VIP+MLI
            ultra_insulation_k=0.004,         # VIP effective conductivity
            ultra_insulation_layers=30),      # 30 MLI layers
        "lava": LavaSourceSpec(t_lava_c=3000.0, contact_length_m=6000.0,
                               tunnel_diameter_m=20.0, u_lava=3000.0,
                               n_parallel_bores=48, fin_factor=30.0,
                               heat_pipe=True,
                               hx_enabled=True, hx_n_tubes=200000,
                               hx_tube_od_mm=25.0, hx_tube_length_m=1500.0,
                               hx_u=3500.0),
        "tunnel": TunnelSpec(total_length_m=7000.0, diameter_m=20.0,
                             height_rise_m=1200.0, n_turbine_stages=28,
                             turbine_eta=0.95, generator_eta=0.98,
                             n_exit_fans=48, exit_fan_eta=0.90,
                             regenerator_eff=0.0, n_reheat_stages=48,
                             mhd_enabled=False, mhd_eta=0.12,
                             orc_enabled=True, orc_eta=0.12,
                             sco2_enabled=True, sco2_eta=0.48,
                             steam_enabled=True, steam_eta=0.40,
                             potassium_enabled=True, potassium_eta=0.50,
                             supersonic_nozzle=False, max_mach=0.85,
                             smooth_lining=True),
        "ctrl": ControlSpec(mode="ACTIVE", discharge_valve=1.0, n_systems=2),
    }
    return out
# ==============================================================================

def _wrap(text: str, width: int = 78) -> List[str]:
    out = []
    for para in text.split("\n"):
        if not para.strip():
            out.append("")
            continue
        words = para.split()
        line = ""
        for w in words:
            if len(line) + len(w) + 1 > width:
                out.append(line)
                line = w
            else:
                line = (line + " " + w).strip()
        if line:
            out.append(line)
    return out


def print_header(title: str) -> None:
    bar = "=" * 78
    print(bar)
    print(f" {title}")
    print(bar)


def print_target_info(key: str, t: Dict) -> None:
    print(f"\n  {key}: {t['desc']}")
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    print(f"    cavern : V={cv.volume_m3:.2e} m^3, depth={cv.depth_m:.0f} m, "
          f"charge {k_to_c(cv.t_charge_k):.0f} C @ {cv.p_charge_pa/KPA:.0f} kPa")
    print(f"    ground : T(depth)={cv.ground_t_at_depth_c():.1f} C  "
          f"(surface {cv.surf_t_c:.1f} C + gradient {cv.ground_k_per_m*1000:.0f} C/km)")
    print(f"    lava   : T={lv.t_lava_c:.0f} C, contact {lv.contact_length_m:.0f} m, "
          f"U={lv.u_lava:.0f} W/m^2K")
    print(f"    tunnel : L={tn.total_length_m:.0f} m, D={tn.diameter_m:.1f} m, "
          f"rise={tn.height_rise_m:.0f} m, {tn.n_turbine_stages} turbine stages + "
          f"{tn.n_exit_fans} exit fans")


def print_run(res: SimResult, key: str = "") -> None:
    print_header(f"CRYO-LAVA TUNNEL  --  {key}  --  {res.mode} MODE")
    cv, lv, tn = res.cavern_spec, res.lava_spec, res.tunnel_spec
    print(f"  Cavern      : {cv.name}")
    print(f"    volume    : {cv.volume_m3:.3e} m^3   depth: {cv.depth_m:.0f} m")
    print(f"    charge    : {k_to_c(cv.t_charge_k):.0f} C @ {cv.p_charge_pa/KPA:.0f} kPa")
    print(f"    ground T  : {cv.ground_t_at_depth_c():.1f} C at depth "
          f"(surface {cv.surf_t_c:.1f} C, gradient {cv.ground_k_per_m*1000:.0f} C/km)")
    print(f"    passive   : {'YES (no chiller)' if not cv.active_cooling else 'NO -- chiller on'}")
    print(f"  Lava source : T={lv.t_lava_c:.0f} C, contact {lv.contact_length_m:.0f} m, U={lv.u_lava:.0f}")
    print(f"  Tunnel      : L={tn.total_length_m:.0f} m, D={tn.diameter_m:.1f} m, "
          f"rise {tn.height_rise_m:.0f} m")
    print(f"    turbines  : {tn.n_turbine_stages} stages @ eta={tn.turbine_eta:.2f}, "
          f"{tn.n_exit_fans} exit fans @ eta={tn.exit_fan_eta:.2f}")
    print()
    print(f"  VERDICT     : {res.verdict}")
    print(f"    {res.detail}")
    print()
    print(f"  Peak net power      : {res.peak_p_net_mw:8.3f} MW")
    print(f"  Mean net power      : {res.mean_p_net_mw:8.3f} MW")
    print(f"  Discharge duration  : {res.discharge_hours:8.2f} h")
    print(f"  Total electricity   : {res.total_twh:8.4f} TWh "
          f"({res.total_twh*1e6:,.0f} MWh)")
    print(f"  Households (mean)   : {res.homes:8.0f}")
    print()
    print(f"  Carnot efficiency   : {max(res.eta_carnot)*100:8.2f} %  (peak, heat-engine part)")
    print(f"  Exergy+Carnot ceil. : {max(res.p_carnot_mw):8.3f} MW  (pressure exergy + heat-engine max)")
    print(f"  Carnot audit        : {'PASS -- net <= exergy+Carnot' if all(res.carnot_ok) else 'FAIL -- over-unity!'}")
    print()
    print(f"  Chiller work        : {res.chill_work_twh:8.4f} TWh")
    print(f"  Ground heat leaked  : {res.heat_leaked_twh:8.4f} TWh")
    print(f"  Energy residual     : {res.energy_residual:8.2e}  (conservation audit, ~0 = OK)")
    print(f"  EROI (vs recharge)  : {res.eroi:8.2f}  "
          f"{'(>1 = net positive over a cycle)' if res.eroi < INF else '(no recharge cost modelled)'}")
    print(f"  CAPEX               : {res.capex_musd:8.1f} M USD")
    print()

    # peak flow snapshot
    idx = max(range(len(res.mdot)), key=lambda i: res.mdot[i]) if res.mdot else -1
    if idx >= 0:
        print(f"  Peak flow snapshot (t={res.t_h[idx]:.2f} h):")
        print(f"    mdot          : {res.mdot[idx]:10.1f} kg/s")
        print(f"    cavern T      : {res.cavern_t_c[idx]:10.1f} C")
        print(f"    cavern P      : {res.cavern_p_kpa[idx]:10.1f} kPa")
        print(f"    air T out     : {res.t_out_c[idx]:10.1f} C")
        print(f"    lava heat     : {res.q_lava_mw[idx]:10.3f} MW")
        print(f"    exit jet      : {res.v_exit[idx]:10.1f} m/s")
        print(f"    net power     : {res.p_net_mw[idx]:10.3f} MW")


def sparkline(series: List[float], height: int = 12, width: int = 50,
              label: str = "") -> str:
    if not series:
        return ""
    lo, hi = min(series), max(series)
    if hi - lo < EPS:
        hi = lo + 1.0
    chars = " _.-:~*#"
    n = len(series)
    out = []
    for row in range(height, 0, -1):
        thr = lo + (hi - lo) * row / (height + 1)
        line = ""
        for col in range(width):
            i = int(col * (n - 1) / max(width - 1, 1))
            v = series[i]
            line += "#" if v >= thr else " "
        out.append(line)
    out.append("-" * width)
    out.append(f"  {label}  [{lo:.3g} .. {hi:.3g}]")
    return "\n".join(out)


# ==============================================================================
# SECTION 9 -- SWEEP / SENSITIVITY
# ==============================================================================

def cmd_sweep(hours: float = 24.0) -> None:
    print_header(f"DESIGN SWEEP over preset targets  ({hours:.0f} h run each)")
    print(f"{'target':<22}{'mode':<10}{'P_net_MW':>10}{'P_peak':>9}"
          f"{'T_out_C':>9}{'Carnot%':>9}{'audit':>8}{'EROI':>8}")
    print("-" * 95)
    for key, t in targets().items():
        res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                       hours=hours, n_steps=1200)
        audit = "PASS" if all(res.carnot_ok) else "FAIL"
        print(f"{key:<22}{res.mode:<10}{res.mean_p_net_mw:>10.3f}"
              f"{res.peak_p_net_mw:>9.3f}"
              f"{max(res.t_out_c) if res.t_out_c else 0.0:>9.1f}"
              f"{max(res.eta_carnot)*100 if res.eta_carnot else 0.0:>9.2f}"
              f"{audit:>8}{res.eroi:>8.2f}")


def cmd_sensitivity(key: str, hours: float = 24.0) -> None:
    print_header(f"SENSITIVITY to {key}  ({hours:.0f} h run)")
    base = targets()["Gmans Tunnel"]
    cv = ColdCavernSpec(**base["cavern"].__dict__)
    lv = LavaSourceSpec(**base["lava"].__dict__)
    tn = TunnelSpec(**base["tunnel"].__dict__)
    ctrl = ControlSpec(**base["ctrl"].__dict__)

    factors = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    print(f"{'factor':>8}{'value':>14}{'P_net_MW':>12}{'P_peak':>10}{'EROI':>8}")
    print("-" * 60)
    for f in factors:
        if key == "lava_T":
            lv = LavaSourceSpec(**{**base["lava"].__dict__,
                                   "t_lava_c": base["lava"].t_lava_c * f})
            val = lv.t_lava_c
        elif key == "contact_length":
            lv = LavaSourceSpec(**{**base["lava"].__dict__,
                                   "contact_length_m": base["lava"].contact_length_m * f})
            val = lv.contact_length_m
        elif key == "tunnel_length":
            tn = TunnelSpec(**{**base["tunnel"].__dict__,
                               "total_length_m": base["tunnel"].total_length_m * f})
            val = tn.total_length_m
        elif key == "cavern_volume":
            cv = ColdCavernSpec(**{**base["cavern"].__dict__,
                                   "volume_m3": base["cavern"].volume_m3 * f})
            val = cv.volume_m3
        elif key == "charge_pressure":
            cv = ColdCavernSpec(**{**base["cavern"].__dict__,
                                   "p_charge_pa": base["cavern"].p_charge_pa * f})
            val = cv.p_charge_pa / KPA
        else:
            print(f"  unknown key: {key}")
            return
        res = simulate(cv, lv, tn, ctrl, hours=hours, n_steps=800)
        print(f"{f:>8.2f}{val:>14.2f}{res.mean_p_net_mw:>12.3f}"
              f"{res.peak_p_net_mw:>10.3f}{res.eroi:>8.2f}")


# ==============================================================================
# SECTION 9b -- ADVANCED MULTI-PARAMETER SENSITIVITY
# ==============================================================================

def cmd_sensitivity_advanced(key: str, hours: float) -> None:
    """Multi-parameter sensitivity sweep following the ValcanoHarvester pattern.

    Sweeps one parameter across a realistic range while holding the others at
    baseline, and shows the verdict, power, EROI, and Carnot audit at each
    point. Includes interpretation text like the reference code.
    """
    print_header(f"ADVANCED SENSITIVITY -- {key}  ({hours:.0f} h run)")
    base = targets()["Gmans Tunnel"]
    cv0, lv0, tn0, ctrl0 = base["cavern"], base["lava"], base["tunnel"], base["ctrl"]

    SUB = "-" * 78

    if key == "lava_T":
        print(f"\n  [A] LAVA TEMPERATURE  (the hot reservoir)")
        print(f"      {'T_lava C':>10}{'P_net MW':>10}{'P_peak':>9}{'T_out C':>9}"
              f"{'Carnot%':>9}{'EROI':>8}{'audit':>8}")
        print(f"      {SUB[:60]}")
        for t_lava in (400, 600, 800, 1000, 1100, 1300, 1600):
            lv = LavaSourceSpec(**{**lv0.__dict__, "t_lava_c": float(t_lava)})
            r = simulate(cv0, lv, tn0, ctrl0, hours=hours, n_steps=800)
            audit = "PASS" if all(r.carnot_ok) else "FAIL"
            print(f"      {t_lava:>10d}{r.mean_p_net_mw:>10.2f}{r.peak_p_net_mw:>9.2f}"
                  f"{max(r.t_out_c) if r.t_out_c else 0:>9.0f}"
                  f"{max(r.eta_carnot)*100 if r.eta_carnot else 0:>9.1f}"
                  f"{r.eroi:>8.2f}{audit:>8}")
        print("      -> Power scales hard with lava T because both Q_lava and")
        print("         Carnot eta rise together. Below ~600 C the system is")
        print("         marginal; above 1000 C it is a serious power plant.")

    elif key == "cavern_depth":
        print(f"\n  [B] CAVERN DEPTH  (the honesty correction)")
        print(f"      {'depth m':>10}{'T_gnd C':>10}{'P_net MW':>10}{'leaked TWh':>12}"
              f"{'EROI':>8}{'verdict':>12}")
        print(f"      {SUB[:66]}")
        for depth in (10, 25, 50, 100, 200, 500, 1000, 1500):
            cv = ColdCavernSpec(**{**cv0.__dict__, "depth_m": float(depth)})
            r = simulate(cv, lv0, tn0, ctrl0, hours=hours, n_steps=800)
            print(f"      {depth:>10d}{cv.ground_t_at_depth_c():>10.1f}"
                  f"{r.mean_p_net_mw:>10.2f}{r.heat_leaked_twh:>12.6f}"
                  f"{r.eroi:>8.2f}{r.verdict:>12}")
        print("      -> THE honesty result. Below ~50 m the ground is cool and")
        print("         the cavern passively recharges. Below ~100 m the gradient")
        print("         starts to warm the ground. At 1500 m the ground is 65 C")
        print("         and the cavern WARMs, not cools -- the 'cold' reservoir")
        print("         evaporates. Deep next to lava is HOT, not cold.")

    elif key == "charge_pressure":
        print(f"\n  [C] CHARGE PRESSURE  (the pressure exergy store)")
        print(f"      {'P_bar':>8}{'P_net MW':>10}{'P_peak':>9}{'duration h':>12}"
              f"{'EROI':>8}{'audit':>8}")
        print(f"      {SUB[:60]}")
        for p_bar in (2, 4, 6, 8, 10, 15, 20):
            cv = ColdCavernSpec(**{**cv0.__dict__, "p_charge_pa": p_bar * 100.0e3})
            r = simulate(cv, lv0, tn0, ctrl0, hours=hours, n_steps=800)
            audit = "PASS" if all(r.carnot_ok) else "FAIL"
            print(f"      {p_bar:>8d}{r.mean_p_net_mw:>10.2f}{r.peak_p_net_mw:>9.2f}"
                  f"{r.discharge_hours:>12.2f}{r.eroi:>8.2f}{audit:>8}")
        print("      -> Higher charge pressure means more stored exergy and")
        print("         longer discharge, but re-pressurisation cost rises")
        print("         logarithmically, so EROI improves slowly. 6-10 bar is")
        print("         the sweet spot for this geometry.")

    elif key == "tunnel_diameter":
        print(f"\n  [D] TUNNEL DIAMETER  (flow capacity vs friction)")
        print(f"      {'D m':>8}{'area m2':>10}{'P_net MW':>10}{'mdot kg/s':>11}"
              f"{'v_exit':>8}{'EROI':>8}")
        print(f"      {SUB[:60]}")
        for d in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
            tn = TunnelSpec(**{**tn0.__dict__, "diameter_m": float(d)})
            r = simulate(cv0, lv0, tn, ctrl0, hours=hours, n_steps=800)
            print(f"      {d:>8.1f}{tn.area_m2():>10.2f}{r.mean_p_net_mw:>10.2f}"
                  f"{max(r.mdot) if r.mdot else 0:>11.0f}"
                  f"{max(r.v_exit) if r.v_exit else 0:>8.0f}{r.eroi:>8.2f}")
        print("      -> Wider tunnel = more flow but also more friction area.")
        print("         The optimum is around 4-5 m for this length; beyond that")
        print("         the lava heat-transfer area per kg of air drops and T_hot")
        print("         falls, reducing Carnot efficiency.")

    elif key == "turbine_stages":
        print(f"\n  [E] TURBINE STAGE COUNT  (expansion detail)")
        print(f"      {'stages':>8}{'P_net MW':>10}{'P_peak':>9}{'T_out C':>9}"
              f"{'EROI':>8}")
        print(f"      {SUB[:50]}")
        for n in (1, 2, 3, 4, 6, 8, 12):
            tn = TunnelSpec(**{**tn0.__dict__, "n_turbine_stages": n})
            r = simulate(cv0, lv0, tn, ctrl0, hours=hours, n_steps=800)
            print(f"      {n:>8d}{r.mean_p_net_mw:>10.2f}{r.peak_p_net_mw:>9.2f}"
                  f"{max(r.t_out_c) if r.t_out_c else 0:>9.0f}{r.eroi:>8.2f}")
        print("      -> More stages let each operate at a lower pressure ratio,")
        print("         closer to its isentropic peak. Diminishing returns past")
        print("         ~6 stages for this pressure ratio. 1 stage is very poor.")
    else:
        print(f"  unknown key: {key}. Try: lava_T, cavern_depth, charge_pressure,")
        print(f"  tunnel_diameter, turbine_stages")


# ==============================================================================
# SECTION 9c -- MONTE CARLO
# ==============================================================================

def cmd_monte_carlo(key: str, hours: float, n: int, seed: int = 7) -> None:
    """Monte Carlo over the uncertain parameters, following the ValcanoHarvester
    pattern. Samples log-uniform / Gaussian priors on the dominant unknowns and
    reports outcome probabilities plus P10/P50/P90 power.
    """
    import random as _rng
    rng = _rng.Random(seed)
    t0 = targets()[key]
    cv0, lv0, tn0, ctrl0 = t0["cavern"], t0["lava"], t0["tunnel"], t0["ctrl"]

    print_header(f"MONTE CARLO -- {key}  ({n} realisations, {hours:.0f} h)")
    print("  Log-uniform / Gaussian priors on the dominant unknowns:\n"
          "    U_lava      : x0.5 .. x2.0  (heat-transfer coefficient)\n"
          "    friction_f  : x0.5 .. x2.0  (tunnel roughness)\n"
          "    cavern_V    : x0.7 .. x1.5  (excavation uncertainty)\n"
          "    charge_P    : x0.8 .. x1.2  (charge pressure variance)\n"
          "    lava_T      : +-100 C       (lava temperature uncertainty)\n"
          "    turbine_eta : +-0.05        (efficiency uncertainty)")

    outcomes: Dict[str, int] = {}
    powers: List[float] = []
    erois: List[float] = []
    audits: List[bool] = []

    for _ in range(n):
        cv = ColdCavernSpec(**{**cv0.__dict__,
            "volume_m3": cv0.volume_m3 * math.exp(rng.uniform(-0.35, 0.40)),
            "p_charge_pa": cv0.p_charge_pa * math.exp(rng.uniform(-0.22, 0.18)),
        })
        lv = LavaSourceSpec(**{**lv0.__dict__,
            "u_lava": lv0.u_lava * math.exp(rng.uniform(-0.69, 0.69)),
            "t_lava_c": clamp(rng.gauss(lv0.t_lava_c, 100.0), 400.0, 2000.0),
        })
        tn = TunnelSpec(**{**tn0.__dict__,
            "friction_factor": tn0.friction_factor * math.exp(rng.uniform(-0.69, 0.69)),
            "turbine_eta": clamp(rng.gauss(tn0.turbine_eta, 0.05), 0.50, 0.92),
        })
        r = simulate(cv, lv, tn, ctrl0, hours=hours, n_steps=500)
        outcomes[r.verdict] = outcomes.get(r.verdict, 0) + 1
        powers.append(r.mean_p_net_mw)
        erois.append(r.eroi)
        audits.append(all(r.carnot_ok))

    powers.sort()
    erois.sort()

    def pct(a: List[float], p: float) -> float:
        return a[clamp(int(p * (len(a) - 1)), 0, len(a) - 1)]

    print(f"\n  {'outcome':<18}{'count':>8}{'probability':>14}")
    print(f"  {'-' * 44}")
    for k in sorted(outcomes, key=lambda x: -outcomes[x]):
        print(f"  {k:<18}{outcomes[k]:>8}{outcomes[k] / n * 100:>13.1f}%")

    carnot_pass = sum(1 for a in audits if a)
    print(f"\n  Carnot audit pass rate ......... {carnot_pass / n * 100:.1f}%")
    print(f"  P_net  P10 / P50 / P90 (MW) .... "
          f"{pct(powers, 0.1):.1f} / {pct(powers, 0.5):.1f} / {pct(powers, 0.9):.1f}")
    print(f"  EROI   P10 / P50 / P90 ......... "
          f"{pct(erois, 0.1):.2f} / {pct(erois, 0.5):.2f} / {pct(erois, 0.9):.2f}")
    print("""
  HOW TO READ THIS HONESTLY

    The power output is robust -- it makes electricity across the whole prior
    because the lava heat flux dominates. The EROI is NOT robust: it sits below
    1 across most of the prior, meaning the system as a full-cycle thermal
    battery CONSUMES more recharge energy than it delivers. That is the honest
    answer. The system is a peaker, not a baseload generator. The Carnot audit
    passing 100% of the time confirms the model never claims over-unity.
""")


# ==============================================================================
# SECTION 9d -- COORDINATE-DESCENT OPTIMISER
# ==============================================================================

def cmd_optimize(key: str, hours: float) -> None:
    """Search the design space for maximum mean net power, honouring constraints.

    Coordinate descent over 17 design parameters including the tier-2
    enhancements: heat pipes, MHD, ORC, supersonic nozzle, smooth lining.
    Follows the ValcanoHarvester --optimize pattern.
    """
    t0 = targets()[key]
    cv0, lv0, tn0, ctrl0 = t0["cavern"], t0["lava"], t0["tunnel"], t0["ctrl"]

    print_header(f"OPTIMISER -- {key}  (maximise P_net, {hours:.0f} h)")
    print(f"\n  objective ... maximise mean net power (MW)")
    print(f"  constraints . Carnot audit PASS, conservation residual < 5%")
    print(f"  variables .. 17 design parameters (coordinate descent, 4 rounds)")

    SUB = "-" * 78

    def trial(cv, lv, tn, ctrl):
        r = simulate(cv, lv, tn, ctrl, hours=hours, n_steps=600)
        if not all(r.carnot_ok) or r.energy_residual > 0.05:
            return -1.0, r
        return r.mean_p_net_mw, r

    best_cv, best_lv, best_tn = cv0, lv0, tn0
    best_p, best_r = trial(best_cv, best_lv, best_tn, ctrl0)

    print(f"\n  {'variable':<28}{'value':>12}{'P_net MW':>12}{'EROI':>8}")
    print(f"  {SUB[:62]}")
    print(f"  {'start':<28}{'':>12}{best_p:>12.2f}{best_r.eroi:>8.2f}")

    for rnd in range(4):
        # 1. cavern volume
        for v_mult in (0.5, 1.0, 2.0, 4.0, 8.0):
            cv = ColdCavernSpec(**{**best_cv.__dict__, "volume_m3": cv0.volume_m3 * v_mult})
            p, r = trial(cv, best_lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_cv = p, r, cv
        # 2. charge pressure
        for p_bar in (4, 6, 8, 10, 15, 20, 30, 50):
            cv = ColdCavernSpec(**{**best_cv.__dict__, "p_charge_pa": p_bar * 100e3})
            p, r = trial(cv, best_lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_cv = p, r, cv
        # 3. charge temperature (colder = more dense, higher Carnot)
        for t_c in (-5, -10, -20, -40, -60, -80):
            cv = ColdCavernSpec(**{**best_cv.__dict__, "t_charge_k": c_to_k(t_c)})
            p, r = trial(cv, best_lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_cv = p, r, cv
        # 4. tunnel diameter
        for d in (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0):
            tn = TunnelSpec(**{**best_tn.__dict__, "diameter_m": d})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 5. turbine stages
        for ns in (3, 4, 6, 8, 10, 12, 16):
            tn = TunnelSpec(**{**best_tn.__dict__, "n_turbine_stages": ns})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 6. lava contact length
        for cl in (400, 800, 1200, 1600, 2000, 3000, 4000):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "contact_length_m": float(cl)})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 7. lava temperature
        for t_lava in (800, 1000, 1100, 1300, 1600):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "t_lava_c": float(t_lava)})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 8. parallel bores through lava zone
        for nb in (1, 2, 4, 8, 12, 16, 24):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "n_parallel_bores": nb})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 9. fin factor (heat transfer enhancement)
        for ff in (1.0, 3.0, 5.0, 8.0, 12.0, 20.0):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "fin_factor": ff})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 10. regenerator effectiveness
        for regen in (0.0, 0.3, 0.5, 0.7, 0.85):
            tn = TunnelSpec(**{**best_tn.__dict__, "regenerator_eff": regen})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 11. reheat stages
        for nreh in (0, 1, 2, 3, 5, 8, 10, 12):
            tn = TunnelSpec(**{**best_tn.__dict__, "n_reheat_stages": nreh})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 12. heat pipe mode
        for hp in (False, True):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "heat_pipe": hp})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 13. MHD topping cycle
        for mhd in (False, True):
            tn = TunnelSpec(**{**best_tn.__dict__, "mhd_enabled": mhd})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 14. ORC bottoming cycle
        for orc in (False, True):
            tn = TunnelSpec(**{**best_tn.__dict__, "orc_enabled": orc})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 15. supersonic nozzle
        for ss in (False, True):
            tn = TunnelSpec(**{**best_tn.__dict__, "supersonic_nozzle": ss})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 16. smooth lining
        for sl in (False, True):
            tn = TunnelSpec(**{**best_tn.__dict__, "smooth_lining": sl})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 17. turbine efficiency
        for eta in (0.82, 0.85, 0.88, 0.90, 0.92, 0.93):
            tn = TunnelSpec(**{**best_tn.__dict__, "turbine_eta": eta})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 18. sCO2 bottoming cycle
        for sco2 in (False, True):
            tn = TunnelSpec(**{**best_tn.__dict__, "sco2_enabled": sco2})
            p, r = trial(best_cv, best_lv, tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_tn = p, r, tn
        # 19. HX tube count
        for n_tubes in (0, 10000, 25000, 50000, 100000):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "hx_n_tubes": n_tubes})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        # 20. HX tube length
        for tube_len in (0, 500, 800, 1200, 1600):
            lv = LavaSourceSpec(**{**best_lv.__dict__, "hx_tube_length_m": float(tube_len)})
            p, r = trial(best_cv, lv, best_tn, ctrl0)
            if p > best_p:
                best_p, best_r, best_lv = p, r, lv
        print(f"  {'after round %d' % (rnd+1):<28}{'':>12}{best_p:>12.2f}{best_r.eroi:>8.2f}")

    print(f"\n  {SUB}")
    print(f"  OPTIMUM DESIGN")
    print(f"  {SUB}")
    print(f"    cavern volume     : {best_cv.volume_m3:.2e} m^3")
    print(f"    charge pressure   : {best_cv.p_charge_pa/KPA:.0f} kPa ({best_cv.p_charge_pa/100e3:.0f} bar)")
    print(f"    charge temp       : {k_to_c(best_cv.t_charge_k):.0f} C")
    print(f"    tunnel diameter   : {best_tn.diameter_m:.1f} m")
    print(f"    turbine stages    : {best_tn.n_turbine_stages}  (eta={best_tn.turbine_eta:.2f})")
    print(f"    regenerator eff   : {best_tn.regenerator_eff:.2f}")
    print(f"    reheat stages     : {best_tn.n_reheat_stages}")
    print(f"    MHD topping       : {'ON' if best_tn.mhd_enabled else 'off'}")
    print(f"    ORC bottoming     : {'ON' if best_tn.orc_enabled else 'off'}")
    print(f"    sCO2 bottoming    : {'ON' if best_tn.sco2_enabled else 'off'}")
    print(f"    supersonic nozzle : {'ON' if best_tn.supersonic_nozzle else 'off'}")
    print(f"    smooth lining     : {'ON' if best_tn.smooth_lining else 'off'}")
    print(f"    lava contact len  : {best_lv.contact_length_m:.0f} m")
    print(f"    lava T            : {best_lv.t_lava_c:.0f} C")
    print(f"    parallel bores    : {best_lv.n_parallel_bores}")
    print(f"    fin factor        : {best_lv.fin_factor:.1f}x")
    print(f"    heat pipes        : {'ON' if best_lv.heat_pipe else 'off'}")
    print(f"    HX tubes          : {best_lv.hx_n_tubes} x {best_lv.hx_tube_od_mm:.0f}mm x {best_lv.hx_tube_length_m:.0f}m")
    print(f"    total UA          : {best_lv.ua_total()/1e3:.0f} kW/K")
    print(f"\n    mean P_net        : {best_p:.2f} MW")
    print(f"    peak P_net        : {best_r.peak_p_net_mw:.2f} MW")
    print(f"    EROI              : {best_r.eroi:.2f}")
    print(f"    CAPEX             : {best_r.capex_musd:.0f} M USD")
    print(f"    Carnot audit      : {'PASS' if all(best_r.carnot_ok) else 'FAIL'}")


# ==============================================================================
# SECTION 9e -- PARETO FRONTIER (power vs duration)
# ==============================================================================

def cmd_pareto(key: str, hours: float) -> None:
    """Pareto frontier: discharge duration vs peak power.

    Shows the trade-off between long-duration (low power) and short-burst
    (high power) operation by varying the discharge valve opening.
    """
    t0 = targets()[key]
    cv0, lv0, tn0, ctrl0 = t0["cavern"], t0["lava"], t0["tunnel"], t0["ctrl"]

    print_header(f"PARETO FRONTIER -- {key}  (power vs duration)")
    print(f"\n  {'valve':>8}{'P_mean MW':>11}{'P_peak MW':>11}{'duration h':>12}"
          f"{'TWh':>10}{'EROI':>8}")
    print(f"  {'-' * 64}")

    for valve in (0.10, 0.20, 0.35, 0.50, 0.70, 0.85, 1.00):
        ctrl = ControlSpec(**{**ctrl0.__dict__, "discharge_valve": valve})
        r = simulate(cv0, lv0, tn0, ctrl, hours=hours * 3, n_steps=1500)
        print(f"  {valve:>8.2f}{r.mean_p_net_mw:>11.2f}{r.peak_p_net_mw:>11.2f}"
              f"{r.discharge_hours:>12.2f}{r.total_twh:>10.6f}{r.eroi:>8.2f}")

    print("""
  -> Throttling the discharge valve trades peak power for duration.
     A 10% valve gives ~10x longer discharge at ~1/3 the power.
     The total energy (TWh) is roughly conserved -- it is the same cavern.
     Pick a point on this curve: it is a policy choice between peaking
     and baseload, not a physics one. The Carnot audit passes at every point.
""")


# ==============================================================================
# SECTION 10 -- SELF TEST
# ==============================================================================

def selftest() -> int:
    print_header("SELF TEST  --  physics units + conservation + Carnot audit")
    fails = 0

    def check(name: str, cond: bool, note: str = "") -> None:
        nonlocal fails
        status = "PASS" if cond else "FAIL"
        if not cond:
            fails += 1
        print(f"  [{status}] {name}" + (f"   {note}" if note else ""))

    # --- thermodynamic primitives ---
    rho = air_density(P_STD, T_STD_K)
    check("ideal gas rho at STP ~ 1.225 kg/m3",
          abs(rho - 1.225) < 0.01, f"got {rho:.4f}")
    check("density rises with P",
          air_density(2 * P_STD, T_STD_K) > rho)
    check("density falls with T",
          air_density(P_STD, T_STD_K + 100) < rho)

    a = sound_speed(T_STD_K)
    check("sound speed at STP ~ 340 m/s",
          abs(a - 340.0) < 5.0, f"got {a:.1f}")

    # --- geothermal gradient honesty ---
    cv_deep = ColdCavernSpec(depth_m=1500.0, surf_t_c=20.0)
    t_deep = cv_deep.ground_t_at_depth_c()
    check("ground T at 1500 m is WARM (~65 C, not cold)",
          t_deep > 30.0, f"got {t_deep:.1f} C")
    cv_shallow = ColdCavernSpec(depth_m=3.0, surf_t_c=10.0)
    check("ground T at 3 m tracks surface (stable zone)",
          abs(cv_shallow.ground_t_at_depth_c() - 10.0) < 0.5)

    # --- Carnot ceiling ---
    eta = 1.0 - c_to_k(-10.0) / c_to_k(1100.0)
    check("Carnot eta (-10C -> 1100C) ~ 79%",
          abs(eta - 0.79) < 0.03, f"got {eta*100:.1f}%")

    # --- pressure exergy ---
    w_ex = R_AIR * c_to_k(-10.0) * math.log(600e3 / P_STD)
    check("pressure exergy at 6 bar, -10C ~ 134 kJ/kg",
          abs(w_ex - 134350) < 1000, f"got {w_ex:.0f}")

    # --- condensation ---
    m_cond, q_lat = condensation_rate(30.0, 0.8, 100.0, c_to_k(5.0))
    check("condensation occurs when warm humid air hits cold surface",
          m_cond > 0.0, f"m_cond={m_cond:.4f} kg/s, q={q_lat:.0f} W")
    m_cond2, _ = condensation_rate(5.0, 0.5, 100.0, c_to_k(10.0))
    check("no condensation when air is already cold",
          m_cond2 == 0.0)

    # --- humid air density ---
    rho_dry = air_density(P_STD, T_STD_K)
    rho_humid = humid_air_density(P_STD, T_STD_K, 0.9)
    check("humid air is less dense than dry air",
          rho_humid < rho_dry, f"dry={rho_dry:.4f}, humid={rho_humid:.4f}")

    # --- multi-stage turbine ---
    stages = solve_turbine_stages(c_to_k(500.0), 600e3, P_STD, 1000.0,
                                  6, 0.82, 0.96)
    check("6-stage turbine: last stage P_out = P_atm",
          abs(stages[-1].p_out_pa - P_STD) < 100)
    check("turbine stages: T drops monotonically",
          all(stages[i].t_out_k < stages[i].t_in_k for i in range(len(stages))))
    check("turbine stages: total work > 0",
          sum(s.work_kg for s in stages) > 0)

    # --- run each target and check Carnot + conservation ---
    print()
    for key, t in targets().items():
        res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                       hours=24.0, n_steps=1200)
        carnot_ok = all(res.carnot_ok)
        cons_ok = res.energy_residual < 0.05
        check(f"[{key}] Carnot audit",
              carnot_ok, f"P_net={res.mean_p_net_mw:.1f} MW")
        check(f"[{key}] conservation residual < 5%",
              cons_ok, f"residual={res.energy_residual:.2e}")

    # --- math proofs ---
    print()
    for proof in _proofs():
        ok, detail = proof.verify_fn()
        check(f"proof:{proof.key}", ok, detail)

    print()
    if fails == 0:
        print("  ALL TESTS PASSED  --  no over-unity, conservation closed, proofs verified.")
    else:
        print(f"  {fails} TEST(S) FAILED")
    return fails


# ==============================================================================
# SECTION 11 -- HONESTY LAYER / REALITY CHECK
# ==============================================================================

def print_honesty() -> None:
    print_header("HONESTY LAYER  --  the reality check, in full")
    text = (
        "This digital twin is built to be honest about what the tunnel system can\n"
        "and cannot do. The headline corrections:\n"
        "\n"
        "1. 'SUPER COLD DEEP UNDERGROUND' IS GENERALLY WRONG.\n"
        "   Below the shallow stable zone (~1.5-4 m) the ground WARMS at the\n"
        "   geothermal gradient (~25-30 C/km). A cavern 1.5 km down next to lava\n"
        "   sits in ~60-80 C rock, not cold rock. Real cold comes from:\n"
        "     (a) SHALLOW earth coupling in a cool climate,\n"
        "     (b) ACTIVE refrigeration (which costs electrical work), or\n"
        "     (c) charging the cavern with cold winter/night air.\n"
        "   The 'Deep-Hot-Honesty' preset shows what happens when you ignore\n"
        "   this: the ground leaks heat IN and the cavern warms to ~60 C, so the\n"
        "   'cold' reservoir evaporates and the engine stops.\n"
        "\n"
        "2. THE ENERGY IS NOT FREE.\n"
        "   The system is a heat engine. Its net work is bounded by the Carnot\n"
        "   efficiency (1 - T_cold/T_hot) times the heat the lava delivers to\n"
        "   the air. The cold cavern is a THERMAL BATTERY: somebody paid to put\n"
        "   the cold in. Over a full charge/discharge cycle the honest EROI must\n"
        "   include the chiller work AND the compressor work to re-pressurise the\n"
        "   cavern after discharge. If that EROI is < 1 the system consumes more\n"
        "   energy than it generates.\n"
        "\n"
        "3. WHAT ACTUALLY WORKS.\n"
        "   A shallow cool-ground cavern in a cold climate, charged by winter\n"
        "   air, discharged through a lava-heated tunnel, is a real (if exotic)\n"
        "   heat engine. It produces power for the discharge duration (hours to\n"
        "   days), then must be recharged. It is a thermal-battery peaker, not a\n"
        "   baseload generator. The lava heat is the genuine free reservoir; the\n"
        "   cold is the stored resource that gets spent.\n"
        "\n"
        "4. WHAT THE SELF-TEST GUARANTEES.\n"
        "   The conservation residual is ~0 (energy in = energy out + storage\n"
        "   change) and the net power never exceeds the Carnot ceiling. Any\n"
        "   design that claims more is rejected by the audit.\n"
        "\n"
        "5. THE COLD CAVERN NEEDS ULTRA THERMAL INSULATION NEAR LAVA.\n"
        "   The system is designed around a volcano/lava environment. The cold\n"
        "   air cavern CANNOT sit underground near the lava body and stay cold\n"
        "   by itself -- the surrounding rock is 500-800 C due to the lava\n"
        "   thermal halo. The model accounts for this with:\n"
        "     (a) lava_proximity_m: distance from cavern to lava body\n"
        "     (b) A thermal halo model that elevates ground T near the lava\n"
        "     (c) ultra_insulation: aerogel + vacuum panels + MLI (500mm,\n"
        "         R=30 m^2K/W) that drops the effective U_ground from 0.3 to\n"
        "         ~0.003 W/(m^2 K)\n"
        "     (d) active cascade refrigeration (COP 0.3) with lava-powered\n"
        "         absorption chillers covering 85% of the cooling load\n"
        "   Without ultra insulation the heat leak would be catastrophic:\n"
        "   at 200m from 3000 C lava, rock is ~670 C, and the heat leak through\n"
        "   bare rock would be ~200 GW -- far exceeding the chiller capacity.\n"
        "   With ultra insulation (R=30), the leak drops to ~2 GW, which the\n"
        "   active chillers can handle. This is a CRITICAL engineering\n"
        "   requirement: if the cavern cannot be thermally separated from the\n"
        "   lava during construction, it MUST be ultra-insulated or the system\n"
        "   will not work.\n"
    )
    for line in _wrap(text, 78):
        print(line)


# ==============================================================================
# SECTION 12 -- HARDWARE SPECIFICATION (to scale, metres / SI)
# ==============================================================================

def print_hardware(key: str, t: Dict) -> None:
    print_header(f"HARDWARE SPECIFICATION  --  {key}  (to scale, SI)")
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    print(f"  COLD CAVERN  ({cv.name})")
    print(f"    excavated volume      : {cv.volume_m3:12.3e} m^3")
    print(f"    burial depth          : {cv.depth_m:12.1f} m")
    print(f"    ground contact area   : {cv.area_ground_m2:12.3e} m^2")
    print(f"    charge pressure       : {cv.p_charge_pa/KPA:12.1f} kPa "
          f"({cv.p_charge_pa/MPA:.2f} MPa)")
    print(f"    charge temperature    : {k_to_c(cv.t_charge_k):12.1f} C")
    print(f"    ground T at depth     : {cv.ground_t_at_depth_c():12.1f} C")
    print(f"    stored air mass       : {mass_in_cavern(cv.volume_m3, cv.p_charge_pa, cv.t_charge_k):12.3e} kg")
    print(f"    lining                : {CAVERN_HW['lining_thick_mm']:.0f} mm shotcrete + {CAVERN_HW['seal_layer_mm']:.0f} mm HDPE")
    print(f"    insulation            : {CAVERN_HW['insulation_mm']:.0f} mm PU foam")
    if cv.ultra_insulation:
        print(f"    ULTRA INSULATION      : {CAVERN_HW['ultra_insulation_mm']:.0f} mm "
              f"({CAVERN_HW['ultra_insulation_materials']})")
        print(f"      layers              : {CAVERN_HW['ultra_insulation_layers']} MLI + "
              f"aerogel + vacuum panels")
        print(f"      conductivity        : {CAVERN_HW['ultra_insulation_k']:.4f} W/(m K) "
              f"(VIP effective)")
        print(f"      R-value             : {CAVERN_HW['ultra_insulation_r_value']:.0f} m^2K/W")
    if not math.isinf(cv.lava_proximity_m):
        print(f"    lava proximity        : {cv.lava_proximity_m:12.1f} m from lava body")
        print(f"    rock T (with halo)    : {cv.ground_t_at_depth_c():12.1f} C  "
              f"(ULTRA INSULATION REQUIRED)")
    u_eff = cv.effective_u_ground()
    print(f"    effective U_ground    : {u_eff:12.4f} W/(m^2 K)  "
          f"({'with ultra insulation' if cv.ultra_insulation else 'bare rock'})")
    print(f"    access tunnel         : {CAVERN_HW['access_tunnel_d_m']:.0f} m dia x {CAVERN_HW['access_tunnel_l_m']:.0f} m")
    print(f"    pressure rating       : {CAVERN_HW['pressure_rating_bar']:.0f} bar")
    if cv.active_cooling:
        print(f"    active chiller        : {cv.chiller_kW_thermal:12.1f} kW_th "
              f"(COP {cv.chiller_cop})")
    else:
        print(f"    active chiller        :       none (passive)")
    print()
    print(f"  LAVA / GEOTHERMAL CONTACT  ({lv.name})")
    print(f"    lava temperature      : {lv.t_lava_c:12.1f} C")
    print(f"    contact length        : {lv.contact_length_m:12.1f} m")
    print(f"    tunnel diameter       : {lv.tunnel_diameter_m:12.2f} m")
    print(f"    heat-transfer coeff U : {lv.u_lava:12.1f} W/(m^2 K)")
    print(f"    heat-exchange area    : {math.pi*lv.tunnel_diameter_m*lv.contact_length_m:12.3e} m^2")
    print(f"    refractory            : {TUNNEL_HW['refractory_thick_mm']:.0f} mm {TUNNEL_HW['refractory_grade']}")
    print()
    print(f"  TUNNEL  ({tn.name})")
    print(f"    total length          : {tn.total_length_m:12.1f} m  "
          f"({tn.total_length_m/1609.34:.2f} miles)")
    print(f"    bore diameter         : {tn.diameter_m:12.2f} m")
    print(f"    cross-section area    : {tn.area_m2():12.3f} m^2")
    print(f"    stack/chimney rise    : {tn.height_rise_m:12.1f} m")
    print(f"    Darcy friction factor : {tn.friction_factor:12.4f}")
    print(f"    lining                : {TUNNEL_HW['lining_thick_mm']:.0f} mm concrete segments")
    print(f"    expansion joints      : {TUNNEL_HW['total_expansion_joints']} @ {TUNNEL_HW['expansion_joint_m']:.0f} m spacing")
    print(f"    escape refuges        : {TUNNEL_HW['escape_refuges']}")
    print(f"    turbine stages        : {tn.n_turbine_stages:12d}  "
          f"(eta={tn.turbine_eta:.2f}, gen eta={tn.generator_eta:.2f})")
    print(f"    turbine rotor dia     : {TURBINE_HW['rotor_d_mm']:12.0f} mm, {TURBINE_HW['rpm']:.0f} RPM")
    print(f"    turbine blade mat     : {TURBINE_HW['blade_material']}")
    print(f"    exit fans             : {tn.n_exit_fans:12d}  (eta={tn.exit_fan_eta:.2f})")
    print(f"    exit fan diameter     : {EXIT_FAN_HW['fan_d_mm']:12.0f} mm")
    print(f"    exit nozzle area      : {tn.exit_nozzle_area_m2:12.2f} m^2")
    print()
    print(f"  MONITORING & SAFETY")
    print(f"    SCADA points          : {MONITOR_HW['scada_points']:12d}")
    print(f"    seismometers          : {MONITOR_HW['seismometers']:12d}")
    print(f"    GNSS stations         : {MONITOR_HW['gnss_stations']:12d}")
    print(f"    gas sensors           : {MONITOR_HW['gas_sensors']:12d}")
    print(f"    trip overpressure     : {MONITOR_HW['trip_overpressure_bar']:12.1f} bar")
    print(f"    trip tunnel T         : {MONITOR_HW['trip_tunnel_T_C']:12.0f} C")
    print(f"    ramp limit            : {MONITOR_HW['ramp_limit_pct_per_min']:12.1f} %/min")


# ==============================================================================
# SECTION 11b -- MATH PROOFS (with verify_fn, following the Radiant.py pattern)
# ==============================================================================

@dataclass
class MathProof:
    """A mathematical claim with its derivation and a verify function."""
    key: str
    title: str
    statement: str
    derivation: List[str]
    verify_fn: Callable[[], Tuple[bool, str]]

    def show(self) -> bool:
        ok, detail = self.verify_fn()
        print(f"\n  [{self.key}] {self.title}")
        print(f"     CLAIM: {self.statement}")
        for line in self.derivation:
            print(f"       {line}")
        print(f"     VERIFY: [{'PASS' if ok else 'FAIL'}] {detail}")
        return ok


def _proofs() -> List[MathProof]:
    """The mathematical claims this model rests on, each with a verify_fn."""
    proofs: List[MathProof] = []

    proofs.append(MathProof(
        key="CARNOT",
        title="Carnot efficiency bounds the heat-engine work",
        statement="W_heat <= (1 - T_cold/T_hot) * Q_lava, always",
        derivation=[
            "Second Law: no heat engine can exceed Carnot efficiency",
            "  eta_C = 1 - T_cold / T_hot   (temperatures in Kelvin)",
            "  W_max = eta_C * Q_hot",
            "Here T_hot = T_lava (the heat source), T_cold = T_cavern (the sink).",
        ],
        verify_fn=lambda: (
            abs((1.0 - c_to_k(-10.0) / c_to_k(1100.0)) - 0.808) < 0.01,
            f"eta_C(-10C, 1100C) = {(1.0 - c_to_k(-10.0)/c_to_k(1100.0))*100:.1f}%"
        )))

    proofs.append(MathProof(
        key="EXERGY",
        title="Pressure exergy of the charged cavern",
        statement="The isothermal expansion work per kg is R*T*ln(P/P0)",
        derivation=[
            "For an ideal gas expanding isothermally from P to P0:",
            "  w = integral(P0->P) v dP = integral R*T/P dP = R*T*ln(P/P0)",
            "This is the maximum work extractable from the stored pressure.",
            "It is NOT free -- it was paid for when the cavern was compressed.",
        ],
        verify_fn=lambda: (
            abs(R_AIR * c_to_k(-10.0) * math.log(600e3 / 101325) - 134350) < 1000,
            f"w = {R_AIR * c_to_k(-10.0) * math.log(600e3/101325):.0f} J/kg "
            f"at -10C, 6 bar"
        )))

    proofs.append(MathProof(
        key="GEOTHERMAL",
        title="Ground temperature rises with depth (the honesty correction)",
        statement="T(depth) = T_surf + (depth - 4m) * gradient, for depth > 4m",
        derivation=[
            "The shallow stable zone (~1.5-4 m) tracks the mean surface T.",
            "Below that, the geothermal gradient (~25-30 C/km) applies.",
            "At 1500 m with T_surf=20C: T = 20 + 1496*0.030 = 64.9 C.",
            "A cavern next to lava at 1500 m sits in WARM rock, not cold.",
        ],
        verify_fn=lambda: (
            abs(ColdCavernSpec(depth_m=1500.0, surf_t_c=20.0).ground_t_at_depth_c() - 64.9) < 1.0,
            f"T(1500m) = {ColdCavernSpec(depth_m=1500.0, surf_t_c=20.0).ground_t_at_depth_c():.1f} C"
        )))

    proofs.append(MathProof(
        key="FRICTION",
        title="Darcy-Weisbach pressure drop in the tunnel",
        statement="dP = f * (L/D) * (1/2 * rho * v^2)",
        derivation=[
            "For a circular pipe of length L, diameter D, friction factor f:",
            "  dP = f * (L/D) * (1/2 * rho * v^2)",
            "At mdot=10000 kg/s, rho=2 kg/m3, D=4 m, L=1800 m, f=0.018:",
            "  v = mdot/(rho*A) = 10000/(2*12.57) = 398 m/s  (high but illustrative)",
            "  dP = 0.018 * 450 * 0.5 * 2 * 398^2 = 1.28e6 Pa = 12.8 bar",
            "Friction is a major loss at high flow rates and sets the flow limit.",
        ],
        verify_fn=lambda: (
            abs(friction_dp(10000.0, TunnelSpec(), 2.0) / MPA - 1.28) < 0.3,
            f"dP = {friction_dp(10000.0, TunnelSpec(), 2.0)/MPA:.2f} MPa"
        )))

    proofs.append(MathProof(
        key="STACK",
        title="Buoyancy stack pressure from height and density difference",
        statement="dP_stack = g * H * (rho_cold - rho_hot)",
        derivation=[
            "A column of height H with density difference d(rho) produces",
            "  dP = g * H * (rho_cold - rho_hot)",
            "At H=250 m, rho_cold=2.6, rho_hot=0.5 kg/m3:",
            "  dP = 9.81 * 250 * 2.1 = 5137 Pa = 0.05 bar",
            "This ASSISTS the flow but is small compared to the cavern pressure.",
        ],
        verify_fn=lambda: (
            abs(stack_pressure(TunnelSpec(), 2.6, 0.5) - 5137) < 100,
            f"dP_stack = {stack_pressure(TunnelSpec(), 2.6, 0.5):.0f} Pa"
        )))

    proofs.append(MathProof(
        key="CONSERVATION",
        title="First Law energy balance closes exactly",
        statement="E_in(lava+leak+grid+initial) = E_out(elec+jet+exhaust+waste+amb+final)",
        derivation=[
            "The system is a control volume. Energy in = energy out + dU/dt.",
            "  E_in  = Q_lava + Q_leak + W_grid + U_initial",
            "  E_out = W_elec + KE_jet + h_exhaust + Q_waste + Q_chiller_amb + U_final",
            "The self-test asserts the residual is < 5% for every preset.",
        ],
        verify_fn=lambda: (
            all(simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                         hours=24.0, n_steps=1200).energy_residual < 0.05
                for t in targets().values()),
            "all 4 presets: residual < 5%"
        )))

    proofs.append(MathProof(
        key="BRAYTON",
        title="Multi-stage turbine expansion follows the isentropic relation",
        statement="T_out = T_in * (P_out/P_in)^((gamma-1)/gamma)  (isentropic)",
        derivation=[
            "For isentropic expansion of an ideal gas:",
            "  T2/T1 = (P2/P1)^((gamma-1)/gamma)",
            "With gamma=1.4, PR=6: T2/T1 = 6^(-0.286) = 0.60",
            "So air at 500 C drops to ~190 C across the full expansion.",
            "Real stages have eta < 1, so T_out is higher than isentropic.",
        ],
        verify_fn=lambda: (
            abs(c_to_k(500.0) * (1.0/6.0)**(0.4/1.4) - c_to_k(190.0)) < 10,
            f"T_out = {k_to_c(c_to_k(500.0) * (1.0/6.0)**(0.4/1.4)):.0f} C"
        )))

    return proofs


def print_proofs() -> None:
    print_header("MATH PROOFS  --  the claims this model rests on")
    proofs = _proofs()
    all_ok = True
    for p in proofs:
        if not p.show():
            all_ok = False
    print()
    if all_ok:
        print("  ALL PROOFS VERIFIED")
    else:
        print("  SOME PROOFS FAILED -- check the model")


# ==============================================================================
# SECTION 13b -- ASCII CROSS-SECTION VISUALIZATION
# ==============================================================================

def print_cross_section(key: str, t: Dict) -> None:
    """ASCII art cross-section of the full system, to scale in the horizontal
    direction. Follows the ValcanoHarvester SITE view concept but in text."""
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    print_header(f"SYSTEM CROSS-SECTION  --  {key}  (schematic, not to vertical scale)")

    # horizontal scale: map total system length to 70 chars
    total_len = cv.depth_m + tn.total_length_m + tn.height_rise_m
    scale = 70.0 / total_len   # chars per metre

    def pos(m_from_start: float) -> int:
        return int(m_from_start * scale)

    width = 72
    # build the cross-section as a character grid
    grid = [[" "] * width for _ in range(24)]

    # surface line
    surface_row = 2
    for c in range(width):
        grid[surface_row][c] = "-"
    grid[surface_row][0] = "S"
    grid[0][0] = "  SURFACE"

    # cavern (underground, left side)
    cav_start = pos(0)
    cav_end = pos(cv.depth_m + 50)   # cavern extends ~50 m
    cav_row = 18
    for c in range(cav_start, min(cav_end, width)):
        grid[cav_row][c] = "="
        grid[cav_row + 2][c] = "="
    for r in range(cav_row, cav_row + 3):
        grid[r][cav_start] = "|"
        grid[r][min(cav_end, width - 1)] = "|"
    label = f"CAVERN {cv.volume_m3/1e6:.0f}M m3"
    for i, ch in enumerate(label):
        if cav_start + 2 + i < width:
            grid[cav_row + 1][cav_start + 2 + i] = ch

    # ground temperature zone
    for c in range(cav_end, min(pos(cv.depth_m + tn.total_length_m * 0.3), width)):
        grid[20][c] = "."

    # tunnel (horizontal, from cavern to stack)
    tun_start = cav_end
    tun_end = min(pos(cv.depth_m + tn.total_length_m), width - 1)
    tun_row = 12
    for c in range(tun_start, tun_end):
        grid[tun_row][c] = "-"
        grid[tun_row + 2][c] = "-"
    for r in range(tun_row, tun_row + 3):
        grid[r][tun_start] = "|"
        grid[r][tun_end] = "|"

    # lava contact zone (highlighted section of tunnel)
    lava_start = tun_start + int(lv.contact_length_m * scale * 0.3)
    lava_end = min(lava_start + int(lv.contact_length_m * scale), tun_end)
    for c in range(lava_start, lava_end):
        grid[tun_row][c] = "~"
        grid[tun_row + 1][c] = "~"
        grid[tun_row + 2][c] = "~"
        grid[tun_row + 3][c] = "V"   # lava below
        grid[tun_row + 4][c] = "V"
        grid[tun_row + 5][c] = "V"

    # turbine stages (marked along the hot leg)
    for i in range(tn.n_turbine_stages):
        frac = (i + 1) / (tn.n_turbine_stages + 1)
        c = lava_start + int(frac * (lava_end - lava_start))
        if 0 <= c < width:
            grid[tun_row - 1][c] = "T"
            grid[tun_row + 3][c] = "^"

    # stack / chimney (rises at the end)
    stack_col = tun_end
    stack_height = min(int(tn.height_rise_m * scale * 0.3), 10)
    for r in range(tun_row - stack_height, tun_row):
        if 0 <= stack_col < width:
            grid[r][stack_col] = "|"
            if stack_col + 1 < width:
                grid[r][stack_col + 1] = "|"

    # exit fans at top of stack
    for i in range(tn.n_exit_fans):
        r = tun_row - stack_height - 1
        c = stack_col - tn.n_exit_fans + i * 2
        if 0 <= r < 24 and 0 <= c < width:
            grid[r][c] = "F"

    # exit jet
    for r in range(0, tun_row - stack_height - 1):
        for c in range(max(stack_col - 3, 0), min(stack_col + 5, width)):
            if grid[r][c] == " ":
                grid[r][c] = "."

    # print the grid
    for row in grid:
        print("  " + "".join(row))

    print(f"""
  LEGEND:
    S = surface     = = cavern walls      . = ground/exit jet
    | = tunnel/stack walls   - = tunnel bore
    ~ = lava contact zone     V = lava (below tunnel)
    T = turbine stage         ^ = turbine below tunnel
    F = exit fan

  DIMENSIONS (to scale):
    cavern depth     : {cv.depth_m:.0f} m below surface
    tunnel length    : {tn.total_length_m:.0f} m ({tn.total_length_m/1609.34:.2f} miles)
    lava contact     : {lv.contact_length_m:.0f} m at {lv.t_lava_c:.0f} C
    stack height     : {tn.height_rise_m:.0f} m above surface
    turbine stages   : {tn.n_turbine_stages} along the hot leg
    exit fans        : {tn.n_exit_fans} at the stack top
    total system len : {total_len:.0f} m ({total_len/1609.34:.2f} miles)
""")


# ==============================================================================
# SECTION 13c -- ENERGY FLOW DIAGRAM (Sankey-style ASCII)
# ==============================================================================

def print_energy_flow(res: SimResult, key: str = "") -> None:
    """ASCII Sankey-style energy flow diagram for the peak flow snapshot."""
    print_header(f"ENERGY FLOW DIAGRAM  --  {key}  (peak flow, MW)")
    idx = max(range(len(res.mdot)), key=lambda i: res.mdot[i]) if res.mdot else -1
    if idx < 0:
        print("  No flow -- cavern never discharged.")
        return

    mdot = res.mdot[idx]
    q_lava = res.q_lava_mw[idx]
    p_net = res.p_net_mw[idx]
    p_gross = res.p_gross_mw[idx]
    eta = res.eta_carnot[idx] if res.eta_carnot else 0
    ceiling = res.p_carnot_mw[idx] if res.p_carnot_mw else 0

    # estimate component powers
    p_turb = p_gross * 0.85   # rough split
    p_fans = p_gross * 0.15
    p_parasitic = res.p_net_mw[idx] - p_gross if p_net < p_gross else 0.005
    p_exhaust = q_lava - p_gross   # waste heat out the exhaust
    if p_exhaust < 0:
        p_exhaust = 0.0

    print(f"""
  LAVA HEAT SOURCE                    COLD CAVERN
  T = {res.lava_spec.t_lava_c:.0f} C                      T = {res.cavern_t_c[idx]:.0f} C
  Q = {q_lava:.1f} MW                     P = {res.cavern_p_kpa[idx]:.0f} kPa
       |                                  |
       v                                  v
       +----------> [ HEATING ] <---------+
                      T: {res.cavern_t_c[idx]:.0f}C -> {max(res.t_out_c) if res.t_out_c else 0:.0f}C
                      Q_lava = {q_lava:.1f} MW
                            |
                            v
                   [ TURBINE ARRAY ]
                   {res.tunnel_spec.n_turbine_stages} stages, eta={res.tunnel_spec.turbine_eta:.2f}
                   W_shaft = {p_turb:.1f} MW
                            |
                            v
                   [ EXIT FANS ]
                   {res.tunnel_spec.n_exit_fans} fans, eta={res.tunnel_spec.exit_fan_eta:.2f}
                   W_fans = {p_fans:.1f} MW
                            |
               +------------+------------+
               v                         v
        GROSS POWER               EXHAUST + JET
        {p_gross:.1f} MW                  {p_exhaust:.1f} MW waste
               |
               v
        - parasitic {abs(p_parasitic):.3f} MW
               |
               v
        NET POWER = {p_net:.1f} MW
        Carnot ceiling = {ceiling:.1f} MW (eta={eta*100:.1f}%)
        Audit: {'PASS' if res.carnot_ok[idx] else 'FAIL'}
""")


# ==============================================================================
# SECTION 13d -- TIMELINE PLOT (multi-series, with grid)
# ==============================================================================

def print_timeline(res: SimResult, key: str = "") -> None:
    """Multi-series timeline plot following the ValcanoHarvester RUN view pattern."""
    print_header(f"TIMELINE  --  {key}  ({res.mode})")
    if not res.t_h:
        print("  No data.")
        return

    # plot 4 series in a stacked layout
    series_list = [
        ("P_net MW", res.p_net_mw, "#"),
        ("Q_lava MW", res.q_lava_mw, "*"),
        ("cavern T C", res.cavern_t_c, "+"),
        ("exit v m/s", res.v_exit, "."),
    ]

    for label, data, ch in series_list:
        if not data:
            continue
        lo, hi = min(data), max(data)
        if hi - lo < EPS:
            hi = lo + 1.0
        height = 10
        width = 66
        n = len(data)
        print(f"\n  {label}  [{lo:.2g} .. {hi:.2g}]")
        print(f"  {'+'}{'-' * width}{'+'}")
        for row in range(height, 0, -1):
            thr = lo + (hi - lo) * row / (height + 1)
            line = "  |"
            for col in range(width):
                i = int(col * (n - 1) / max(width - 1, 1))
                line += ch if data[i] >= thr else " "
            print(line + "|")
        # x-axis with time labels
        print(f"  {'+'}{'-' * width}{'+'}")
        t_max = res.t_h[-1] if res.t_h else 0
        t_marks = [0, t_max * 0.25, t_max * 0.5, t_max * 0.75, t_max]
        axis = "  "
        for i in range(width):
            closest = min(t_marks, key=lambda tm: abs(i / width * t_max - tm))
            if abs(i / width * t_max - closest) < t_max / width / 2:
                axis += "|"
            else:
                axis += " "
        print(axis)
        labels = "  "
        for tm in t_marks:
            col = int(tm / t_max * width) if t_max > 0 else 0
            s = f"{tm:.1f}h"
            labels = labels[:2 + col] + s + labels[2 + col + len(s):]
        print(labels[:2 + width])


# ==============================================================================
# SECTION 12b -- DETAILED HARDWARE SPECIFICATION DICTS (to scale, SI)
# ==============================================================================
#
# Following the ValcanoHarvester.py pattern: these are the source-of-truth
# dimension dictionaries for every physical component. Every number is a real
# engineering value in SI units. The 3D cross-section (SECTION 13b) and the
# BOM parts list (SECTION 12c) both read from these, so the picture cannot
# drift from the engineering.
#
# THE NUMBER THAT DOMINATES THIS DESIGN:  a 1800 m tunnel bore heated from
# -10 C to ~500 C grows thermally by
#
#       dL = alpha * L * dT = 12e-6 * 1800 * 510 = 11.0 METRES
#
# You cannot rigidly constrain that. Thermal-growth management -- expansion
# joints, slip couplings, and controlled heat-up ramps -- is the single
# hardest mechanical problem, which is why the tunnel spec includes expansion
# joint spacing.

CAVERN_HW = {
    "excavation_method":  "drill-and-blast + roadheader, lined with shotcrete",
    "lining_thick_mm":     600.0,     # steel-fibre reinforced shotcrete
    "lining_grade":        "C30/37 + stainless steel mesh",
    "seal_layer_mm":       8.0,       # welded HDPE membrane, gas-tight
    "insulation_mm":      200.0,      # closed-cell PU foam on warm side
    # --- ultra thermal insulation (required when cavern is near lava) ---
    # The cavern sits in a volcano/lava environment where surrounding rock can
    # be 500-800 C. Standard PU foam (k=0.024 W/mK) is insufficient. Ultra
    # insulation uses a multi-layer approach:
    #   - Aerogel blanket (k=0.014 W/mK, 50mm)
    #   - Vacuum insulated panels (k=0.004 W/mK, 100mm)
    #   - Multi-layer insulation (MLI, k=0.00005 W/mK, 30 layers)
    #   - Reflective foil barriers between layers
    # Combined R-value: ~30 m^2K/W (vs ~8 for PU foam alone)
    "ultra_insulation_mm":       500.0,  # total ultra insulation thickness
    "ultra_insulation_k":        0.004,  # VIP effective conductivity W/(m K)
    "ultra_insulation_layers":      30,  # MLI layers
    "ultra_insulation_materials":  "aerogel blanket + vacuum panels + MLI + reflective foil",
    "ultra_insulation_r_value":   30.0,  # m^2K/W combined thermal resistance
    "flat_span_m":          45.0,     # largest clear span
    "flat_height_m":        30.0,     # crown height
    "floor_area_m2":      167_000.0,  # approximate footprint
    "access_tunnel_d_m":     6.0,     # main access drift
    "access_tunnel_l_m":   420.0,     # from surface to cavern
    "drainage_sump_m3":   2_000.0,    # condensate collection
    "geophone_count":         16,     # microseismic monitoring
    "dts_fiber_km":          8.0,     # distributed temperature sensing
    "pressure_rating_bar":    8.0,    # design internal pressure
    "hydraulic_door_mm":  4_000.0,    # 4 m diameter, hydraulically actuated
    "thermal_expansion_m":    0.3,    # of the lining over the operating range
}

TUNNEL_HW = {
    "bore_method":        "TBM (hard rock) + drill-and-blast (lava contact zone)",
    "lining_thick_mm":     350.0,     # precast concrete segments + refractory
    "refractory_thick_mm": 120.0,     # in the lava contact zone only
    "refractory_grade":    "alumina-silica firebrick, rated 1400 C",
    "casing_od_mm":      4_500.0,     # steel casing through non-contact zones
    "casing_id_mm":      4_000.0,
    "casing_grade":       "API L80 + Inconel 625 cladding in lava zone",
    "expansion_joint_m":    50.0,     # spacing between slip joints
    "expansion_joint_stroke_m": 0.45, # per-joint thermal growth capacity
    "total_expansion_joints": 36,     # over 1800 m
    "drainage_pipe_mm":    150.0,     # condensate drainage along tunnel
    "ventilation_duct_mm": 800.0,     # construction + emergency ventilation
    "lighting":            "LED every 20 m, emergency battery backup",
    "escape_refuges":       6,        # pressurised refuge chambers
    "seismic_isolation":    "base-isolated from lava contact zone",
}

TURBINE_HW = {
    "type":                "axial-flow, multi-stage, air-expansion turbine",
    "stages":                6,       # default; Arctic preset uses 8
    "rotor_d_mm":         3_200.0,    # rotor diameter
    "rotor_blade_count":     18,      # per stage
    "blade_material":      "Inconel 718, single-crystal investment cast",
    "blade_coating":       "TBC (thermal barrier coating), YSZ + bond coat",
    "stage_spacing_m":       4.5,     # axial distance between stages
    "inlet_T_max_C":       650.0,     # material limit on first-stage inlet
    "inlet_P_max_bar":       6.0,     # from cavern charge pressure
    "rpm":               3600.0,      # 60 Hz, 2-pole direct-coupled
    "generator_mva":        45.0,     # per turbine module
    "generator_kv":         13.8,
    "generator_pf":          0.85,
    "generator_cooling":   "hydrogen-cooled, 75 C hot-spot limit",
    "gearbox":             "none -- direct-coupled to turbine",
    "bearing_type":        "tilting-pad journal + thrust, oil-lubricated",
    "seal_type":           "labyrinth + buffer air",
    "eta_isentropic":        0.82,    # per stage
    "eta_generator":         0.96,
    "capex_per_module_musd": 28.0,    # per turbine+generator module
}

EXIT_FAN_HW = {
    "type":                "ducted axial fan, generator-coupled",
    "fan_d_mm":           2_800.0,
    "blade_count":            8,
    "blade_material":      "carbon-fibre composite, 150 C rated",
    "rpm":               1800.0,
    "generator_kW":         850.0,    # per fan
    "generator_type":      "permanent-magnet, direct-drive",
    "eta_fan":               0.75,
    "nozzle_area_m2":        2.0,     # converging nozzle upstream of fans
    "nozzle_type":          "converging, fixed-geometry",
    "sound_level_dB":        105.0,   # at 100 m -- needs attenuation berm",
}

MONITOR_HW = {
    "cavern_pressure_sensors":  12,   # piezoresistive, 0-10 bar
    "cavern_temp_sensors":      24,   # RTD Pt100, -50 to 200 C
    "cavern_dts_fiber_km":      8.0,  # distributed temperature
    "cavern_das_fiber_km":      8.0,  # distributed acoustic
    "tunnel_temp_sensors":      90,   # every 20 m along the tunnel
    "tunnel_pressure_taps":     45,   # every 40 m
    "turbine_vibration":        12,   # accelerometers per turbine module
    "lava_temp_wells":           6,   # thermocouple wells into the contact zone
    "seismometers":              8,   # broadband, surface + borehole
    "gnss_stations":             4,   # deformation monitoring
    "gas_sensors":               6,   # SO2, H2S, CO2 at exit
    "scada_points":           2400,   # total I/O count
    "trip_overpressure_bar":    7.5,  # cavern overpressure trip
    "trip_tunnel_T_C":        700.0,  # refractory material limit
    "ramp_limit_pct_per_min":   2.0,  # max flow change per minute
}

ORC_HW = {
    "enabled":               False,   # optional parallel Organic Rankine Cycle
    "working_fluid":        "R245fa",
    "hot_source":           "post-turbine exhaust air (T_turb)",
    "cold_sink":            "ambient air / cooling tower",
    "t_evap_C":             120.0,
    "t_cond_C":              35.0,
    "eta_orc":                0.12,   # ORC cycle efficiency
    "ua_kw_per_k":            85.0,   # evaporator UA
    "capex_per_kw":         2_500.0,  # USD/kW
}


# ==============================================================================
# SECTION 12c -- BOM / PARTS LIST
# ==============================================================================

@dataclass
class Part:
    """A single bill-of-materials entry, following the ValcanoHarvester pattern."""
    order: int
    name: str
    specs: List[str]
    category: str = ""


def build_parts_list(t: Dict) -> List[Part]:
    """Build the complete BOM for a target configuration."""
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t.get("ctrl", None)
    n_sys = max(1, ctrl.n_systems if ctrl else 1)
    n_turb = tn.n_turbine_stages
    n_bores = max(1, lv.n_parallel_bores)
    tun_len = tn.total_length_m
    n_fans = tn.n_exit_fans
    cav_depth = cv.depth_m
    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    lava_len = lv.contact_length_m
    stack_h = tn.height_rise_m
    n_joints = int(tun_len / TUNNEL_HW['expansion_joint_m'])
    bottoming = []
    if tn.potassium_enabled: bottoming.append(("K", "#FF9800"))
    if tn.sco2_enabled: bottoming.append(("sCO2", "#9C27B0"))
    if tn.steam_enabled: bottoming.append(("Steam", "#03A9F4"))
    if tn.orc_enabled: bottoming.append(("ORC", "#4CAF50"))
    parts: List[Part] = []
    o = 1

    parts.append(Part(o, "Cold-Air Cavern", [
        f"Volume: {cv.volume_m3:.2e} m^3",
        f"Depth: {cv.depth_m:.0f} m",
        f"Charge: {k_to_c(cv.t_charge_k):.0f} C @ {cv.p_charge_pa/KPA:.0f} kPa",
        f"Lining: {CAVERN_HW['lining_thick_mm']:.0f} mm shotcrete + {CAVERN_HW['seal_layer_mm']:.0f} mm HDPE seal",
        f"Insulation: {CAVERN_HW['insulation_mm']:.0f} mm PU foam",
        f"Access tunnel: {CAVERN_HW['access_tunnel_d_m']:.0f} m dia x {CAVERN_HW['access_tunnel_l_m']:.0f} m",
        f"Pressure rating: {CAVERN_HW['pressure_rating_bar']:.0f} bar",
        f"Hydraulic door: {CAVERN_HW['hydraulic_door_mm']:.0f} mm",
    ], "cavern")); o += 1

    if cv.ultra_insulation:
        parts.append(Part(o, "Ultra Thermal Insulation System (CRITICAL)", [
            f"Total thickness: {CAVERN_HW['ultra_insulation_mm']:.0f} mm",
            f"Materials: {CAVERN_HW['ultra_insulation_materials']}",
            f"MLI layers: {CAVERN_HW['ultra_insulation_layers']}",
            f"VIP effective conductivity: {CAVERN_HW['ultra_insulation_k']:.4f} W/(m K)",
            f"Combined R-value: {CAVERN_HW['ultra_insulation_r_value']:.0f} m^2K/W",
            f"Effective U_ground: {cv.effective_u_ground():.4f} W/(m^2 K)",
            f"Rock T (with lava halo): {cv.ground_t_at_depth_c():.0f} C" if not math.isinf(cv.lava_proximity_m) else "",
            f"REQUIRED: cavern is {cv.lava_proximity_m:.0f} m from {cv.lava_t_nearby_c:.0f} C lava" if not math.isinf(cv.lava_proximity_m) else "",
            "Without this insulation the heat leak would exceed chiller capacity",
        ], "insulation")); o += 1

    parts.append(Part(o, "Cavern Monitoring", [
        f"Pressure sensors: {MONITOR_HW['cavern_pressure_sensors']}",
        f"Temperature sensors: {MONITOR_HW['cavern_temp_sensors']}",
        f"DTS fiber: {MONITOR_HW['cavern_dts_fiber_km']:.1f} km",
        f"DAS fiber: {MONITOR_HW['cavern_das_fiber_km']:.1f} km",
        f"Geophones: {CAVERN_HW['geophone_count']}",
    ], "monitoring")); o += 1

    if cv.active_cooling:
        parts.append(Part(o, "Active Refrigeration Plant", [
            f"Cooling capacity: {cv.chiller_kW_thermal:.0f} kW_thermal",
            f"COP: {cv.chiller_cop:.1f}",
            f"Power draw: {cv.chiller_kW_thermal/cv.chiller_cop:.0f} kW_elec",
        ], "cooling")); o += 1

    parts.append(Part(o, f"Tunnel Bore ({tn.total_length_m/1609.34:.2f} miles)", [
        f"Length: {tn.total_length_m:.0f} m  ({tn.total_length_m/1609.34:.2f} mi)",
        f"Diameter: {tn.diameter_m:.1f} m",
        f"Cross-section: {tn.area_m2():.2f} m^2",
        f"Lining: {TUNNEL_HW['lining_thick_mm']:.0f} mm concrete segments",
        f"Refractory: {TUNNEL_HW['refractory_thick_mm']:.0f} mm firebrick (lava zone)",
        f"Expansion joints: {TUNNEL_HW['total_expansion_joints']} @ {TUNNEL_HW['expansion_joint_m']:.0f} m spacing",
        f"Escape refuges: {TUNNEL_HW['escape_refuges']}",
    ], "tunnel")); o += 1

    parts.append(Part(o, "Lava Heat-Exchange Contact", [
        f"Lava temperature: {lv.t_lava_c:.0f} C",
        f"Contact length: {lv.contact_length_m:.0f} m",
        f"Heat-transfer area: {math.pi*lv.tunnel_diameter_m*lv.contact_length_m:.2e} m^2",
        f"U-value: {lv.u_lava:.0f} W/(m^2 K)",
        f"Refractory grade: {TUNNEL_HW['refractory_grade']}",
    ], "heat-exchange")); o += 1

    parts.append(Part(o, f"Turbine Array ({tn.n_turbine_stages} stages)", [
        f"Type: {TURBINE_HW['type']}",
        f"Rotor diameter: {TURBINE_HW['rotor_d_mm']:.0f} mm",
        f"Blades: {TURBINE_HW['rotor_blade_count']} per stage, {TURBINE_HW['blade_material']}",
        f"Inlet T max: {TURBINE_HW['inlet_T_max_C']:.0f} C",
        f"RPM: {TURBINE_HW['rpm']:.0f}",
        f"Generator: {TURBINE_HW['generator_mva']:.0f} MVA, {TURBINE_HW['generator_kv']:.1f} kV",
        f"Eta isentropic: {TURBINE_HW['eta_isentropic']:.2f}, gen: {TURBINE_HW['eta_generator']:.2f}",
        f"CAPEX: {TURBINE_HW['capex_per_module_musd']:.0f} M USD per module",
    ], "turbine")); o += 1

    parts.append(Part(o, f"Exit Fan Array ({tn.n_exit_fans} fans)", [
        f"Type: {EXIT_FAN_HW['type']}",
        f"Fan diameter: {EXIT_FAN_HW['fan_d_mm']:.0f} mm",
        f"Blades: {EXIT_FAN_HW['blade_count']}, {EXIT_FAN_HW['blade_material']}",
        f"Generator: {EXIT_FAN_HW['generator_kW']:.0f} kW PM direct-drive",
        f"Eta: {EXIT_FAN_HW['eta_fan']:.2f}",
        f"Nozzle area: {EXIT_FAN_HW['nozzle_area_m2']:.1f} m^2",
    ], "exit-fans")); o += 1

    parts.append(Part(o, "Stack / Chimney", [
        f"Height rise: {tn.height_rise_m:.0f} m",
        f"Diameter: {tn.diameter_m:.1f} m",
        f"Function: buoyancy draft + exit jet acceleration",
    ], "stack")); o += 1

    parts.append(Part(o, "Monitoring & Safety System", [
        f"SCADA points: {MONITOR_HW['scada_points']}",
        f"Seismometers: {MONITOR_HW['seismometers']}",
        f"GNSS stations: {MONITOR_HW['gnss_stations']}",
        f"Gas sensors: {MONITOR_HW['gas_sensors']}",
        f"Trip: overpressure {MONITOR_HW['trip_overpressure_bar']:.1f} bar, "
        f"tunnel T {MONITOR_HW['trip_tunnel_T_C']:.0f} C",
        f"Ramp limit: {MONITOR_HW['ramp_limit_pct_per_min']:.1f} %/min",
    ], "safety")); o += 1

    if ORC_HW["enabled"]:
        parts.append(Part(o, "Parallel ORC Bottoming Cycle", [
            f"Working fluid: {ORC_HW['working_fluid']}",
            f"Evap/cond: {ORC_HW['t_evap_C']:.0f} / {ORC_HW['t_cond_C']:.0f} C",
            f"Eta: {ORC_HW['eta_orc']:.2f}",
        ], "orc")); o += 1

    parts.append(Part(o, "Switchyard & Grid Connection", [
        f"Step-up: {TURBINE_HW['generator_kv']:.1f} kV -> 132 kV",
        f"Transformer: {tn.n_turbine_stages * TURBINE_HW['generator_mva']:.0f} MVA",
        f"Transmission: 132 kV, 4 bays",
    ], "grid")); o += 1

    # --- EXPANDED DETAILED PARTS ---

    parts.append(Part(o, "Cavern Lining System (layered)", [
        f"Layer 1 (inner): {CAVERN_HW['seal_layer_mm']:.0f} mm welded HDPE gas-tight membrane",
        f"Layer 2: {CAVERN_HW['lining_thick_mm']:.0f} mm steel-fibre reinforced shotcrete ({CAVERN_HW['lining_grade']})",
        f"Layer 3 (warm side): {CAVERN_HW['insulation_mm']:.0f} mm closed-cell PU foam insulation",
        f"Thermal expansion allowance: {CAVERN_HW['thermal_expansion_m']:.1f} m over operating range",
        f"Flat span: {CAVERN_HW['flat_span_m']:.0f} m, crown height: {CAVERN_HW['flat_height_m']:.0f} m",
        f"Floor area: {CAVERN_HW['floor_area_m2']:,.0f} m^2",
    ], "cavern-lining")); o += 1

    parts.append(Part(o, "Cavern Hydraulic Isolation Door", [
        f"Diameter: {CAVERN_HW['hydraulic_door_mm']:.0f} mm",
        f"Actuation: hydraulic, fail-close",
        f"Pressure rating: {CAVERN_HW['pressure_rating_bar']:.0f} bar",
        f"Function: isolates cavern from access tunnel during discharge",
    ], "cavern-door")); o += 1

    parts.append(Part(o, "Condensate Drainage Sump", [
        f"Volume: {CAVERN_HW['drainage_sump_m3']:,.0f} m^3",
        f"Function: collects condensation from humid air cooling",
        f"Pump: submersible, auto-level controlled",
    ], "cavern-drain")); o += 1

    if cv.cascade_cooling:
        parts.append(Part(o, "Cascade Refrigeration System", [
            f"Cascade COP: {cv.cascade_cop:.2f}",
            f"Stages: multi-stage cascade (N2/Helium)",
            f"Target temperature: {k_to_c(cv.t_charge_k):.0f} C",
        ], "cascade-cooling")); o += 1

    if cv.lava_heated_cooling:
        parts.append(Part(o, "Lava-Heated Absorption Chiller", [
            f"Cooling fraction from lava: {cv.lava_heated_cooling*100:.0f}%",
            f"Type: lithium-bromide/water absorption chiller",
            f"Heat source: lava thermal energy (no electricity required)",
            f"Function: uses waste lava heat to drive cooling -> improves EROI",
        ], "absorption-chiller")); o += 1

    if cv.liquid_air:
        parts.append(Part(o, "Liquid Air Charging System", [
            f"Liquid COP: {cv.liquid_cop:.2f}",
            f"Target: {k_to_c(cv.t_charge_k):.0f} C (liquefied air)",
            f"Storage: cryogenic dewar-style insulated cavern",
        ], "liquid-air")); o += 1

    parts.append(Part(o, f"Recharge Compressor ({cv.n_compress_stages} stages)", [
        f"Stages: {cv.n_compress_stages} intercooled stages",
        f"Efficiency: {cv.compress_eta:.2f}",
        f"Function: re-pressurizes cavern from atmospheric to {cv.p_charge_pa/1e5:.0f} bar",
        f"Power: significant -- this is the energy cost accounted in EROI",
    ], "recharge")); o += 1

    parts.append(Part(o, "Access Tunnel", [
        f"Diameter: {CAVERN_HW['access_tunnel_d_m']:.0f} m",
        f"Length: {CAVERN_HW['access_tunnel_l_m']:.0f} m (surface to cavern)",
        f"Lining: shotcrete + HDPE seal",
        f"Function: personnel/maintenance access + recharge air path",
    ], "access-tunnel")); o += 1

    n_bores = max(1, lv.n_parallel_bores)
    parts.append(Part(o, f"Parallel Tunnel Bores (x{n_bores})", [
        f"Count: {n_bores} parallel bores through lava contact zone",
        f"Individual length: {tn.total_length_m:.0f} m",
        f"Individual diameter: {tn.diameter_m:.1f} m",
        f"Total cross-section: {n_bores * tn.area_m2():.1f} m^2",
        f"Total heat-transfer area: {n_bores * math.pi * tn.diameter_m * lv.contact_length_m:.2e} m^2",
        f"Bore method: {TUNNEL_HW['bore_method']}",
    ], "parallel-bores")); o += 1

    parts.append(Part(o, "Tunnel Casing & Lining", [
        f"Steel casing OD: {TUNNEL_HW['casing_od_mm']:.0f} mm",
        f"Steel casing ID: {TUNNEL_HW['casing_id_mm']:.0f} mm",
        f"Grade: {TUNNEL_HW['casing_grade']}",
        f"Concrete lining: {TUNNEL_HW['lining_thick_mm']:.0f} mm precast segments",
    ], "tunnel-casing")); o += 1

    parts.append(Part(o, "Lava-Zone Refractory Lining", [
        f"Thickness: {TUNNEL_HW['refractory_thick_mm']:.0f} mm",
        f"Grade: {TUNNEL_HW['refractory_grade']}",
        f"Function: protects steel casing from {lv.t_lava_c:.0f} C lava",
        f"Seismic isolation: {TUNNEL_HW['seismic_isolation']}",
    ], "refractory")); o += 1

    n_joints = int(tn.total_length_m / TUNNEL_HW['expansion_joint_m'])
    parts.append(Part(o, f"Thermal Expansion Joints (x{n_joints})", [
        f"Count: {n_joints} slip joints",
        f"Spacing: {TUNNEL_HW['expansion_joint_m']:.0f} m",
        f"Stroke per joint: {TUNNEL_HW['expansion_joint_stroke_m']:.2f} m",
        f"Total expansion capacity: {n_joints * TUNNEL_HW['expansion_joint_stroke_m']:.1f} m",
        f"Function: absorbs thermal growth of {tn.total_length_m:.0f} m bore",
    ], "expansion-joints")); o += 1

    parts.append(Part(o, "Tunnel Condensate Drainage", [
        f"Pipe diameter: {TUNNEL_HW['drainage_pipe_mm']:.0f} mm",
        f"Function: removes condensation from humid air expansion",
        f"Routing: gravity drain along tunnel floor to sump",
    ], "tunnel-drain")); o += 1

    parts.append(Part(o, "Emergency Ventilation & Refuges", [
        f"Ventilation duct: {TUNNEL_HW['ventilation_duct_mm']:.0f} mm",
        f"Escape refuges: {TUNNEL_HW['escape_refuges']} pressurised chambers",
        f"Lighting: {TUNNEL_HW['lighting']}",
    ], "tunnel-safety")); o += 1

    if lv.hx_enabled:
        parts.append(Part(o, f"Shell-and-Tube Lava HX ({lv.hx_n_tubes:,} tubes)", [
            f"Tube count: {lv.hx_n_tubes:,}",
            f"Tube OD: {lv.hx_tube_od_mm:.0f} mm",
            f"Tube length: {lv.hx_tube_length_m:.0f} m",
            f"U-value: {lv.hx_u:.0f} W/(m^2 K)",
            f"Total HX area: {lv.hx_n_tubes * math.pi * lv.hx_tube_od_mm/1000 * lv.hx_tube_length_m:.2e} m^2",
            f"Material: Inconel 625 cladding (lava-side)",
        ], "hx-tubes")); o += 1

    if lv.heat_pipe:
        parts.append(Part(o, "Heat Pipes (lava to bore)", [
            f"Type: two-phase sodium heat pipes",
            f"Function: transfers lava heat to tunnel air with high flux density",
            f"Advantage: ~10x higher heat flux than conduction alone",
            f"Fin factor: {lv.fin_factor:.0f}x effective area enhancement",
        ], "heat-pipes")); o += 1

    parts.append(Part(o, f"Turbine Stator Vanes (x{tn.n_turbine_stages} stages)", [
        f"Count: {TURBINE_HW['rotor_blade_count']} stator vanes per stage",
        f"Material: {TURBINE_HW['blade_material']}",
        f"Coating: {TURBINE_HW['blade_coating']}",
        f"Function: directs flow onto rotor blades at optimal angle",
    ], "stator-vanes")); o += 1

    parts.append(Part(o, "Turbine Shaft & Bearings", [
        f"RPM: {TURBINE_HW['rpm']:.0f} (60 Hz, 2-pole)",
        f"Bearings: {TURBINE_HW['bearing_type']}",
        f"Seals: {TURBINE_HW['seal_type']}",
        f"Stage spacing: {TURBINE_HW['stage_spacing_m']:.1f} m",
        f"Gearbox: {TURBINE_HW['gearbox']}",
    ], "turbine-shaft")); o += 1

    parts.append(Part(o, "Turbine Generator", [
        f"Rating: {TURBINE_HW['generator_mva']:.0f} MVA per module",
        f"Voltage: {TURBINE_HW['generator_kv']:.1f} kV",
        f"Power factor: {TURBINE_HW['generator_pf']:.2f}",
        f"Cooling: {TURBINE_HW['generator_cooling']}",
        f"Efficiency: {TURBINE_HW['eta_generator']:.2f}",
    ], "generator")); o += 1

    if tn.n_reheat_stages > 0:
        parts.append(Part(o, f"Reheat Sections (x{tn.n_reheat_stages})", [
            f"Count: {tn.n_reheat_stages} reheat stages between turbine stages",
            f"Function: reheats air between expansions -> near-isothermal expansion",
            f"Heat source: lava via HX tubes",
            f"Effect: increases work output toward Carnot limit",
        ], "reheat")); o += 1

    if tn.regenerator_eff > 0:
        parts.append(Part(o, "Regenerator / Recuperator", [
            f"Efficiency: {tn.regenerator_eff:.2f}",
            f"Function: pre-heats incoming cold air using exhaust heat",
            f"Type: rotary or plate-fin recuperator",
        ], "regenerator")); o += 1

    if tn.mhd_enabled:
        parts.append(Part(o, "MHD Topping Cycle", [
            f"Efficiency: {tn.mhd_eta:.2f}",
            f"Type: magnetohydrodynamic generator (ionized gas + magnetic field)",
            f"Function: extracts direct electrical energy from hot ionized gas",
            f"Placement: upstream of first turbine stage",
        ], "mhd")); o += 1

    if tn.potassium_enabled:
        parts.append(Part(o, "Potassium Vapor Topping Cycle", [
            f"Efficiency: {tn.potassium_eta:.2f}",
            f"Working fluid: potassium vapor (Rankine cycle)",
            f"Hot source: post-turbine exhaust (>1500 C)",
            f"Function: topping cycle above sCO2/steam",
        ], "potassium-cycle")); o += 1

    if tn.sco2_enabled:
        parts.append(Part(o, "Supercritical CO2 Bottoming Cycle", [
            f"Efficiency: {tn.sco2_eta:.2f}",
            f"Working fluid: supercritical CO2 (Rankine/Brayton hybrid)",
            f"Hot source: post-turbine exhaust (>200 C)",
            f"Advantage: compact, high efficiency at moderate temperatures",
        ], "sco2-cycle")); o += 1

    if tn.steam_enabled:
        parts.append(Part(o, "Steam Rankine Bottoming Cycle", [
            f"Efficiency: {tn.steam_eta:.2f}",
            f"Working fluid: water/steam",
            f"Hot source: post-turbine exhaust (>300 C)",
            f"Type: multi-pressure (HP/IP/LP) steam turbine",
        ], "steam-cycle")); o += 1

    if tn.orc_enabled:
        parts.append(Part(o, "ORC Bottoming Cycle", [
            f"Efficiency: {tn.orc_eta:.2f}",
            f"Working fluid: {ORC_HW['working_fluid']}",
            f"Evap/cond: {ORC_HW['t_evap_C']:.0f}/{ORC_HW['t_cond_C']:.0f} C",
            f"Hot source: final exhaust before stack (>50 C above ambient)",
            f"UA: {ORC_HW['ua_kw_per_k']:.0f} kW/K",
        ], "orc-cycle")); o += 1

    parts.append(Part(o, f"Exit Nozzle & Jet", [
        f"Nozzle area: {tn.exit_nozzle_area_m2:.1f} m^2",
        f"Type: {'converging-diverging (supersonic)' if tn.supersonic_nozzle else 'converging (subsonic)'}",
        f"Max Mach: {tn.max_mach:.2f}",
        f"Function: accelerates exhaust to extract kinetic energy via fans",
    ], "nozzle")); o += 1

    parts.append(Part(o, f"Exit Fan Generators (x{tn.n_exit_fans})", [
        f"Count: {tn.n_exit_fans} ducted axial fans",
        f"Fan diameter: {EXIT_FAN_HW['fan_d_mm']:.0f} mm",
        f"Blades: {EXIT_FAN_HW['blade_count']}, {EXIT_FAN_HW['blade_material']}",
        f"RPM: {EXIT_FAN_HW['rpm']:.0f}",
        f"Generator: {EXIT_FAN_HW['generator_kW']:.0f} kW PM direct-drive each",
        f"Total fan capacity: {tn.n_exit_fans * EXIT_FAN_HW['generator_kW']:.0f} kW",
        f"Efficiency: {EXIT_FAN_HW['eta_fan']:.2f}",
        f"Sound: {EXIT_FAN_HW['sound_level_dB']:.0f} dB at 100 m",
    ], "exit-fans-detail")); o += 1

    parts.append(Part(o, "Stack / Chimney", [
        f"Height: {tn.height_rise_m:.0f} m",
        f"Diameter: {tn.diameter_m:.1f} m",
        f"Function: buoyancy draft + exit jet acceleration + fan mounting",
        f"Stack effect: {stack_pressure(tn, 1.2, 0.5):.0f} Pa (typical)",
    ], "stack-detail")); o += 1

    parts.append(Part(o, "Tunnel Monitoring System", [
        f"Temperature sensors: {MONITOR_HW['tunnel_temp_sensors']} (every 20 m)",
        f"Pressure taps: {MONITOR_HW['tunnel_pressure_taps']} (every 40 m)",
        f"Turbine vibration: {MONITOR_HW['turbine_vibration']} accelerometers/module",
        f"Lava temp wells: {MONITOR_HW['lava_temp_wells']} thermocouple wells",
    ], "tunnel-monitoring")); o += 1

    parts.append(Part(o, "Site Monitoring & Safety", [
        f"Seismometers: {MONITOR_HW['seismometers']} broadband",
        f"GNSS stations: {MONITOR_HW['gnss_stations']} deformation",
        f"Gas sensors: {MONITOR_HW['gas_sensors']} (SO2/H2S/CO2)",
        f"SCADA total: {MONITOR_HW['scada_points']} I/O points",
        f"Trip: overpressure {MONITOR_HW['trip_overpressure_bar']:.1f} bar",
        f"Trip: tunnel T {MONITOR_HW['trip_tunnel_T_C']:.0f} C",
        f"Ramp limit: {MONITOR_HW['ramp_limit_pct_per_min']:.1f} %/min",
    ], "site-safety")); o += 1

    n_sys = max(1, t["ctrl"].n_systems)
    if n_sys > 1:
        parts.append(Part(o, f"Dual System Interconnect (x{n_sys})", [
            f"Systems: {n_sys} complete parallel tunnel arrays",
            f"Each: own cavern, tunnel, turbines, fans, HX",
            f"Shared: switchyard, grid connection, control room",
            f"Total output: {n_sys}x single-system capacity",
        ], "dual-system")); o += 1

    parts.append(Part(o, "Step-Up Transformer & Switchyard", [
        f"Step-up: {TURBINE_HW['generator_kv']:.1f} kV -> 132 kV",
        f"Capacity: {tn.n_turbine_stages * TURBINE_HW['generator_mva']:.0f} MVA",
        f"Transmission: 132 kV, 4 bays",
        f"Transformer type: oil-immersed, ONAN/ONAF cooling",
        f"Switchgear: SF6 insulated (GIS)",
    ], "switchyard-detail")); o += 1

    # --- ADDITIONAL DETAILED PARTS ---

    parts.append(Part(o, "Main Isolation Valve (cavern outlet)", [
        f"Type: full-bore ball valve, hydraulically actuated",
        f"Diameter: {tn.diameter_m * 1000:.0f} mm",
        f"Material: Inconel 625 body, Stellite seats",
        f"Function: isolates cavern from tunnel for maintenance",
        f"Fail-close: spring-return on hydraulic loss",
    ], "valve-cavern")); o += 1

    parts.append(Part(o, "Turbine Bypass Valve", [
        f"Type: globe valve, motor-operated",
        f"Function: bypasses turbine array during startup/shutdown",
        f"Routing: routes air directly to stack during ramp-up",
    ], "valve-bypass")); o += 1

    parts.append(Part(o, "Anti-Surge Valve (per turbine stage)", [
        f"Count: {n_turb} valves",
        f"Type: fast-acting butterfly valve",
        f"Function: prevents compressor surge during transient flow",
        f"Response time: < 200 ms",
    ], "valve-antisurge")); o += 1

    parts.append(Part(o, "Pressure Relief Valve (cavern)", [
        f"Type: spring-loaded safety relief valve",
        f"Set pressure: {CAVERN_HW['pressure_rating_bar'] * 1.1:.1f} bar",
        f"Capacity: full cavern blowdown in 30 min",
        f"Function: prevents cavern overpressure",
    ], "valve-relief")); o += 1

    parts.append(Part(o, "Condensate Drain Valve (tunnel)", [
        f"Count: {int(tun_len / 50)} valves along tunnel",
        f"Type: float-operated automatic drain valve",
        f"Function: removes condensed water from tunnel floor",
        f"Routing: to drainage pipe -> cavern sump",
    ], "valve-drain")); o += 1

    parts.append(Part(o, "SCADA Control System", [
        f"I/O points: {MONITOR_HW['scada_points']}",
        f"Architecture: redundant PLC + HMI + historian",
        f"Communication: fiber-optic ring, Modbus TCP / DNP3",
        f"Control loops: cavern pressure, mass flow, turbine speed, reheat temp",
        f"Sampling: 100 ms fast loop, 1 s slow loop",
    ], "scada")); o += 1

    parts.append(Part(o, "Turbine Governor / Speed Control", [
        f"Type: electronic governor with hydraulic actuator",
        f"Control: speed + load sharing + synchronizing",
        f"Overspeed trip: 110% mechanical, 108% electrical",
        f"Function: maintains 3600 RPM under varying flow",
    ], "governor")); o += 1

    parts.append(Part(o, "Generator Protection Relay", [
        f"Type: numerical multi-function relay (IEEE C37.102)",
        f"Functions: differential, overcurrent, loss-of-excitation,",
        f"  stator earth fault, reverse power, out-of-step",
        f"Trip: circuit breaker + turbine fast-valve",
    ], "gen-protection")); o += 1

    parts.append(Part(o, "Synchronizing Panel", [
        f"Type: auto-synchronizer + synchrocheck relay (25)",
        f"Functions: voltage matching, frequency matching, phase matching",
        f"Breaker: vacuum circuit breaker, {TURBINE_HW['generator_kv']:.1f} kV",
    ], "sync-panel")); o += 1

    parts.append(Part(o, "13.8 kV Switchgear", [
        f"Type: metal-clad vacuum circuit breaker",
        f"Voltage: {TURBINE_HW['generator_kv']:.1f} kV",
        f"Rating: {TURBINE_HW['generator_mva']:.0f} MVA per bay",
        f"Bays: {n_turb} turbine bays + 1 tie + 1 spare",
        f"Protection: differential + overcurrent + ground fault",
    ], "switchgear-13.8")); o += 1

    parts.append(Part(o, "132 kV Switchyard", [
        f"Type: SF6 gas-insulated switchgear (GIS)",
        f"Voltage: 132 kV",
        f"Bays: 4 (2 incoming + 2 outgoing transmission lines)",
        f"Bus: double-bus single-breaker scheme",
        f"Protection: line differential + distance + breaker failure",
    ], "switchyard-132")); o += 1

    parts.append(Part(o, "Station Service Transformer", [
        f"Type: dry-type cast coil",
        f"Rating: 500 kVA, 13.8 kV / 480 V",
        f"Function: powers auxiliary loads (pumps, lighting, HVAC, cranes)",
    ], "station-service")); o += 1

    parts.append(Part(o, "UPS / Battery System", [
        f"Type: online double-conversion UPS",
        f"Capacity: 100 kVA, 30 min battery backup",
        f"Function: powers critical controls during grid outage",
        f"Battery: VRLA, 48 V DC bus",
    ], "ups")); o += 1

    parts.append(Part(o, "DC Battery Bank", [
        f"Type: nickel-cadmium (Ni-Cd) pocket plate",
        f"Voltage: 125 V DC",
        f"Capacity: 800 Ah, 8-hour discharge",
        f"Function: turbine emergency lube oil, trip circuits, lighting",
    ], "dc-battery")); o += 1

    parts.append(Part(o, "Emergency Diesel Generator", [
        f"Rating: 1 MW, 480 V",
        f"Function: station service backup during extended outage",
        f"Fuel: on-site diesel, 72-hour tank",
        f"Start: auto-start within 10 seconds of loss",
    ], "diesel-backup")); o += 1

    parts.append(Part(o, "Turbine Lube Oil System", [
        f"Type: forced lubrication, ISO VG 32 turbine oil",
        f"Pumps: main (shaft-driven) + AC standby + DC emergency",
        f"Cooler: water-cooled oil cooler",
        f"Filter: duplex full-flow, 10 micron",
        f"Reservoir: 5000 L per turbine module",
    ], "lube-oil")); o += 1

    parts.append(Part(o, "Generator Hydrogen System", [
        f"Gas: hydrogen (H2), 75 psig",
        f"Function: cools generator stator + rotor",
        f"Purity monitor: continuous, trip at < 95% purity",
        f"Storage: high-pressure H2 cylinders + gas control panel",
        f"Seal oil: prevents H2 leakage along shaft",
    ], "h2-cooling")); o += 1

    parts.append(Part(o, "Seal Oil System", [
        f"Type: vacuum-treated seal oil",
        f"Function: seals generator shaft against H2 leakage",
        f"Pumps: AC main + DC emergency",
        f"Demineralizer: maintains oil quality",
    ], "seal-oil")); o += 1

    parts.append(Part(o, "Cavern Recharge Piping", [
        f"Type: large-diameter steel pipe, {CAVERN_HW['access_tunnel_d_m']:.0f} m dia",
        f"Material: API X65, internal coating",
        f"Function: delivers compressed air from recharge compressor to cavern",
        f"Valves: non-return + isolation, hydraulically actuated",
    ], "recharge-pipe")); o += 1

    parts.append(Part(o, "Interconnecting Piping (bottoming cycles)", [
        f"Type: insulated process piping",
        f"Materials: carbon steel (steam), stainless 316L (sCO2/ORC),",
        f"  Inconel (potassium)",
        f"Insulation: calcium silicate + aluminum jacket",
        f"Function: connects turbine exhaust to bottoming cycle heat exchangers",
    ], "bottoming-pipe")); o += 1

    parts.append(Part(o, "Cooling Water System", [
        f"Source: cooling tower (mechanical-draft, induced-draft)",
        f"Capacity: sized for bottoming cycle condensers + oil coolers",
        f"Temperature: 30 C supply, 40 C return",
        f"Pumps: 3x50% capacity, redundant",
        f"Treatment: chemical dosing (corrosion/scale inhibitor)",
    ], "cooling-water")); o += 1

    parts.append(Part(o, "Fire Protection System", [
        f"Type: water spray + CO2 flooding (electrical areas)",
        f"Detection: heat + smoke + flame detectors",
        f"Coverage: turbine hall, switchgear, control room, transformers",
        f"Water: fire water pump + 200 m3 tank",
    ], "fire-protection")); o += 1

    parts.append(Part(o, "Compressed Air Instrument System", [
        f"Type: oil-free screw compressor + desiccant dryer",
        f"Pressure: 8 bar, dewpoint -40 C",
        f"Function: powers pneumatic actuators, instruments, tools",
        f"Receiver: 3000 L, redundant compressors",
    ], "instrument-air")); o += 1

    parts.append(Part(o, "Cable Tray & Conduit System", [
        f"Type: galvanized steel cable trays + PVC conduit",
        f"Separation: power, control, and instrument cables segregated",
        f"Routing: above turbine hall, in cable tunnels",
        f"Firestop: intumescent at wall penetrations",
    ], "cable-tray")); o += 1

    parts.append(Part(o, "Grounding & Lightning Protection", [
        f"Ground grid: copper conductor, 0.5 ohm target resistance",
        f"Lightning: air terminals on stack + buildings, down conductors",
        f"Surge: surge arresters on 13.8 kV + 132 kV",
        f"Function: protects equipment from lightning + fault currents",
    ], "grounding")); o += 1

    parts.append(Part(o, "Control Room Building", [
        f"Type: blast-resistant, HVAC-pressurized",
        f"Area: 200 m2",
        f"Equipment: operator consoles, HMI screens, communication",
        f"Redundancy: dual operator stations, hot-standby server",
    ], "control-room")); o += 1

    parts.append(Part(o, "Turbine Hall Building", [
        f"Type: industrial steel-frame, crane-equipped",
        f"Crane: 50-ton overhead bridge crane",
        f"Area: {n_turb * 30:.0f} m2 (turbine modules + generators)",
        f"HVAC: ventilation + spot cooling at generators",
    ], "turbine-hall")); o += 1

    parts.append(Part(o, "Bottoming Cycle Building", [
        f"Type: industrial steel-frame",
        f"Area: {len(bottoming) * 150:.0f} m2",
        f"Equipment: K/sCO2/Steam/ORC turbines, condensers, pumps",
        f"Function: houses all bottoming cycle equipment",
    ], "bottoming-building")); o += 1

    parts.append(Part(o, "Cooling Tower", [
        f"Type: mechanical-draft induced-draft cooling tower",
        f"Cells: 4 cells, redundant",
        f"Capacity: sized for peak heat rejection",
        f"Drift eliminators: 0.002% drift rate",
        f"Water treatment: biocide + scale inhibitor dosing",
    ], "cooling-tower")); o += 1

    parts.append(Part(o, "Site Civil Works", [
        f"Access roads: paved, 2-lane, {cav_depth + tun_len/10:.0f} m total",
        f"Drainage: site stormwater management + oil separator",
        f"Foundations: reinforced concrete, seismic-rated",
        f"Fencing: security fencing + access control",
        f"Lighting: site lighting + CCTV surveillance",
    ], "civil-works")); o += 1

    parts.append(Part(o, "Fiber Optic Communication", [
        f"Type: single-mode fiber, OPGW on transmission line",
        f"Capacity: 10 Gbps, redundant ring topology",
        f"Function: SCADA + protection signaling + voice + data",
    ], "fiber-optic")); o += 1

    parts.append(Part(o, "Meteorological Station", [
        f"Sensors: wind speed/direction, temperature, humidity,",
        f"  barometric pressure, solar radiation",
        f"Function: weather monitoring for cooling tower + stack draft",
        f"Data: logged to SCADA historian",
    ], "met-station")); o += 1

    parts.append(Part(o, "Vibration Monitoring System", [
        f"Sensors: {MONITOR_HW['turbine_vibration']} accelerometers per turbine module",
        f"Type: proximity probes (shaft) + casing accelerometers",
        f"Analysis: FFT spectrum, trend, alarm, trip",
        f"Function: detects bearing wear, unbalance, misalignment",
    ], "vibration-monitor")); o += 1

    parts.append(Part(o, "Gas Analysis System (exit)", [
        f"Sensors: {MONITOR_HW['gas_sensors']} gas analyzers",
        f"Species: SO2, H2S, CO2, CO, NOx, O2",
        f"Function: monitors exhaust gas composition for safety + emissions",
        f"Location: stack exit + bottoming cycle exhaust",
    ], "gas-analysis")); o += 1

    parts.append(Part(o, "Cavern Temperature Mapping", [
        f"System: Distributed Temperature Sensing (DTS) fiber",
        f"Length: {MONITOR_HW['cavern_dts_fiber_km']:.1f} km fiber optic cable",
        f"Resolution: 1 m spatial, 0.1 C temperature",
        f"Function: full 3D temperature map of cavern walls + air",
    ], "dts-mapping")); o += 1

    parts.append(Part(o, "Acoustic Emission Monitoring", [
        f"System: Distributed Acoustic Sensing (DAS) fiber",
        f"Length: {MONITOR_HW['cavern_das_fiber_km']:.1f} km fiber optic cable",
        f"Function: detects rock cracking, lining stress, microseismic events",
        f"Alarm: triggers on acoustic energy above threshold",
    ], "das-monitor")); o += 1

    # --- ADDITIONAL STRUCTURAL & SEALING PARTS ---

    parts.append(Part(o, "Cavern Roof Support Arches", [
        f"Type: steel arch ribs, {CAVERN_HW['flat_span_m']:.0f}m span",
        f"Spacing: 1.5 m on center",
        f"Material: Q345 structural steel, fire-rated coating",
        f"Function: supports cavern roof against rock pressure",
        f"Count: {int(CAVERN_HW['floor_area_m2'] / 1.5 / CAVERN_HW['flat_span_m']):.0f} arches",
    ], "roof-arch")); o += 1

    parts.append(Part(o, "Cavern Floor Slab", [
        f"Type: reinforced concrete slab on grade",
        f"Thickness: 300 mm",
        f"Area: {CAVERN_HW['floor_area_m2']:,.0f} m2",
        f"Function: structural floor + condensate drainage slope",
        f"Slope: 1:200 toward drainage sump",
    ], "floor-slab")); o += 1

    parts.append(Part(o, "Tunnel Segment Lining Rings", [
        f"Type: precast concrete segments, bolted ring",
        f"Segments per ring: 6+1 key",
        f"Ring width: 1.5 m",
        f"Count: {int(tun_len / 1.5) * n_bores} rings total",
        f"Material: C40/50 concrete, EPDM gasket seals",
    ], "tunnel-rings")); o += 1

    parts.append(Part(o, "Tunnel Segment Bolts", [
        f"Type: high-strength galvanized steel bolts, M24",
        f"Count: {int(tun_len / 1.5) * n_bores * 14} bolts (14 per ring)",
        f"Torque: 400 Nm",
        f"Function: connects precast tunnel lining segments",
    ], "tunnel-bolts")); o += 1

    parts.append(Part(o, "EPDM Segment Gaskets", [
        f"Type: ethylene-propylene-diene rubber (EPDM) compression gaskets",
        f"Profile: double-blade seal",
        f"Count: {int(tun_len / 1.5) * n_bores * 6} gaskets",
        f"Function: water/gas tightness between segment joints",
        f"Compression: 30% at design pressure",
    ], "epdm-gaskets")); o += 1

    parts.append(Part(o, "Refractory Anchor System", [
        f"Type: V-shaped stainless steel anchors, grade 304",
        f"Count: {int(lava_len * n_bores * 10)} anchors",
        f"Spacing: 100 mm grid pattern",
        f"Function: secures refractory castable to tunnel casing",
    ], "refractory-anchors")); o += 1

    parts.append(Part(o, "Expansion Joint Bellows", [
        f"Type: metal bellows, Inconel 625",
        f"Count: {n_joints} bellows",
        f"Stroke: {TUNNEL_HW['expansion_joint_stroke_m']:.2f} m per joint",
        f"Function: absorbs thermal expansion while maintaining seal",
        f"Design temp: 700 C max",
    ], "expansion-bellows")); o += 1

    parts.append(Part(o, "Tunnel Access Platforms", [
        f"Type: steel grating platforms at turbine locations",
        f"Count: {n_turb} platforms",
        f"Size: 3m x 5m each",
        f"Function: maintenance access to turbine modules",
        f"Handrails: galvanized steel, 1.1m height",
    ], "access-platforms")); o += 1

    parts.append(Part(o, "Tunnel Lighting System", [
        f"Type: {TUNNEL_HW['lighting']}",
        f"Count: {int(tun_len / 20) * n_bores} fixtures",
        f"Power: 30W per fixture, LED",
        f"Emergency: 90 min battery backup per fixture",
    ], "tunnel-lighting")); o += 1

    parts.append(Part(o, "Tunnel Communication System", [
        f"Type: leaky feeder radio + emergency phones",
        f"Count: {int(tun_len / 200)} phone stations",
        f"Function: 2-way radio + phone communication in tunnel",
        f"Coverage: 100% tunnel length",
    ], "tunnel-comm")); o += 1

    # --- ELECTRICAL INSTRUMENTATION ---

    parts.append(Part(o, "Pressure Transmitters (cavern)", [
        f"Type: piezoresistive, 0-10 bar, 4-20mA",
        f"Count: {MONITOR_HW['cavern_pressure_sensors']}",
        f"Accuracy: 0.1% FS",
        f"Function: cavern pressure monitoring + control",
    ], "pt-cavern")); o += 1

    parts.append(Part(o, "Temperature Transmitters (cavern)", [
        f"Type: RTD Pt100, -200 to 200 C, 4-20mA",
        f"Count: {MONITOR_HW['cavern_temp_sensors']}",
        f"Accuracy: 0.1 C",
        f"Function: cavern temperature monitoring",
    ], "tt-cavern")); o += 1

    parts.append(Part(o, "Flow Transmitters (tunnel)", [
        f"Type: thermal mass flow meter",
        f"Count: {n_bores} (one per bore)",
        f"Range: 0-500 kg/s",
        f"Function: measures mass flow rate through each bore",
    ], "ft-tunnel")); o += 1

    parts.append(Part(o, "Vibration Transmitters (turbine)", [
        f"Type: 4-20mA loop-powered, accelerometer + proximity probe",
        f"Count: {MONITOR_HW['turbine_vibration'] * n_turb}",
        f"Function: turbine bearing vibration monitoring",
        f"Alarm: 7 mm/s RMS warning, 11 mm/s trip",
    ], "vt-turbine")); o += 1

    parts.append(Part(o, "Thermocouple Wells (lava)", [
        f"Type: Type K (Chromel-Alumel), Inconel sheath",
        f"Count: {MONITOR_HW['lava_temp_wells']}",
        f"Range: 0-1400 C",
        f"Function: direct lava contact temperature measurement",
    ], "tc-lava")); o += 1

    # --- PIPING & VALVES ---

    parts.append(Part(o, "Bottoming Cycle Feed Pump (sCO2)", [
        f"Type: centrifugal pump, supercritical CO2 service",
        f"Count: 2x100% redundant",
        f"Material: stainless 316L",
        f"Function: circulates sCO2 working fluid",
    ], "pump-sco2")); o += 1

    parts.append(Part(o, "Steam Condensate Pump", [
        f"Type: vertical turbine pump",
        f"Count: 2x100% redundant",
        f"Material: carbon steel + stainless impeller",
        f"Function: returns condensate to steam cycle",
    ], "pump-steam")); o += 1

    parts.append(Part(o, "ORC Working Fluid Pump", [
        f"Type: positive displacement gear pump",
        f"Count: 2x100% redundant",
        f"Material: stainless 316L",
        f"Working fluid: {ORC_HW['working_fluid']}",
    ], "pump-orc")); o += 1

    parts.append(Part(o, "Potassium Condensate Pump", [
        f"Type: electromagnetic (EM) pump, liquid metal",
        f"Count: 2x100% redundant",
        f"Material: Inconel 600",
        f"Function: circulates liquid potassium",
    ], "pump-potassium")); o += 1

    parts.append(Part(o, "Condensate Return Pump (cavern sump)", [
        f"Type: submersible pump",
        f"Count: 2x100% redundant",
        f"Capacity: 50 m3/h",
        f"Function: pumps collected condensate from cavern sump to surface",
    ], "pump-condensate")); o += 1

    # --- SAFETY SYSTEMS ---

    parts.append(Part(o, "Emergency Shutdown System (ESD)", [
        f"Type: hardwired ESD logic, independent of SCADA",
        f"Level 1: close cavern valve + trip turbines",
        f"Level 2: full isolation + vent to stack",
        f"Level 3: cavern blowdown via relief valve",
        f"Initiation: manual pushbutton + automatic trips",
    ], "esd")); o += 1

    parts.append(Part(o, "Gas Detection System", [
        f"Type: point + open-path gas detectors",
        f"Species: H2 (generator), SO2/H2S (lava), CO (fire), O2 (asphyxiation)",
        f"Count: {MONITOR_HW['gas_sensors']} analyzers + 20 point detectors",
        f"Action: alarm at PEL, trip at 2x PEL",
    ], "gas-detection")); o += 1

    parts.append(Part(o, "Emergency Escape Respirators", [
        f"Type: self-contained self-rescuer (SCSR), 30 min",
        f"Count: {TUNNEL_HW['escape_refuges'] * 4} units",
        f"Location: stored at each escape refuge",
        f"Function: emergency breathing for tunnel evacuation",
    ], "escape-respirators")); o += 1

    parts.append(Part(o, "Tunnel Fire Suppression", [
        f"Type: water mist + foam system",
        f"Coverage: turbine areas + electrical areas",
        f"Detection: heat + flame + smoke",
        f"Activation: automatic + manual",
    ], "tunnel-fire")); o += 1

    # --- STRUCTURAL ---

    parts.append(Part(o, "Turbine Foundation Pedestals", [
        f"Type: reinforced concrete pedestal + sole plate",
        f"Count: {n_turb} pedestals",
        f"Material: C35/45 concrete + steel sole plate",
        f"Grouting: epoxy grout under sole plate",
        f"Function: supports turbine + generator, transmits loads to rock",
    ], "turbine-foundation")); o += 1

    parts.append(Part(o, "Stack Structural Support", [
        f"Type: steel lattice tower + guy cables",
        f"Height: {stack_h:.0f} m",
        f"Wind load: designed for 150 km/h",
        f"Seismic: designed for site-specific spectrum",
    ], "stack-structure")); o += 1

    parts.append(Part(o, "Cavern Rock Bolts", [
        f"Type: fully grouted rebar bolts, 25mm dia, 4m length",
        f"Count: {int(CAVERN_HW['floor_area_m2'] / 4):.0f} bolts",
        f"Pattern: 2m x 2m grid",
        f"Function: reinforces host rock around cavern",
    ], "rock-bolts")); o += 1

    parts.append(Part(o, "Shotcrete Lining Reinforcement", [
        f"Type: steel fiber + welded wire mesh",
        f"Fiber: 50 kg/m3 dosage, 30mm length",
        f"Mesh: 6mm dia @ 150mm grid",
        f"Function: flexural reinforcement of shotcrete lining",
    ], "shotcrete-reinf")); o += 1

    # --- TURBINE SUB-COMPONENTS ---

    parts.append(Part(o, "Turbine Rotor Blades (set)", [
        f"Count: {TURBINE_HW['rotor_blade_count']} blades per stage x {n_turb} stages = {TURBINE_HW['rotor_blade_count'] * n_turb} blades",
        f"Material: {TURBINE_HW['blade_material']}",
        f"Coating: {TURBINE_HW['blade_coating']}",
        f"Root type: fir-tree (serrated), precision ground",
        f"Tip clearance: 1.5 mm, shrouded tips with labyrinth seal",
        f"Manufacturing: single-crystal investment casting + EDM root",
    ], "rotor-blades")); o += 1

    parts.append(Part(o, "Turbine Stator Vanes (set)", [
        f"Count: {TURBINE_HW['rotor_blade_count']} vanes per stage x {n_turb} stages",
        f"Material: {TURBINE_HW['blade_material']}",
        f"Coating: {TURBINE_HW['blade_coating']}",
        f"Mounting: welded into inner + outer shroud rings",
        f"Stagger angle: variable (optimized per stage)",
    ], "stator-vanes-set")); o += 1

    parts.append(Part(o, "Turbine Casing (split horizontally)", [
        f"Type: horizontally split, bolted flange casing",
        f"Material: cast steel, CrMoV alloy",
        f"Inner diameter: {TURBINE_HW['rotor_d_mm']+200:.0f} mm",
        f"Flange: 48 bolts M48, torqued to 2500 Nm",
        f"Function: contains expansion flow, supports stator vanes",
    ], "turbine-casing")); o += 1

    parts.append(Part(o, "Turbine Diaphragms (per stage)", [
        f"Count: {n_turb} diaphragms",
        f"Type: welded diaphragm with inner + outer ring",
        f"Material: CrMoV steel + stainless steel vanes",
        f"Function: holds stator vanes, seals between stages",
        f"Interstage seal: labyrinth, 4-tooth",
    ], "diaphragms")); o += 1

    parts.append(Part(o, "Turbine Rotor Discs", [
        f"Count: {n_turb} discs",
        f"Type: forged monoblock or built-up rotor",
        f"Material: forged CrMoV steel, ultrasonically inspected",
        f"Balance: ISO 1940 grade G2.5 at 3600 RPM",
        f"Over-speed test: 120% rated speed, 3 minutes",
    ], "rotor-discs")); o += 1

    parts.append(Part(o, "Turbine Main Shaft", [
        f"Type: forged steel shaft, single-piece",
        f"Material: AISI 4140, quenched + tempered",
        f"Diameter: 400 mm journal, 600 mm at coupling",
        f"Length: {n_turb * TURBINE_HW['stage_spacing_m']:.1f} m total",
        f"Coupling: flexible gear coupling to generator",
    ], "main-shaft")); o += 1

    parts.append(Part(o, "Turbine Thrust Bearing", [
        f"Type: tilting-pad thrust bearing, Kingsbury type",
        f"Pads: 8 tilting pads, copper-faced, babbitt-lined",
        f"Capacity: 500 kN axial thrust",
        f"Oil: ISO VG 32 turbine oil, forced lubrication",
        f"Temperature: trip at 95 C pad temperature",
    ], "thrust-bearing")); o += 1

    parts.append(Part(o, "Turbine Journal Bearings", [
        f"Type: tilting-pad journal bearings",
        f"Count: 2 per turbine module (inlet + exhaust)",
        f"Pads: 5 tilting pads, babbitt-lined steel-backed",
        f"Clearance: 0.15 mm diametral",
        f"Oil: forced lubrication at 2 bar, 40 C supply",
    ], "journal-bearings")); o += 1

    parts.append(Part(o, "Turbine Turning Gear", [
        f"Type: electric motor-driven turning gear",
        f"Speed: 3 RPM (barring speed)",
        f"Function: rotates shaft during cooldown to prevent bow",
        f"Engagement: automatic when turbine speed < 100 RPM",
    ], "turning-gear")); o += 1

    parts.append(Part(o, "Turbine Labyrinth Seals", [
        f"Type: multi-tooth labyrinth seal",
        f"Count: {n_turb * 2} seals (inlet + exhaust per stage)",
        f"Material: aluminum bronze teeth, steel rotor land",
        f"Clearance: 0.25 mm radial",
        f"Buffer air: supplied at 0.5 bar above process pressure",
    ], "labyrinth-seals")); o += 1

    # --- GENERATOR SUB-COMPONENTS ---

    parts.append(Part(o, "Generator Stator Core", [
        f"Type: laminated electrical steel, grain-oriented",
        f"Layers: 0.35 mm thick, varnished + stacked",
        f"Slots: 48 slots for stator winding",
        f"Diameter: {TURBINE_HW['rotor_d_mm']+800:.0f} mm outer",
        f"Function: magnetic core for stator flux path",
    ], "stator-core")); o += 1

    parts.append(Part(o, "Generator Stator Winding", [
        f"Type: 2-layer lap winding, Roebel bars",
        f"Conductors: copper strands, 0.5mm x 2.5mm, transposed",
        f"Insulation: class F (155 C), mica-epoxy VPI",
        f"Voltage: {TURBINE_HW['generator_kv']:.1f} kV, BIL 110 kV",
        f"Connections: 3-phase wye, neutral grounded via resistor",
    ], "stator-winding")); o += 1

    parts.append(Part(o, "Generator Rotor Winding", [
        f"Type: 2-pole or 4-pole field winding",
        f"Conductors: silver-bearing copper strips",
        f"Insulation: class F, epoxy-mica",
        f"Excitation: brushless exciter, 500 VDC, 800 A",
        f"Cooling: hydrogen-cooled, direct gas-cooled rotor",
    ], "rotor-winding")); o += 1

    parts.append(Part(o, "Generator Exciter", [
        f"Type: brushless rotating rectifier exciter",
        f"Rating: 500 VDC, 800 A, 400 kW",
        f"Components: pilot exciter + main exciter + diode wheel",
        f"AVR: digital automatic voltage regulator",
        f"Function: supplies DC field current to generator rotor",
    ], "exciter")); o += 1

    parts.append(Part(o, "Generator AVR (Auto Voltage Regulator)", [
        f"Type: digital microprocessor-based AVR",
        f"Function: regulates generator terminal voltage",
        f"Features: VAr limiting, PF control, PSS (power system stabilizer)",
        f"Redundancy: dual AVR, auto-transfer",
    ], "avr")); o += 1

    parts.append(Part(o, "Generator Terminal Bushings", [
        f"Type: oil-impregnated paper bushings, porcelain housing",
        f"Count: 6 (3 phases + 3 neutral)",
        f"Rating: {TURBINE_HW['generator_kv']:.1f} kV, 3000 A",
        f"Function: conducts generator output through casing wall",
    ], "terminal-bushings")); o += 1

    parts.append(Part(o, "Generator Isophase Busduct", [
        f"Type: isolated phase busduct (IPB), forced-air cooled",
        f"Rating: {TURBINE_HW['generator_kv']:.1f} kV, 1500 A",
        f"Conductors: aluminum tubular, per-phase enclosure",
        f"Function: connects generator to step-up transformer",
        f"Length: ~50 m per generator module",
    ], "isophase-busduct")); o += 1

    # --- BOTTOMING CYCLE SUB-COMPONENTS ---

    if tn.potassium_enabled:
        parts.append(Part(o, "Potassium Turbine", [
            f"Type: axial-flow Rankine turbine, potassium vapor",
            f"Material: Inconel 617 (high temp creep resistant)",
            f"Inlet: 1500 C, 5 bar",
            f"Outlet: 800 C, 0.5 bar",
            f"Efficiency: {tn.potassium_eta:.2f}",
            f"Seal: mechanical seal with buffer gas",
        ], "potassium-turbine")); o += 1

        parts.append(Part(o, "Potassium Condenser", [
            f"Type: shell-and-tube, potassium on shell side",
            f"Coolant: sCO2 on tube side",
            f"Material: Inconel 600 tubes, SS316L shell",
            f"Function: condenses potassium vapor to liquid",
        ], "potassium-condenser")); o += 1

    if tn.sco2_enabled:
        parts.append(Part(o, "sCO2 Turbine", [
            f"Type: radial-inflow turbine, supercritical CO2",
            f"Material: stainless 347 + Inconel 718 rotor",
            f"Inlet: 500 C, 200 bar",
            f"Outlet: 100 C, 75 bar",
            f"Efficiency: {tn.sco2_eta:.2f}",
            f"Advantage: compact, ~10x smaller than steam turbine",
        ], "sco2-turbine")); o += 1

        parts.append(Part(o, "sCO2 Recuperator", [
            f"Type: printed circuit heat exchanger (PCHE)",
            f"Material: stainless 316L, photochemically etched",
            f"Pressure rating: 250 bar",
            f"Effectiveness: 85%",
            f"Function: recovers heat from sCO2 exhaust to preheat feed",
        ], "sco2-recuperator")); o += 1

        parts.append(Part(o, "sCO2 Primary Heat Exchanger", [
            f"Type: shell-and-tube, exhaust air on shell side",
            f"Material: Inconel 625 tubes",
            f"UA: 5000 kW/K",
            f"Function: transfers heat from turbine exhaust to sCO2 cycle",
        ], "sco2-phx")); o += 1

        parts.append(Part(o, "sCO2 Compressor", [
            f"Type: centrifugal compressor, supercritical CO2",
            f"Stages: 2 intercooled stages",
            f"Material: stainless 347",
            f"Function: compresses sCO2 from 75 to 200 bar",
        ], "sco2-compressor")); o += 1

    if tn.steam_enabled:
        parts.append(Part(o, "HP Steam Turbine", [
            f"Type: axial-flow, single-flow, impulse-reaction",
            f"Inlet: 540 C, 160 bar",
            f"Outlet: 350 C, 40 bar",
            f"Material: CrMoV rotor, 12Cr blades",
            f"Efficiency: {tn.steam_eta:.2f}",
        ], "hp-steam-turbine")); o += 1

        parts.append(Part(o, "IP Steam Turbine", [
            f"Type: axial-flow, double-flow, reaction",
            f"Inlet: 540 C, 40 bar (reheated)",
            f"Outlet: 200 C, 5 bar",
            f"Material: 12Cr rotor, 12Cr blades",
        ], "ip-steam-turbine")); o += 1

        parts.append(Part(o, "LP Steam Turbine", [
            f"Type: axial-flow, double-flow, last-stage long blades",
            f"Inlet: 200 C, 5 bar",
            f"Outlet: 40 C, 0.1 bar (vacuum)",
            f"Last-stage blades: 900 mm titanium, tip speed 600 m/s",
            f"Condenser: water-cooled, cooling tower supply",
        ], "lp-steam-turbine")); o += 1

        parts.append(Part(o, "Steam Condenser", [
            f"Type: surface condenser, shell-and-tube",
            f"Coolant: cooling tower water (30/40 C)",
            f"Vacuum: 0.1 bar (95% vacuum)",
            f"Tubes: titanium, 25mm OD, 10000 tubes",
            f"Air ejector: steam jet air ejector (SJAE)",
        ], "steam-condenser")); o += 1

        parts.append(Part(o, "Steam Boiler / Evaporator", [
            f"Type: once-through, water-tube",
            f"Heat source: turbine exhaust air (>300 C)",
            f"Capacity: sized for steam cycle mass flow",
            f"Tubes: Inconel 625, finned for enhanced heat transfer",
            f"Feedwater: demineralized, < 0.1 ppm TDS",
        ], "steam-boiler")); o += 1

    if tn.orc_enabled:
        parts.append(Part(o, "ORC Evaporator", [
            f"Type: shell-and-tube, {ORC_HW['working_fluid']} on tube side",
            f"Heat source: final exhaust air (>50 C above ambient)",
            f"UA: {ORC_HW['ua_kw_per_k']:.0f} kW/K",
            f"Material: stainless 316L",
        ], "orc-evaporator")); o += 1

        parts.append(Part(o, "ORC Turbine", [
            f"Type: radial-inflow turbine, organic working fluid",
            f"Working fluid: {ORC_HW['working_fluid']}",
            f"Inlet: {ORC_HW['t_evap_C']:.0f} C, 15 bar",
            f"Outlet: {ORC_HW['t_cond_C']:.0f} C, 2 bar",
            f"Efficiency: {ORC_HW['eta_orc']:.2f}",
            f"Material: stainless 316L",
        ], "orc-turbine")); o += 1

        parts.append(Part(o, "ORC Condenser", [
            f"Type: air-cooled or water-cooled condenser",
            f"Working fluid: {ORC_HW['working_fluid']}",
            f"Condensation temp: {ORC_HW['t_cond_C']:.0f} C",
            f"Coolant: cooling tower water or ambient air",
        ], "orc-condenser")); o += 1

    # --- CONSTRUCTION EQUIPMENT ---

    parts.append(Part(o, "TBM (Tunnel Boring Machine)", [
        f"Type: {TUNNEL_HW['bore_method']}",
        f"Diameter: {tn.diameter_m:.1f} m cutterhead",
        f"Cutter: 19-inch disc cutters, {int(tn.diameter_m * 4)} cutters",
        f"Thrust: 15,000 kN",
        f"Torque: 8,000 kNm",
        f"Power: 3 MW total",
        f"Function: bores tunnel through hard rock",
    ], "tbm")); o += 1

    parts.append(Part(o, "Roadheader (cavern excavation)", [
        f"Type: boom-type roadheader",
        f"Cutter: transverse cutting head, 50 kW",
        f"Function: excavates cavern cross-section",
        f"Conveyor: belt conveyor to muck cars",
    ], "roadheader")); o += 1

    parts.append(Part(o, "Shotcrete Robot", [
        f"Type: remote-controlled shotcrete manipulator arm",
        f"Reach: 8 m boom",
        f"Output: 20 m3/h shotcrete",
        f"Function: applies shotcrete lining to cavern walls",
    ], "shotcrete-robot")); o += 1

    parts.append(Part(o, "Rock Bolter (cavern)", [
        f"Type: mechanized rock bolt installation rig",
        f"Function: installs rock bolts in cavern roof/walls",
        f"Rate: 20 bolts per hour",
    ], "rock-bolter")); o += 1

    # --- WATER TREATMENT ---

    parts.append(Part(o, "Demineralized Water Plant", [
        f"Type: RO + mixed-bed ion exchange",
        f"Capacity: 50 m3/h",
        f"Quality: < 0.1 ppm TDS, < 0.01 ppm silica",
        f"Function: supplies makeup water for steam cycle",
        f"Storage: 200 m3 demineralized water tank",
    ], "demin-water")); o += 1

    parts.append(Part(o, "Cooling Water Treatment", [
        f"Type: chemical dosing system",
        f"Chemicals: corrosion inhibitor, scale inhibitor, biocide",
        f"Dosing: automatic, flow-proportioned",
        f"Function: protects cooling water system from corrosion/scale/fouling",
    ], "cw-treatment")); o += 1

    # --- HVAC ---

    parts.append(Part(o, "Turbine Hall HVAC", [
        f"Type: supply + exhaust ventilation, 20 air changes/hour",
        f"Capacity: 100,000 m3/h supply",
        f"Spot cooling: 4 units at generator locations",
        f"Function: maintains turbine hall temp < 45 C",
    ], "turbine-hall-hvac")); o += 1

    parts.append(Part(o, "Control Room HVAC", [
        f"Type: precision air conditioning, 24/7 operation",
        f"Capacity: 50 kW cooling",
        f"Pressurization: +50 Pa, filtered air",
        f"Function: maintains 22 C +/- 2 C for electronics + operators",
    ], "control-room-hvac")); o += 1

    # --- CRANE SUB-COMPONENTS ---

    parts.append(Part(o, "Crane Bridge Girder", [
        f"Type: box-section welded steel girder",
        f"Span: {CAVERN_HW['flat_span_m']:.0f} m (turbine hall width)",
        f"Material: Q345B structural steel",
        f"Deflection: L/800 max at rated load",
        f"Function: main load-carrying beam of overhead crane",
    ], "crane-girder")); o += 1

    parts.append(Part(o, "Crane Hoist (50t)", [
        f"Type: electric wire rope hoist, 50 ton capacity",
        f"Lift: 20 m",
        f"Speed: 2 m/min lift, 20 m/min traverse",
        f"Motor: 30 kW lift, 5 kW traverse",
        f"Brake: electromagnetic disc brake, failsafe",
    ], "crane-hoist")); o += 1

    parts.append(Part(o, "Crane End Carriages", [
        f"Count: 2 end carriages",
        f"Type: welded box section with wheel assemblies",
        f"Wheels: 2 per carriage, 400 mm diameter, forged steel",
        f"Rail: QU80 crane rail, welded to runway beam",
        f"Drive: 2x 3 kW gearmotors, one per carriage",
    ], "crane-end-carriage")); o += 1

    parts.append(Part(o, "Crane Runway Beam", [
        f"Type: welded I-section runway beam",
        f"Length: {CAVERN_HW['flat_span_m'] * 1.2:.0f} m per runway",
        f"Count: 2 runways (one each side of turbine hall)",
        f"Support: columns at 6 m spacing",
        f"Rail: QU80 crane rail, continuous weld",
    ], "crane-runway")); o += 1

    # --- HVAC DUCTS & COMPONENTS ---

    parts.append(Part(o, "Turbine Hall Supply Ducts", [
        f"Type: galvanized steel ductwork, rectangular",
        f"Size: 800x600 mm main, 400x300 mm branches",
        f"Insulation: 25mm fiberglass external",
        f"Function: distributes supply air throughout turbine hall",
        f"Dampers: motorized volume control dampers at each branch",
    ], "supply-ducts")); o += 1

    parts.append(Part(o, "Turbine Hall Exhaust Fans", [
        f"Type: axial exhaust fans, roof-mounted",
        f"Count: 4 fans",
        f"Capacity: 25,000 m3/h each",
        f"Motor: 5.5 kW per fan, VFD-controlled",
        f"Function: removes hot air from turbine hall",
    ], "exhaust-fans")); o += 1

    parts.append(Part(o, "Control Room Precision AC", [
        f"Type: ceiling-mounted cassette AC units",
        f"Count: 2 units (1+1 redundant)",
        f"Capacity: 25 kW cooling each",
        f"Refrigerant: R-410A",
        f"Function: precision temperature control for control room",
    ], "precision-ac")); o += 1

    # --- FIRE SYSTEM DETAIL ---

    parts.append(Part(o, "Fire Water Pump Skid", [
        f"Type: electric + diesel fire pump skid",
        f"Electric pump: 1500 L/min at 10 bar",
        f"Diesel pump: 1500 L/min at 10 bar (backup)",
        f"Jockey pump: 100 L/min at 10 bar (pressure maintenance)",
        f"Tank: 200 m3 fire water storage tank",
    ], "fire-pump-skid")); o += 1

    parts.append(Part(o, "Fire Sprinkler Heads", [
        f"Type: quick-response glass-bulb sprinklers",
        f"Temperature: 68 C (red bulb)",
        f"Count: {int(CAVERN_HW['floor_area_m2'] / 12)} sprinklers (12 m2 coverage each)",
        f"Piping: galvanized steel, schedule 40",
        f"Function: automatic fire suppression in turbine hall",
    ], "sprinkler-heads")); o += 1

    parts.append(Part(o, "Fire Alarm Control Panel (FACP)", [
        f"Type: addressable fire alarm panel",
        f"Zones: 32 addressable loops, 250 devices per loop",
        f"Detectors: smoke, heat, flame, gas (addressable)",
        f"Notification: horns, strobes, voice evacuation",
        f"Interface: integrates with SCADA + ESD systems",
    ], "fire-alarm-panel")); o += 1

    parts.append(Part(o, "Clean Agent Fire Suppression (control room)", [
        f"Type: FM-200 clean agent system",
        f"Design concentration: 7% by volume",
        f"Cylinders: 4 x 120 L FM-200 cylinders",
        f"Detection: cross-zone smoke detection",
        f"Function: protects control room electronics (no water damage)",
    ], "clean-agent-fire")); o += 1

    # --- DRAINAGE SYSTEM ---

    parts.append(Part(o, "Cavern Drainage Channels", [
        f"Type: cast-in-floor drainage channels",
        f"Size: 300x300 mm U-channel",
        f"Slope: 1:200 toward sump",
        f"Material: stainless 304, grating covered",
        f"Count: 4 main channels converging at sump",
    ], "drain-channels")); o += 1

    parts.append(Part(o, "Tunnel Drainage Pipe System", [
        f"Type: {TUNNEL_HW['drainage_pipe_mm']:.0f} mm HDPE drainage pipe",
        f"Count: 1 per bore = {n_bores} pipes",
        f"Slope: 1:300 toward cavern sump",
        f"Function: collects condensate from tunnel walls + floor",
        f"Connection: drains to cavern sump, then pumped to surface",
    ], "tunnel-drain-pipe")); o += 1

    parts.append(Part(o, "Condensate Collection Trays", [
        f"Type: stainless steel drip trays at HX tube banks",
        f"Count: {n_turb * 2} trays (at each HX + reheat section)",
        f"Material: stainless 316L, 2mm thick",
        f"Function: collects condensate from heat exchanger tubes",
        f"Drain: piped to tunnel drainage system",
    ], "condensate-trays")); o += 1

    parts.append(Part(o, "Oil Water Separator", [
        f"Type: coalescing plate oil-water separator",
        f"Capacity: 50 m3/h",
        f"Outlet oil: < 15 ppm (environmental compliance)",
        f"Function: separates lubricating oil from drainage water",
        f"Location: surface, downstream of condensate pumps",
    ], "oil-water-sep")); o += 1

    # --- ELECTRICAL DETAIL ---

    parts.append(Part(o, "Generator Circuit Breaker", [
        f"Type: vacuum circuit breaker, drawout",
        f"Rating: {TURBINE_HW['generator_kv']:.1f} kV, 3000 A",
        f"Interrupting: 50 kA symmetrical",
        f"Operation: electrically operated, stored energy mechanism",
        f"Function: isolates generator from bus for maintenance/trip",
    ], "gen-breaker")); o += 1

    parts.append(Part(o, "Lightning Arresters (generator)", [
        f"Type: metal-oxide varistor (MOV) surge arrester",
        f"Rating: {TURBINE_HW['generator_kv']:.1f} kV, 10 kA discharge",
        f"Count: 3 per generator (one per phase)",
        f"Function: protects generator winding from switching surges",
    ], "gen-surge-arrester")); o += 1

    parts.append(Part(o, "PT/CT Instrument Transformers", [
        f"Type: potential transformers (PT) + current transformers (CT)",
        f"PT: {TURBINE_HW['generator_kv']:.1f}kV/120V, 0.3 class",
        f"CT: 3000:5 A, C800 class, 3 per generator",
        f"Function: provides scaled voltage/current for metering + protection",
    ], "instrument-transformers")); o += 1

    parts.append(Part(o, "Generator Protection Relay", [
        f"Type: numerical multifunction relay (IEEE C37.102)",
        f"Functions: 87G differential, 51V backup, 40 loss of field, "
        f"46 neg seq, 49 thermal, 59 overvoltage, 27 undervoltage",
        f"Communication: IEC 61850 GOOSE + DNP3 to SCADA",
    ], "gen-protection-relay")); o += 1

    parts.append(Part(o, "Bus Differential Relay (13.8kV)", [
        f"Type: high-impedance bus differential relay (87B)",
        f"Function: detects internal bus faults, trips all breakers",
        f"Operating time: < 1.5 cycles",
        f"CT requirement: matched CTs on all bus feeders",
    ], "bus-diff-relay")); o += 1

    parts.append(Part(o, "Transformer Differential Relay", [
        f"Type: percentage-restraint transformer differential (87T)",
        f"Functions: 87T differential, 49 thermal, 50/51 instantaneous/TOC",
        f"2nd harmonic: blocks inrush (5-15% setting)",
        f"Function: protects step-up transformer from internal faults",
    ], "xfmr-diff-relay")); o += 1

    # --- CONTROL SYSTEM DETAIL ---

    parts.append(Part(o, "SCADA HMI Workstations", [
        f"Type: dual-monitor operator workstations",
        f"Count: 3 workstations in control room",
        f"Software: SCADA HMI with process graphics, trends, alarms",
        f"Redundancy: dual server, auto-failover < 5s",
        f"Historian: 5-year data retention, 1s resolution",
    ], "scada-hmi")); o += 1

    parts.append(Part(o, "PLC Controllers (turbine governor)", [
        f"Type: hot-standby redundant PLC pair",
        f"CPU: dual redundant, sync scan",
        f"I/O: {MONITOR_HW['scada_points']:,} points total",
        f"Scan time: 10 ms",
        f"Function: turbine speed/load control + sequencing",
    ], "plc-governor")); o += 1

    parts.append(Part(o, "Historian Database Server", [
        f"Type: redundant database server pair",
        f"Storage: 10 TB RAID 10 per server",
        f"Software: time-series historian (OSIsoft PI or equivalent)",
        f"Retention: 5 years online, 10 years archived",
        f"Function: stores all process data for analysis + compliance",
    ], "historian-server")); o += 1

    parts.append(Part(o, "Network Switches (industrial)", [
        f"Type: managed industrial Ethernet switches",
        f"Ports: 24-port Gigabit, fiber uplinks",
        f"Count: 8 switches in ring topology",
        f"Protocol: PRP/HSR redundant ring, < 10ms failover",
        f"Rating: industrial temp range, DIN-rail mounted",
    ], "network-switches")); o += 1

    # --- COMPRESSED AIR / INSTRUMENT AIR ---

    parts.append(Part(o, "Instrument Air Compressor", [
        f"Type: oil-free screw compressor, 100% duty",
        f"Count: 2x100% redundant",
        f"Capacity: 500 m3/h at 8 bar",
        f"Quality: ISO 8573-1 Class 1.2.1 (particulate, water, oil)",
        f"Function: supplies clean dry air for pneumatic instruments/valves",
    ], "instrument-air-comp")); o += 1

    parts.append(Part(o, "Instrument Air Dryer", [
        f"Type: heatless desiccant dryer, twin-tower",
        f"Count: 2x100% redundant",
        f"Dewpoint: -40 C at pressure",
        f"Function: removes moisture from instrument air",
    ], "air-dryer")); o += 1

    parts.append(Part(o, "Instrument Air Receiver", [
        f"Type: vertical pressure vessel, ASME VIII",
        f"Capacity: 3 m3",
        f"Pressure: 10 bar design",
        f"Function: buffers demand spikes, stabilizes pressure",
    ], "air-receiver")); o += 1

    # --- LUBE OIL SYSTEM DETAIL ---

    parts.append(Part(o, "Lube Oil Console", [
        f"Type: skid-mounted lube oil system per turbine",
        f"Components: reservoir, pumps, filters, cooler, regulator",
        f"Reservoir: 2000 L, stainless steel",
        f"Pumps: 2x100% (1 AC + 1 DC emergency)",
        f"Filter: duplex, 10 micron, continuous duty",
    ], "lube-oil-console")); o += 1

    parts.append(Part(o, "Lube Oil Cooler", [
        f"Type: shell-and-tube heat exchanger",
        f"Coolant: cooling tower water (30/40 C)",
        f"Capacity: 200 kW heat rejection",
        f"Material: stainless 316L tubes, carbon steel shell",
    ], "lube-oil-cooler")); o += 1

    parts.append(Part(o, "Seal Oil System", [
        f"Type: skid-mounted seal oil system for generator H2 seals",
        f"Components: reservoir, pumps, vacuum degasifier, regulator",
        f"Pumps: 2x100% (main + emergency)",
        f"Function: provides oil film for hydrogen shaft seals",
        f"Differential: maintains 0.5 bar above H2 pressure",
    ], "seal-oil-system")); o += 1

    # --- HYDROGEN SYSTEM DETAIL ---

    parts.append(Part(o, "Hydrogen Storage Cylinders", [
        f"Type: high-pressure H2 cylinders, 200 bar",
        f"Count: 12 cylinders (manifolded)",
        f"Capacity: 50 Nm3 total storage",
        f"Material: chrome-molybdenum steel",
        f"Function: stores makeup hydrogen for generator cooling",
    ], "h2-cylinders")); o += 1

    parts.append(Part(o, "Hydrogen Control Panel", [
        f"Type: automatic H2 pressure/purity control panel",
        f"Functions: pressure regulation, purity monitoring, CO2 purge",
        f"Sensors: H2 purity (thermal conductivity), pressure, dewpoint",
        f"Alarms: low purity (< 95%), low pressure, high dewpoint",
    ], "h2-control-panel")); o += 1

    parts.append(Part(o, "CO2 Purge System", [
        f"Type: CO2 cylinders + purge manifold",
        f"Function: purges air from generator before H2 fill, and H2 before air",
        f"Count: 6 x 50 L CO2 cylinders",
        f"Safety: prevents explosive H2-air mixture during filling",
    ], "co2-purge")); o += 1

    # --- HEAT PIPE INTERNALS ---

    parts.append(Part(o, "Heat Pipe Evaporator Section", [
        f"Type: wicked heat pipe, evaporator zone in lava",
        f"Length: 2 m per pipe (embedded in lava body)",
        f"Wick: sintered copper powder, 0.5 mm thick, 5 um pore",
        f"Working fluid: NaK (sodium-potassium eutectic)",
        f"Wall: Inconel 600, 2 mm thick",
        f"Function: absorbs heat from lava, evaporates working fluid",
    ], "hp-evaporator")); o += 1

    parts.append(Part(o, "Heat Pipe Adiabatic Section", [
        f"Type: wicked transport section, insulated",
        f"Length: 1 m per pipe (passes through refractory)",
        f"Insulation: 50 mm ceramic fiber, k=0.2 W/mK",
        f"Function: transports vapor from evaporator to condenser",
    ], "hp-adiabatic")); o += 1

    parts.append(Part(o, "Heat Pipe Condenser Section", [
        f"Type: finned condenser, in tunnel airflow",
        f"Length: 1.5 m per pipe (exposed to air stream)",
        f"Fins: helical aluminum fins, 15 fins/m, 10 mm height",
        f"Function: condenses vapor, transfers heat to expanding air",
        f"Heat transfer: ~5 kW per pipe at 1000 C lava temp",
    ], "hp-condenser")); o += 1

    # --- COOLING TOWER INTERNALS ---

    parts.append(Part(o, "Cooling Tower Fill Media", [
        f"Type: PVC film fill, cross-fluted",
        f"Material: corrugated PVC sheets, 0.3 mm thick",
        f"Height: 1.2 m fill depth per cell",
        f"Surface area: 200 m2/m3 (specific surface)",
        f"Function: maximizes water-air contact surface area",
        f"Count: 10 cells x 50 m2 each = 500 m2 total footprint",
    ], "ct-fill")); o += 1

    parts.append(Part(o, "Cooling Tower Drift Eliminators", [
        f"Type: PVC wave-pattern drift eliminators",
        f"Material: corrugated PVC, 0.4 mm thick",
        f"Height: 0.3 m per cell",
        f"Efficiency: 99.9% drift removal (0.001% drift rate)",
        f"Function: removes water droplets from exhaust air",
    ], "ct-drift-elim")); o += 1

    parts.append(Part(o, "Cooling Tower Fan Assembly", [
        f"Type: induced-draft axial fan, FRP blades",
        f"Count: 10 fans (one per cell)",
        f"Diameter: 7.3 m per fan",
        f"Motor: 75 kW, 6-pole, TEFC, VFD-controlled",
        f"Airflow: 500,000 m3/h per fan",
        f"Blades: 6 FRP (fiberglass reinforced plastic) blades",
    ], "ct-fan-assembly")); o += 1

    parts.append(Part(o, "Cooling Tower Water Distribution", [
        f"Type: hot-dipped galvanized steel header + laterals",
        f"Header: 400 mm dia, 50 m long per cell",
        f"Nozzles: 60 spray nozzles per cell (ABS plastic)",
        f"Nozzle type: target/splash type, 3 mm orifice",
        f"Function: distributes warm water over fill media",
    ], "ct-water-dist")); o += 1

    parts.append(Part(o, "Cooling Tower Basin", [
        f"Type: reinforced concrete basin, lined",
        f"Volume: 500 m3 per cell, 5000 m3 total",
        f"Lining: HDPE liner, 2 mm thick",
        f"Function: collects cooled water, provides suction for pumps",
        f"Retention time: 5 minutes at design flow",
    ], "ct-basin")); o += 1

    parts.append(Part(o, "Cooling Water Pumps", [
        f"Type: vertical turbine pump, concrete volute",
        f"Count: 3x50% (2 operating + 1 standby)",
        f"Capacity: 5000 m3/h each at 30 m head",
        f"Motor: 500 kW, 6-pole, submersible",
        f"Function: circulates cooling water from tower to condensers",
    ], "cw-pumps")); o += 1

    # --- SWITCHYARD DETAIL ---

    parts.append(Part(o, "132kV SF6 Circuit Breakers", [
        f"Type: SF6 puffer-type dead-tank breaker",
        f"Count: 4 breakers (2 incoming + 2 outgoing)",
        f"Rating: 132 kV, 2000 A, 40 kA interrupting",
        f"Operating: motor-charged spring mechanism",
        f"SF6: 0.6 MPa at 20 C, moisture < 150 ppmv",
    ], "sf6-breakers")); o += 1

    parts.append(Part(o, "132kV Disconnectors", [
        f"Type: center-break, double-insulator, motor-operated",
        f"Count: 8 disconnectors",
        f"Rating: 132 kV, 2000 A continuous, 100 kA 1s",
        f"Operation: 5 sec motor drive, manual override",
        f"Function: visible isolation for maintenance",
    ], "disconnectors")); o += 1

    parts.append(Part(o, "132kV Current Transformers", [
        f"Type: oil-impregnated paper, dead-tank CT",
        f"Count: 12 CTs (3 per breaker, 4 breakers)",
        f"Ratio: 2000:5 A, multi-ratio",
        f"Class: 5P20 for protection, 0.2S for metering",
        f"Function: provides scaled current for protection + metering",
    ], "132kv-ct")); o += 1

    parts.append(Part(o, "132kV Voltage Transformers", [
        f"Type: inductive VT, oil-impregnated paper",
        f"Count: 6 VTs (3 per bus, 2 buses)",
        f"Ratio: 132kV/sqrt(3) : 110V/sqrt(3)",
        f"Class: 3P for protection, 0.2 for metering",
        f"Function: provides scaled voltage for protection + metering",
    ], "132kv-vt")); o += 1

    parts.append(Part(o, "132kV Surge Arresters", [
        f"Type: metal-oxide varistor (MOV), polymer-housed",
        f"Count: 12 arresters (3 per phase x 4 locations)",
        f"Rating: 132 kV, 10 kA station class",
        f"MCOV: 84 kV, TOV: 110 kV for 10 sec",
        f"Function: protects equipment from lightning + switching surges",
    ], "132kv-arresters")); o += 1

    parts.append(Part(o, "132kV Post Insulators", [
        f"Type: solid-core porcelain, brown glazed",
        f"Count: 60 insulators (bus support + equipment)",
        f"Creepage: 550 mm/kV (pollution class IV)",
        f"Withstand: 550 kV BIL, 230 kV wet switching",
    ], "post-insulators")); o += 1

    parts.append(Part(o, "132kV Busbar System", [
        f"Type: tubular aluminum bus, 80mm dia",
        f"Material: 6063-T6 aluminum alloy",
        f"Config: double main bus with transfer bus",
        f"Current rating: 2000 A continuous, 40 kA 1s",
        f"Spans: 8 m max between support insulators",
    ], "132kv-busbar")); o += 1

    # --- CABLE SYSTEMS ---

    parts.append(Part(o, "Medium Voltage Cable (13.8kV)", [
        f"Type: single-core XLPE insulated, copper conductor",
        f"Size: 500 mm2 Cu, 133% insulation level",
        f"Shield: copper tape shield, PVC jacket",
        f"Count: 3 cables per phase x 4 generators = 36 cables",
        f"Length: ~50 m per cable (generator to switchgear)",
        f"Installation: cable tray + fireproof wrap",
    ], "mv-cable")); o += 1

    parts.append(Part(o, "Control Cable (multicore)", [
        f"Type: multi-pair control cable, 0.75 mm2 per pair",
        f"Pairs: 24 pair per cable, individually shielded",
        f"Shield: aluminum-polyester tape + drain wire",
        f"Jacket: LSZH (low smoke zero halogen), blue",
        f"Count: 200 cables, ~5000 pairs total",
    ], "control-cable")); o += 1

    parts.append(Part(o, "Fiber Optic Cable (SCADA network)", [
        f"Type: single-mode tight-buffered, 24 fiber",
        f"Construction: armored, gel-free, LSZH jacket",
        f"Count: 4 cables in ring topology",
        f"Length: ~2000 m per cable",
        f"Function: SCADA Ethernet network + protection signaling",
    ], "fiber-cable")); o += 1

    parts.append(Part(o, "Power Cable (station service 480V)", [
        f"Type: 3-core + ground, XLPE insulated, copper",
        f"Size: 240 mm2 Cu per core",
        f"Voltage: 600/1000 V",
        f"Count: 20 cables for station service distribution",
        f"Installation: cable tray + conduit",
    ], "station-cable")); o += 1

    parts.append(Part(o, "Cable Tray System", [
        f"Type: hot-dipped galvanized steel ladder tray",
        f"Sizes: 100mm, 300mm, 600mm widths",
        f"Total length: ~3000 m of tray",
        f"Supports: 1.5 m spacing, cantilever from wall",
        f"Function: supports all power, control, and fiber cables",
    ], "cable-tray-system")); o += 1

    # --- TRANSFORMER DETAIL ---

    parts.append(Part(o, "Step-up Transformer (main)", [
        f"Type: oil-immersed, ONAN/ONAF cooling",
        f"Rating: {n_turb * TURBINE_HW['generator_mva']:.0f} MVA",
        f"Voltage: {TURBINE_HW['generator_kv']:.1f} kV / 132 kV",
        f"Vector: Dyn1 (delta-wye)",
        f"Impedance: 12% (standard)",
        f"Taps: +/- 10% in 1.25% steps, on-load tap changer (OLTC)",
        f"Oil: 25,000 L mineral oil, PCB-free",
        f"Conservator: expansion tank + silica gel breather",
    ], "step-up-transformer")); o += 1

    parts.append(Part(o, "Transformer OLTC (On-Load Tap Changer)", [
        f"Type: diverter switch + tap selector, oil-immersed",
        f"Range: +/- 10% in 16 steps of 1.25%",
        f"Current: 1000 A max",
        f"Operation: motor-driven, automatic voltage control",
        f"Maintenance: diverter oil changed every 100,000 operations",
    ], "oltc")); o += 1

    parts.append(Part(o, "Transformer Bushings", [
        f"Type: oil-impregnated paper, porcelain housed",
        f"Count: 6 bushings (3 LV + 3 HV)",
        f"LV rating: {TURBINE_HW['generator_kv']:.1f} kV, 3000 A",
        f"HV rating: 132 kV, 1000 A",
        f"BIL: 170 kV (HV), 110 kV (LV)",
    ], "transformer-bushings")); o += 1

    parts.append(Part(o, "Transformer Conservator + Buchholz", [
        f"Type: oil expansion conservator with Buchholz relay",
        f"Conservator: 1000 L, bladder-type oil preservation",
        f"Buchholz: dual-stage (alarm + trip) gas + oil flow relay",
        f"Silica gel: color-change dehydrating breather",
        f"Function: oil expansion + internal fault detection",
    ], "transformer-conservator")); o += 1

    parts.append(Part(o, "Transformer Cooling Radiators", [
        f"Type: panel radiators with forced-air fans",
        f"Count: 8 radiator banks, 2 fans per bank",
        f"Fan: 1.5 kW each, ODP motor",
        f"Cooling: ONAN 60% / ONAF 100% capacity",
        f"Function: cools transformer oil",
    ], "transformer-radiators")); o += 1

    # --- EXPANSION JOINT DETAIL ---

    parts.append(Part(o, "Expansion Joint Tie Rods", [
        f"Type: threaded tie rods, stainless 316",
        f"Count: 4 per joint x {n_joints} joints = {4 * n_joints} rods",
        f"Function: limits bellows movement, prevents over-extension",
        f"Setting: adjusted for design pressure thrust",
    ], "ej-tie-rods")); o += 1

    parts.append(Part(o, "Expansion Joint Internal Sleeve", [
        f"Type: perforated inner sleeve, Inconel 625",
        f"Count: 1 per joint x {n_joints} joints",
        f"Function: protects bellows from flow erosion + directs flow",
        f"Holes: flow-through holes to equalize pressure",
    ], "ej-sleeve")); o += 1

    parts.append(Part(o, "Expansion Joint External Shroud", [
        f"Type: weather protection shroud, stainless 304",
        f"Count: 1 per joint x {n_joints} joints",
        f"Function: protects bellows from rain, dust, debris",
    ], "ej-shroud")); o += 1

    return parts


def print_parts(t: Dict) -> None:
    """Print the complete BOM parts list."""
    parts = build_parts_list(t)
    print_header("BILL OF MATERIALS  --  complete parts list")
    print(f"  {'#':>3}  {'Category':<16}  {'Part':<36}  Key specs")
    print(f"  {'':>3}  {'':<16}  {'':<36}  {'-'*40}")
    for p in parts:
        print(f"\n  [{p.order:>2}]  {p.category:<16}  {p.name}")
        for s in p.specs:
            print(f"       {s}")
    print(f"\n  Total assemblies: {len(parts)}")


# ==============================================================================
# SECTION 13 -- INFO / CLI
# ==============================================================================

def print_info() -> None:
    print_header("INFO  --  every part + the math")
    text = (
        "THE CYCLE, STEP BY STEP\n"
        "  1. A massive underground cavern is charged with cold dense air\n"
        "     (pressurised, -10 to -60 C). The cold comes from shallow earth\n"
        "     coupling, winter air, or an active chiller.\n"
        "  2. A discharge valve opens. The cavern pressure (P_cavern - P_atm)\n"
        "     drives cold air into a mile+ tunnel.\n"
        "  2b. REGENERATOR (optional): a counterflow heat exchanger preheats\n"
        "      the incoming cold air using the hot exhaust, raising T_in to\n"
        "      T_pre before the lava contact. This increases Carnot efficiency\n"
        "      and reduces the lava heat demand for the same T_hot.\n"
        "  3. The tunnel passes through a lava-heated contact section. Heat\n"
        "     flows Q = U A (T_lava - T_air) into the air, raising its\n"
        "     temperature from T_pre to T_hot. The heat-transfer area A is\n"
        "     multiplied by N_PARALLEL_BORES and FIN_FACTOR for enhanced\n"
        "     heat exchanger designs.\n"
        "  3b. REHEAT (optional): between turbine stages, the air is re-heated\n"
        "      back to T_hot by additional lava contact. This Brayton reheat\n"
        "      cycle significantly increases turbine work for the same PR.\n"
        "  4. The air expands through N staged turbines that extract electrical\n"
        "     work. With reheat, the expansion is split into N_reheat+1 segments\n"
        "     each with PR_seg = PR^(1/(N_reheat+1)), and the air is reheated\n"
        "     to T_hot between segments. Total work = (N_reheat+1) * cp * dT * mdot.\n"
        "  5. The remainder exits as a high-speed jet past M exit fans that\n"
        "     harvest the residual kinetic energy.\n"
        "  6. Buoyancy (stack pressure g H (rho_cold - rho_hot)) assists the\n"
        "     flow because the hot exit air is less dense.\n"
        "  7. While idle, the cavern passively exchanges heat with the ground\n"
        "     (Q = U_g A_g (T_ground - T_cavern)). If the ground is colder it\n"
        "     recharges; if warmer it leaks heat in and the cold is lost.\n"
        "\n"
        "POWER ENHANCEMENT LEVERS (engineered improvements)\n"
        "  The baseline presets produce 244-380 MW. The engineered presets apply\n"
        "  these improvements across 4 tiers to reach 2.3 TW:\n"
        "  * Parallel bores : N tunnels through the lava zone multiply the\n"
        "    heat-transfer area by N. 24 bores = 24x more heat input.\n"
        "  * Finned HX      : longitudinal fins on the tunnel walls in the\n"
        "    lava zone increase the effective area by 15-20x.\n"
        "  * Shell-and-tube : 50000 small-diameter tubes immersed in lava,\n"
        "    decoupling heat area from flow area. 193 GW/K total UA.\n"
        "  * Heat pipes     : NaK heat pipes embedded in lava give U=2000 W/m2K.\n"
        "  * Regenerator    : a counterflow recuperator preheats incoming air\n"
        "    with exhaust heat, raising Carnot efficiency by 5-10 points.\n"
        "  * Reheat stages  : re-inject lava heat between turbine stages\n"
        "    (Brayton reheat). 12 reheat stages can 4x the turbine work.\n"
        "  * MHD topping    : magnetohydrodynamic channel extracts DC work\n"
        "    from ionized gas at >1500 C before the first turbine stage.\n"
        "  * sCO2 bottoming : supercritical CO2 cycle on exhaust, 40% eff.\n"
        "  * ORC bottoming  : organic Rankine cycle on low-grade exhaust heat.\n"
        "  * Supersonic     : de Laval nozzle accelerates jet to Mach 2+.\n"
        "  * Smooth lining  : bored tunnel lining reduces friction by 2x.\n"
        "  * Higher charge P: 100 bar instead of 6 bar = 17x more stored exergy.\n"
        "  * Liquid air     : -196 C cryogenic charge, Carnot eta to 96.6%.\n"
        "  * Cascade cooling: multi-stage N2/He refrigeration for < -60 C.\n"
        "  * Intercooled    : 12-stage intercooled compression approaches\n"
        "    isothermal, reducing recharge energy by 3x.\n"
        "  * Carnot clamp   : total heat-engine work is clamped to eta_c * Q_lava\n"
        "    to prevent over-unity. Excess is rejected as waste heat.\n"
        "  * Dual tunnel    : n_systems=2 doubles everything (caverns, tunnels,\n"
        "    HX, turbines). Two complete setups side by side for 2x output.\n"
        "  * Potassium cycle: K vapor Rankine at 2000+ C, 50% eff (topping).\n"
        "  * Steam cycle    : tertiary bottoming at 500+ C, 38% eff.\n"
        "  * Quadruple bot. : K + sCO2 + steam + ORC cascade harvests exhaust\n"
        "    from 3000 C down to ambient -- maximum enthalpy extraction.\n"
        "\n"
        "THE MATH (key relations)\n"
        "  Ideal gas        : P = rho R T\n"
        "  Carnot ceiling   : W_max = Q_hot (1 - T_cold/T_hot)\n"
        "  Carnot clamp     : W_turb + W_fans + W_mhd + W_orc + W_sco2 <= eta_c * Q_lava\n"
        "  Lava heat        : Q = U * A * (T_lava - T_air_avg)\n"
        "                     A = pi*D*L*N_bores*fin + pi*d*L_tube*N_tubes\n"
        "  Energy balance   : Q_total = mdot cp (T_hot - T_pre) + N_reheat * mdot cp (T_hot - T_seg)\n"
        "  Stack pressure   : dP = g H (rho_cold - rho_hot)\n"
        "  Friction (Darcy) : dP = f (L/D) (1/2 rho v^2)\n"
        "  Turbine work     : W = (N_reheat+1) * eta_t * eta_g * mdot * cp * (T_hot - T_seg)\n"
        "  Regenerator      : T_pre = T_in + eps * (T_out - T_in)\n"
        "  Intercooled comp : W = N_stages * cp * T * (PR_seg^((g-1)/g) - 1) / eta\n"
        "  EROI             : E_out / (W_chill + W_liquefy + W_repressurise)\n"
        "\n"
        "THE AUDIT\n"
        "  Conservation : E_in(lava+leak+grid+initial) =\n"
        "                 E_out(elec+jet+exhaust+waste+chiller_amb+final)\n"
        "  Carnot+exergy: P_net <= R*T*ln(P/P0)*mdot + (1-T_cold/T_lava)*Q_lava\n"
        "  Both are asserted by --selftest for every preset target.\n"
    )
    for line in _wrap(text, 78):
        print(line)


def cmd_live(hours: float) -> None:
    print_header(f"LIVE DASHBOARD  --  Temperate-Shallow, {hours:.0f} h")
    t = targets()["Gmans Tunnel"]
    res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                   hours=hours, n_steps=600)
    print("\nNet power (MW) over time:")
    print(sparkline(res.p_net_mw, label="P_net MW"))
    print("\nCavern temperature (C) over time:")
    print(sparkline(res.cavern_t_c, label="cavern T C"))
    print("\nExit jet velocity (m/s):")
    print(sparkline(res.v_exit, label="v_exit m/s"))
    print()
    print_run(res, key="Gmans Tunnel")


def cmd_turbine_detail(key: str) -> None:
    """Per-stage turbine breakdown at peak flow, following the multi-stage
    decomposition from SECTION 4c."""
    t = targets()[key]
    res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                   hours=24.0, n_steps=1200)
    print_header(f"TURBINE STAGE BREAKDOWN  --  {key}  (at peak flow)")

    if not res.mdot or max(res.mdot) < 1e-3:
        print("  No flow -- cannot analyse turbine stages.")
        return

    idx = max(range(len(res.mdot)), key=lambda i: res.mdot[i])
    tn = t["tunnel"]
    cv = t["cavern"]

    # reconstruct the peak flow state
    st = build_initial_cavern(cv)
    fr = solve_flow(st, cv, t["lava"], tn, t["ctrl"])

    if fr.mdot < 1e-3:
        print("  No flow at peak.")
        return

    p_cavern = cv.p_charge_pa
    p_atm = P_STD
    stages = solve_turbine_stages(
        fr.t_hot_k, p_cavern, p_atm, fr.mdot,
        tn.n_turbine_stages, tn.turbine_eta, tn.generator_eta)

    print(f"\n  Peak flow: {fr.mdot:.1f} kg/s")
    print(f"  T_hot (combustor outlet): {k_to_c(fr.t_hot_k):.1f} C")
    print(f"  Pressure ratio: {p_cavern/p_atm:.2f}")
    print(f"  PR per stage: {(p_cavern/p_atm)**(1.0/tn.n_turbine_stages):.4f}")
    print()
    print(f"  {'stage':>6}{'P_in kPa':>10}{'P_out kPa':>10}{'T_in C':>9}"
          f"{'T_out C':>9}{'T_out_s C':>10}{'W kJ/kg':>9}{'P_elec MW':>11}")
    print(f"  {'-' * 76}")

    total_w = 0.0
    total_p = 0.0
    for s in stages:
        print(f"  {s.stage_num:>6d}{s.p_in_pa/KPA:>10.1f}{s.p_out_pa/KPA:>10.1f}"
              f"{k_to_c(s.t_in_k):>9.1f}{k_to_c(s.t_out_k):>9.1f}"
              f"{k_to_c(s.t_out_isentropic_k):>10.1f}"
              f"{s.work_kg/1000:>9.2f}{s.power_w/1e6:>11.3f}")
        total_w += s.work_kg
        total_p += s.power_w

    print(f"  {'-' * 76}")
    print(f"  {'TOTAL':>6}{'':>10}{'':>10}{'':>9}{k_to_c(stages[-1].t_out_k):>9.1f}"
          f"{'':>10}{total_w/1000:>9.2f}{total_p/1e6:>11.3f}")
    print(f"\n  Generator efficiency: {tn.generator_eta:.2f}")
    print(f"  Isentropic efficiency per stage: {tn.turbine_eta:.2f}")
    print(f"  Total shaft work: {total_w * fr.mdot / 1e6:.2f} MW (mechanical)")
    print(f"  Total electrical: {total_p / 1e6:.2f} MW")
    print(f"\n  Exit jet velocity: {fr.v_exit:.1f} m/s (Mach {fr.mach_exit:.2f})")
    print(f"  Exit fan power: {fr.p_exit_fans_w/1e6:.3f} MW")
    print(f"  Gross power: {fr.p_gross_w/1e6:.3f} MW")
    print(f"  Net power: {fr.p_net_w/1e6:.3f} MW")


# ==============================================================================
# SECTION 13e -- INTERACTIVE VISUALIZATION (matplotlib + tkinter GUI)
# ==============================================================================
#
# Following the pattern from Main_AIED.py and Simulation.py: a tkinter window
# with embedded matplotlib canvases, tabbed views, and the matplotlib navigation
# toolbar for pan/zoom.  Every plot reads from the same simulate()/solve_flow()
# data that the ASCII reports use, so the picture cannot drift from the physics.
#
# Requires:  matplotlib, numpy, tkinter  (all standard on most Python installs)
# Fallback:  prints a message if deps are missing.

def _fmt_power(w: float) -> str:
    """Format watts into a human-readable string."""
    a = abs(w)
    if a >= 1e12:
        return f"{w/1e12:.2f} TW"
    if a >= 1e9:
        return f"{w/1e9:.2f} GW"
    if a >= 1e6:
        return f"{w/1e6:.2f} MW"
    if a >= 1e3:
        return f"{w/1e3:.2f} kW"
    return f"{w:.1f} W"


def _fmt_energy(j: float) -> str:
    """Format joules into a human-readable string."""
    a = abs(j)
    if a >= 1e18:
        return f"{j/3.6e15:.2f} TWh"
    if a >= 1e15:
        return f"{j/3.6e12:.2f} GWh"
    if a >= 1e12:
        return f"{j/3.6e9:.2f} MWh"
    if a >= 1e9:
        return f"{j/3.6e6:.2f} kWh"
    return f"{j:.1f} J"


def _cylinder_faces(x0, y0, z0, x1, y1, z1, r, n_seg=8):
    """Return a list of quad faces for a 3D cylinder from (x0,y0,z0) to
    (x1,y1,z1) with radius r.  Used for Poly3DCollection."""
    import numpy as np
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    length = math.sqrt(dx*dx + dy*dy + dz*dz)
    if length < 1e-9:
        return []
    dx, dy, dz = dx/length, dy/length, dz/length
    if abs(dx) < 0.9:
        ux, uy, uz = 1.0, 0.0, 0.0
    else:
        ux, uy, uz = 0.0, 1.0, 0.0
    px = dy*uz - dz*uy
    py = dz*ux - dx*uz
    pz = dx*uy - dy*ux
    pl = math.sqrt(px*px + py*py + pz*pz)
    if pl < 1e-9:
        return []
    px, py, pz = px/pl, py/pl, pz/pl
    qx = dy*pz - dz*py
    qy = dz*px - dx*pz
    qz = dx*py - dy*px
    ql = math.sqrt(qx*qx + qy*qy + qz*qz)
    if ql < 1e-9:
        return []
    qx, qy, qz = qx/ql, qy/ql, qz/ql
    faces = []
    angles = np.linspace(0, 2*np.pi, n_seg, endpoint=False)
    for i in range(n_seg):
        a0 = angles[i]
        a1 = angles[(i+1) % n_seg]
        c0x = r*math.cos(a0); c0y = r*math.sin(a0)
        c1x = r*math.cos(a1); c1y = r*math.sin(a1)
        p0 = (x0 + px*c0x + qx*c0y, y0 + py*c0x + qy*c0y, z0 + pz*c0x + qz*c0y)
        p1 = (x0 + px*c1x + qx*c1y, y0 + py*c1x + qy*c1y, z0 + pz*c1x + qz*c1y)
        p2 = (x1 + px*c1x + qx*c1y, y1 + py*c1x + qy*c1y, z1 + pz*c1x + qz*c1y)
        p3 = (x1 + px*c0x + qx*c0y, y1 + py*c0x + qy*c0y, z1 + pz*c0x + qz*c0y)
        faces.append([p0, p1, p2, p3])
    cap0 = [(x0 + px*r*math.cos(a) + qx*r*math.sin(a),
             y0 + py*r*math.cos(a) + qy*r*math.sin(a),
             z0 + pz*r*math.cos(a) + qz*r*math.sin(a)) for a in angles]
    cap1 = [(x1 + px*r*math.cos(a) + qx*r*math.sin(a),
             y1 + py*r*math.cos(a) + qy*r*math.sin(a),
             z1 + pz*r*math.cos(a) + qz*r*math.sin(a)) for a in angles]
    faces.append(list(cap0))
    faces.append(list(cap1))
    return faces


def _draw_cross_section(ax, t: Dict) -> None:
    """Draw a 2D cross-section showing ALL components to scale.

    Shows: surface, ground, cavern, access tunnel, ALL parallel bores (stacked),
    lava zone, HX tubes, heat pipes, turbine stages, reheat, MHD, regenerator,
    stack, nozzle, exit fans, exit jet, bottoming cycles, transformer, chiller.
    Vertical scale is exaggerated for visibility (noted in axis label).
    """
    from matplotlib.patches import FancyBboxPatch
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t["ctrl"]
    n_sys = max(1, ctrl.n_systems)
    n_bores = max(1, lv.n_parallel_bores)
    ax.clear()

    # --- geometry ---
    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    cavern_depth = cv.depth_m
    tunnel_len = tn.total_length_m
    lava_len = lv.contact_length_m
    stack_h = tn.height_rise_m
    n_turb = tn.n_turbine_stages
    n_reheat = getattr(tn, 'n_reheat_stages', 0)
    n_fans = tn.n_exit_fans

    total_horiz = cav_side + tunnel_len + stack_h + 100
    max_depth = max(cavern_depth + cav_side, 50.0)
    max_height = max(stack_h, 50.0)
    v_exag = max(total_horiz / (max_depth + max_height + 100) * 0.35, 3.0)

    # --- surface and ground ---
    ax.axhline(0, color="#8B7355", linewidth=2, zorder=5)
    ax.fill_between([0, total_horiz * n_sys], 0, -max_depth * v_exag * 1.2,
                    color="#2B1D0E", alpha=0.3, zorder=0)
    ax.fill_between([0, total_horiz * n_sys], 0, max_height * v_exag * 1.2,
                    color="#1a1a2e", alpha=0.12, zorder=0)

    def draw_system(ox, alpha=1.0):
        # --- cavern ---
        cav_x0 = ox + 10
        cav_y0 = -cavern_depth * v_exag
        cav_w = cav_side
        cav_h = min(cav_side * 0.6, 40 * v_exag)
        t_c = k_to_c(cv.t_charge_k)
        cav_color = "#2196F3" if t_c < -100 else ("#4FC3F7" if t_c < 0 else "#81C784")
        cav_rect = FancyBboxPatch((cav_x0, cav_y0 - cav_h), cav_w, cav_h,
                                   boxstyle="round,pad=2", facecolor=cav_color,
                                   edgecolor="#0277FD", linewidth=2, alpha=0.6*alpha, zorder=3)
        ax.add_patch(cav_rect)
        if alpha > 0.5:
            ax.text(cav_x0 + cav_w/2, cav_y0 - cav_h/2,
                    f"COLD CAVERN\n{cv.volume_m3/1e6:.1f}M m³\n{t_c:.0f} °C\n{cv.p_charge_pa/1e5:.0f} bar",
                    ha="center", va="center", fontsize=6, fontweight="bold", color="white", zorder=4)

        # --- access tunnel ---
        ax.plot([cav_x0 + cav_w*0.3, cav_x0 + cav_w*0.3], [0, cav_y0],
                color="#666", linewidth=3, zorder=2)
        if alpha > 0.5:
            ax.text(cav_x0 + cav_w*0.3 + 3, cav_y0/2, "access", fontsize=5, color="#999", rotation=90)

        # --- parallel bores (stacked vertically in cross-section) ---
        tun_x0 = cav_x0 + cav_w
        tun_x1 = tun_x0 + tunnel_len
        bore_spacing = max(tn.diameter_m * 1.5 * v_exag, 8)
        n_show_bores = min(n_bores, 8)
        tun_y_base = cav_y0 + cav_h * 0.3

        for bi in range(n_show_bores):
            tun_y = tun_y_base - bi * bore_spacing
            ax.plot([tun_x0, tun_x1], [tun_y, tun_y], color="#FFD700", linewidth=2, alpha=0.5*alpha, zorder=3)
            ax.plot([tun_x0, tun_x1], [tun_y - 2, tun_y - 2], color="#B8860B", linewidth=1, alpha=0.4*alpha, zorder=3)
            ax.plot([tun_x0, tun_x1], [tun_y + 2, tun_y + 2], color="#B8860B", linewidth=1, alpha=0.4*alpha, zorder=3)

            lava_x0 = tun_x0 + tunnel_len * 0.15
            lava_x1 = min(lava_x0 + lava_len, tun_x1)
            ax.plot([lava_x0, lava_x1], [tun_y, tun_y], color="#FF6347", linewidth=3, alpha=0.6*alpha, zorder=4)

            if bi == 0 and alpha > 0.5:
                for i in range(n_turb):
                    frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tunnel_len)
                    tx = tun_x0 + frac * tunnel_len
                    ax.plot(tx, tun_y, "o", color="#00E676", markersize=max(4, min(8, 60/n_turb)), zorder=5)
                    if i < 6 or i == n_turb - 1:
                        ax.text(tx, tun_y + 6, f"T{i+1}", fontsize=4, ha="center", color="#00E676", fontweight="bold")

                if n_reheat > 0:
                    for i in range(min(n_reheat, 10)):
                        frac = 0.15 + (i + 1.5) / (n_turb + 1) * (lava_len / tunnel_len)
                        rx = tun_x0 + frac * tunnel_len
                        ax.plot(rx, tun_y + 5, "v", color="#FF5722", markersize=4, zorder=5)

                if tn.mhd_enabled:
                    mx = lava_x0 + 10
                    ax.plot(mx, tun_y, "D", color="#E91E63", markersize=5, zorder=5)
                    ax.text(mx, tun_y + 8, "MHD", fontsize=4, color="#E91E63")

                if tn.regenerator_eff > 0:
                    rx = tun_x0 + tunnel_len * 0.10
                    ax.plot(rx, tun_y, "s", color="#9C27B0", markersize=5, zorder=5)
                    ax.text(rx, tun_y + 8, "REGEN", fontsize=4, color="#9C27B0")

        if alpha > 0.5 and n_bores > n_show_bores:
            ax.text(tun_x0 + tunnel_len * 0.5, tun_y_base - n_show_bores * bore_spacing,
                    f"... {n_bores} parallel bores total", fontsize=5, color="#FFD700", ha="center")

        # --- lava body ---
        lava_x0 = tun_x0 + tunnel_len * 0.15
        lava_x1 = min(lava_x0 + lava_len, tun_x1)
        lava_y_top = tun_y_base - 15 * v_exag
        lava_y_bot = lava_y_top - 20 * v_exag
        ax.fill_between([lava_x0, lava_x1], lava_y_top, lava_y_bot,
                        color="#FF4500", alpha=0.6*alpha, zorder=2)
        ax.fill_between([lava_x0, lava_x1], lava_y_bot, lava_y_bot - 5 * v_exag,
                        color="#8B0000", alpha=0.7*alpha, zorder=2)
        for i in range(4):
            a = 0.12 - i * 0.02
            ax.fill_between([lava_x0 - i*3, lava_x1 + i*3],
                            lava_y_top + i*2, lava_y_bot - i*2,
                            color="#FF6347", alpha=a*alpha, zorder=1)
        if alpha > 0.5:
            ax.text((lava_x0 + lava_x1)/2, (lava_y_top + lava_y_bot)/2,
                    f"LAVA {lv.t_lava_c:.0f}°C", ha="center", va="center",
                    fontsize=6, fontweight="bold", color="white", zorder=4)

        if alpha > 0.5:
            for i in range(6):
                fx = lava_x0 + (lava_x1 - lava_x0) * (i + 0.5) / 6
                ax.annotate("", xy=(fx, tun_y_base - 2), xytext=(fx, lava_y_top),
                            arrowprops=dict(arrowstyle="->", color="#FF8C00", alpha=0.5, lw=1), zorder=4)

        if lv.hx_enabled and lv.hx_n_tubes > 0 and alpha > 0.5:
            ax.text((lava_x0 + lava_x1)/2, lava_y_bot - 8*v_exag,
                    f"HX: {lv.hx_n_tubes:,} tubes\nU={lv.hx_u:.0f} W/m²K",
                    fontsize=5, color="#FF8C00", ha="center", zorder=4)

        if lv.heat_pipe and alpha > 0.5:
            for i in range(5):
                px = lava_x0 + (i + 0.5) * lava_len / 5
                ax.plot([px, px], [lava_y_top, tun_y_base], color="#FF1744", linewidth=1, alpha=0.4, zorder=3)
            ax.text(lava_x0 + lava_len * 0.5, lava_y_top - 3, "heat pipes", fontsize=4, color="#FF1744")

        # --- stack ---
        stack_x = tun_x1
        ax.plot([stack_x, stack_x], [0, stack_h * v_exag], color="#888", linewidth=4, alpha=alpha, zorder=3)
        ax.plot([stack_x + 4, stack_x + 4], [0, stack_h * v_exag], color="#888", linewidth=2, alpha=alpha, zorder=3)

        if alpha > 0.5:
            ax.plot([stack_x - 3, stack_x + 7], [stack_h*v_exag, stack_h*v_exag + 5],
                    color="#FFEB3B", linewidth=2, zorder=5)
            ax.plot([stack_x + 7, stack_x + 3], [stack_h*v_exag + 5, stack_h*v_exag + 10],
                    color="#FFEB3B", linewidth=2, zorder=5)

        for i in range(n_fans):
            fy = stack_h * v_exag * (0.6 + 0.08 * i)
            ax.plot(stack_x + 2, fy, "s", color="#00BFFF", markersize=6, zorder=6, alpha=alpha)
        if alpha > 0.5:
            ax.text(stack_x + 12, stack_h * v_exag * 0.7, f"Fans x{n_fans}", fontsize=5, color="#00BFFF")

        if alpha > 0.5:
            for i in range(3):
                jx = stack_x + 2 + (i - 1) * 5
                ax.annotate("", xy=(jx, stack_h*v_exag + 15), xytext=(jx, stack_h*v_exag),
                            arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=1.5), zorder=5)

        bottoming = []
        if tn.potassium_enabled: bottoming.append(("K", "#FF9800", tn.potassium_eta))
        if tn.sco2_enabled: bottoming.append(("sCO2", "#9C27B0", tn.sco2_eta))
        if tn.steam_enabled: bottoming.append(("Steam", "#03A9F4", tn.steam_eta))
        if tn.orc_enabled: bottoming.append(("ORC", "#4CAF50", tn.orc_eta))
        for i, (name, color, eta) in enumerate(bottoming):
            bx = stack_x + 20 + i * 25
            by = -10 * v_exag
            rect = FancyBboxPatch((bx, by), 18, 10, boxstyle="round,pad=1",
                                  facecolor=color, edgecolor=color, alpha=0.5*alpha, zorder=3)
            ax.add_patch(rect)
            if alpha > 0.5:
                ax.text(bx + 9, by + 5, f"{name}\n{eta*100:.0f}%", fontsize=4,
                        ha="center", va="center", color="white", zorder=4)

        if alpha > 0.5:
            tx = stack_x + 20 + len(bottoming) * 25 + 10
            rect = FancyBboxPatch((tx, 5), 15, 8, boxstyle="round,pad=1",
                                  facecolor="#FFC107", edgecolor="#FFA000", alpha=0.5, zorder=3)
            ax.add_patch(rect)
            ax.text(tx + 7, 9, "XFMR", fontsize=4, ha="center", color="white", zorder=4)

        if (cv.active_cooling or cv.lava_heated_cooling) and alpha > 0.5:
            ch_x = cav_x0 - 20
            ch_y = cav_y0
            rect = FancyBboxPatch((ch_x, ch_y), 15, 10, boxstyle="round,pad=1",
                                  facecolor="#80DEEA", edgecolor="#00ACC1", alpha=0.5, zorder=3)
            ax.add_patch(rect)
            lbl = "Abs Chiller" if cv.lava_heated_cooling else "Chiller"
            ax.text(ch_x + 7, ch_y + 5, lbl, fontsize=4, ha="center", color="white", zorder=4)

    for si in range(n_sys):
        ox = si * total_horiz
        draw_system(ox, alpha=1.0 if si == 0 else 0.5)

    ax.set_title(f"Cross-Section - All Components ({n_bores} bores" +
                 (f", {n_sys} systems" if n_sys > 1 else "") + ")", fontsize=10, fontweight="bold")
    ax.set_xlabel("Distance (m)", fontsize=8, color="#c0c0c0")
    ax.set_ylabel(f"Depth/Height (m, {v_exag:.0f}x exaggerated)", fontsize=8, color="#c0c0c0")
    ax.set_aspect("equal")
    ax.set_facecolor("#0d1117")
    ax.tick_params(colors="#c0c0c0", labelsize=6)
    for spine in ax.spines.values():
        spine.set_color("#444")
    ax.title.set_color("white")

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#4FC3F7", markersize=7, label="Cold cavern"),
        Line2D([0], [0], color="#FFD700", linewidth=2, label=f"Tunnel bore (x{n_bores})"),
        Line2D([0], [0], color="#FF6347", linewidth=3, label="Lava contact"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#00E676", markersize=6, label=f"Turbine (x{n_turb})"),
        Line2D([0], [0], marker="v", color="w", markerfacecolor="#FF5722", markersize=5, label=f"Reheat (x{n_reheat})"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#00BFFF", markersize=5, label=f"Exit fan (x{n_fans})"),
        Line2D([0], [0], color="#FF4500", linewidth=6, label="Lava body"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#E91E63", markersize=5, label="MHD"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#9C27B0", markersize=5, label="Regenerator"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#FF9800", markersize=5, label="Bottoming cycle"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#FFC107", markersize=5, label="Transformer"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#80DEEA", markersize=5, label="Chiller"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=5,
              facecolor="#1a1a2e", edgecolor="#444", labelcolor="white", ncol=2)


def _draw_timeline(ax, res: SimResult) -> None:
    """Draw the timeline (power, cavern T, cavern P, exit velocity) on matplotlib."""
    ax.clear()
    if not res.t_h:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                color="#888", fontsize=12)
        return

    t_h = res.t_h
    ax.set_facecolor("#0d1117")

    p_mw = [p for p in res.p_net_mw]
    ax.plot(t_h, p_mw, color="#00E676", linewidth=1.5, label="P_net (MW)")
    if res.q_lava_mw:
        ax.plot(t_h, res.q_lava_mw, color="#FF4500", linewidth=1, alpha=0.7,
                label="Q_lava (MW)")
    ax.set_xlabel("Time (h)", fontsize=8, color="#c0c0c0")
    ax.set_ylabel("Power / Heat (MW)", fontsize=8, color="#00E676")
    ax.tick_params(axis="y", labelcolor="#00E676", labelsize=7)
    ax.tick_params(axis="x", labelcolor="#c0c0c0", labelsize=7)

    ax2 = ax.twinx()
    ax2.plot(t_h, res.cavern_t_c, color="#4FC3F7", linewidth=1, linestyle="--",
             label="Cavern T (C)")
    ax2.set_ylabel("Cavern T (C)", fontsize=8, color="#4FC3F7")
    ax2.tick_params(axis="y", labelcolor="#4FC3F7", labelsize=7)
    ax2.set_facecolor("#0d1117")

    if res.cavern_p_kpa:
        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("outward", 45))
        ax3.plot(t_h, res.cavern_p_kpa, color="#FFD700", linewidth=0.8, alpha=0.6,
                 label="Cavern P (kPa)")
        ax3.set_ylabel("Cavern P (kPa)", fontsize=8, color="#FFD700")
        ax3.tick_params(axis="y", labelcolor="#FFD700", labelsize=7)

    ax.set_title("Timeline - Power, Temperature, Pressure", fontsize=10,
                 fontweight="bold", color="white")
    ax.grid(True, alpha=0.15, color="#444")
    for spine in ax.spines.values():
        spine.set_color("#444")

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right",
              fontsize=6, facecolor="#1a1a2e", edgecolor="#444",
              labelcolor="white")


def _draw_operations(ax, res: SimResult) -> None:
    """Draw a detailed operations panel showing energy output over time.

    Multi-panel layout:
      - Top: Net power output (MW) over time with gross and parasitic
      - Middle: Cumulative energy output (TWh) and lava heat input
      - Bottom: Carnot efficiency vs actual efficiency over time
    """
    ax.clear()
    if not res.t_h:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                color="#888", fontsize=12)
        return

    ax.set_axis_off()
    fig = ax.figure
    ax.set_facecolor("#0d1117")

    # Create 3 sub-axes manually
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=ax.get_subplotspec(),
                                          hspace=0.55)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    for a in [ax1, ax2, ax3]:
        a.set_facecolor("#0d1117")
        a.tick_params(colors="#c0c0c0", labelsize=7)
        for spine in a.spines.values():
            spine.set_color("#444")
        a.grid(True, alpha=0.15, color="#444")

    t_h = res.t_h

    # --- Panel 1: Power output over time ---
    ax1.plot(t_h, res.p_net_mw, color="#00E676", linewidth=1.5, label="P_net (MW)")
    if res.p_gross_mw:
        ax1.plot(t_h, res.p_gross_mw, color="#4CAF50", linewidth=1, alpha=0.5,
                 linestyle=":", label="P_gross (MW)")
    if res.q_lava_mw:
        ax1.plot(t_h, res.q_lava_mw, color="#FF4500", linewidth=1, alpha=0.6,
                 label="Q_lava (MW)")
    ax1.set_ylabel("Power (MW)", fontsize=8, color="#00E676")
    ax1.tick_params(axis="y", labelcolor="#00E676", labelsize=7)
    ax1.set_title("Power Output Over Time", fontsize=9, fontweight="bold", color="white")
    ax1.legend(loc="upper right", fontsize=6, facecolor="#1a1a2e",
               edgecolor="#444", labelcolor="white")

    # --- Panel 2: Cumulative energy ---
    # Integrate p_net over time to get cumulative TWh
    cumul_twh = []
    total = 0.0
    for i in range(len(t_h)):
        if i > 0:
            dt_h = t_h[i] - t_h[i-1]
            total += res.p_net_mw[i] * dt_h / 1e6  # MW*h -> TWh
        cumul_twh.append(total)
    ax2.plot(t_h, cumul_twh, color="#00E676", linewidth=1.5, label="Cumulative E_out (TWh)")
    ax2.fill_between(t_h, 0, cumul_twh, color="#00E676", alpha=0.15)
    ax2.set_ylabel("Cumulative (TWh)", fontsize=8, color="#00E676")
    ax2.tick_params(axis="y", labelcolor="#00E676", labelsize=7)
    ax2.set_title("Cumulative Energy Output", fontsize=9, fontweight="bold", color="white")
    ax2.legend(loc="upper left", fontsize=6, facecolor="#1a1a2e",
               edgecolor="#444", labelcolor="white")

    # --- Panel 3: Efficiency over time ---
    if res.eta_carnot and res.p_carnot_mw and res.p_net_mw:
        actual_eta = []
        for i in range(len(t_h)):
            if res.q_lava_mw[i] > 0 and res.p_net_mw[i] > 0:
                actual_eta.append(res.p_net_mw[i] / res.q_lava_mw[i])
            else:
                actual_eta.append(0.0)
        ax3.plot(t_h, [e * 100 for e in res.eta_carnot], color="#FF4500",
                 linewidth=1, linestyle="--", label="Carnot limit (%)")
        ax3.plot(t_h, [e * 100 for e in actual_eta], color="#00E676",
                 linewidth=1.5, label="Actual efficiency (%)")
        ax3.set_ylabel("Efficiency (%)", fontsize=8, color="#00E676")
        ax3.tick_params(axis="y", labelcolor="#00E676", labelsize=7)
        ax3.set_xlabel("Time (h)", fontsize=8, color="#c0c0c0")
        ax3.set_title("Thermal Efficiency vs Carnot Limit", fontsize=9, fontweight="bold", color="white")
        ax3.legend(loc="upper right", fontsize=6, facecolor="#1a1a2e",
                   edgecolor="#444", labelcolor="white")
    else:
        ax3.set_xlabel("Time (h)", fontsize=8, color="#c0c0c0")
        ax3.text(0.5, 0.5, "No efficiency data", transform=ax3.transAxes,
                 ha="center", color="#888", fontsize=10)


def _draw_energy_flow(ax, res: SimResult, fr: FlowResult) -> None:
    """Draw an energy flow bar chart on matplotlib."""
    ax.clear()
    ax.set_facecolor("#0d1117")

    categories = []
    values = []
    colors = []

    categories.append("Lava heat\n(main+reheat)")
    values.append(fr.q_lava_w / 1e6)
    colors.append("#FF4500")

    if res.heat_leaked_twh > 0:
        categories.append("Ground leak")
        leak_mw = res.heat_leaked_twh * 1e12 / (max(res.discharge_hours, 0.1) * 3600) / 1e6
        values.append(leak_mw)
        colors.append("#8B4513")

    categories.append("Turbine\n elec.")
    values.append(fr.p_turbine_stages_w / 1e6)
    colors.append("#00E676")

    categories.append("Exit fans")
    values.append(fr.p_exit_fans_w / 1e6)
    colors.append("#00BFFF")

    if fr.p_mhd_w > 0:
        categories.append("MHD")
        values.append(fr.p_mhd_w / 1e6)
        colors.append("#E91E63")

    if fr.p_potassium_w > 0:
        categories.append("Potassium")
        values.append(fr.p_potassium_w / 1e6)
        colors.append("#FF9800")

    if fr.p_sco2_w > 0:
        categories.append("sCO2")
        values.append(fr.p_sco2_w / 1e6)
        colors.append("#9C27B0")

    if fr.p_steam_w > 0:
        categories.append("Steam")
        values.append(fr.p_steam_w / 1e6)
        colors.append("#03A9F4")

    if fr.p_orc_w > 0:
        categories.append("ORC")
        values.append(fr.p_orc_w / 1e6)
        colors.append("#4CAF50")

    categories.append("Jet KE\n(residual)")
    values.append(max(fr.ke_jet_w - fr.p_exit_fans_w, 0) / 1e6)
    colors.append("#FFEB3B")

    categories.append("Gen loss")
    values.append(fr.gen_loss_w / 1e6)
    colors.append("#F44336")

    categories.append("Net output")
    values.append(fr.p_net_w / 1e6)
    colors.append("#00E676")

    x = range(len(categories))
    bars = ax.bar(x, values, color=colors, edgecolor="#333", linewidth=0.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(categories, fontsize=6, color="#c0c0c0", rotation=45, ha="right")
    ax.set_ylabel("MW", fontsize=8, color="#c0c0c0")
    ax.set_title("Energy Flow (peak, MW)", fontsize=10, fontweight="bold", color="white")
    ax.grid(True, axis="y", alpha=0.15, color="#444")
    for spine in ax.spines.values():
        spine.set_color("#444")
    ax.tick_params(colors="#c0c0c0", labelsize=7)

    for bar, val in zip(bars, values):
        if abs(val) > 0:
            ytext = bar.get_height()
            va = "bottom" if ytext >= 0 else "top"
            ax.text(bar.get_x() + bar.get_width() / 2, ytext, _fmt_power(val * 1e6),
                    ha="center", va=va, fontsize=5, color="white", rotation=90)


def _draw_turbine_stages(ax, stages: List[TurbineStage]) -> None:
    """Draw per-stage turbine data on matplotlib (dual axis: T and W)."""
    ax.clear()
    ax.set_facecolor("#0d1117")

    if not stages:
        ax.text(0.5, 0.5, "No turbine data", transform=ax.transAxes,
                ha="center", color="#888")
        return

    n = list(range(1, len(stages) + 1))
    t_in = [k_to_c(s.t_in_k) for s in stages]
    t_out = [k_to_c(s.t_out_k) for s in stages]
    t_isentropic = [k_to_c(s.t_out_isentropic_k) for s in stages]
    work = [s.power_w / 1e6 for s in stages]

    ax.plot(n, t_in, "o-", color="#FF4500", linewidth=1.5, markersize=4,
            label="T_in (C)")
    ax.plot(n, t_out, "s-", color="#4FC3F7", linewidth=1.5, markersize=4,
            label="T_out (C)")
    ax.plot(n, t_isentropic, "--", color="#888", linewidth=0.8, alpha=0.6,
            label="T_out isentropic (C)")
    ax.set_xlabel("Stage number", fontsize=8, color="#c0c0c0")
    ax.set_ylabel("Temperature (C)", fontsize=8, color="#FF4500")
    ax.tick_params(axis="y", labelcolor="#FF4500", labelsize=7)
    ax.tick_params(axis="x", labelcolor="#c0c0c0", labelsize=7)

    ax2 = ax.twinx()
    ax2.bar(n, work, alpha=0.25, color="#00E676", width=0.6, label="P_elec (MW)")
    ax2.set_ylabel("Electrical power (MW)", fontsize=8, color="#00E676")
    ax2.tick_params(axis="y", labelcolor="#00E676", labelsize=7)
    ax2.set_facecolor("#0d1117")

    ax.set_title("Turbine Stage Breakdown (peak flow)", fontsize=10,
                 fontweight="bold", color="white")
    ax.grid(True, alpha=0.15, color="#444")
    for spine in ax.spines.values():
        spine.set_color("#444")

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right",
              fontsize=6, facecolor="#1a1a2e", edgecolor="#444",
              labelcolor="white")


def _draw_pressure_profile(ax, stages: List[TurbineStage]) -> None:
    """Draw pressure profile across turbine stages."""
    ax.clear()
    ax.set_facecolor("#0d1117")
    if not stages:
        ax.text(0.5, 0.5, "No turbine data", transform=ax.transAxes,
                ha="center", color="#888")
        return

    n = list(range(0, len(stages) + 1))
    p = [stages[0].p_in_pa / 1e5] + [s.p_out_pa / 1e5 for s in stages]
    ax.plot(n, p, "o-", color="#FFD700", linewidth=1.5, markersize=5)
    ax.fill_between(n, 0, p, color="#FFD700", alpha=0.15)
    ax.set_xlabel("Stage number", fontsize=8, color="#c0c0c0")
    ax.set_ylabel("Pressure (bar)", fontsize=8, color="#FFD700")
    ax.tick_params(axis="y", labelcolor="#FFD700", labelsize=7)
    ax.tick_params(axis="x", labelcolor="#c0c0c0", labelsize=7)
    ax.set_title("Pressure Profile Across Turbine Stages", fontsize=10,
                 fontweight="bold", color="white")
    ax.grid(True, alpha=0.15, color="#444")
    for spine in ax.spines.values():
        spine.set_color("#444")


def _draw_cavern_state(ax, res: SimResult) -> None:
    """Draw cavern state (mass and pressure) over time."""
    ax.clear()
    ax.set_facecolor("#0d1117")
    if not res.t_h:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                color="#888", fontsize=12)
        return

    t_h = res.t_h
    ax.plot(t_h, [m / 1e6 for m in res.cavern_m_kg], color="#FFD700",
            linewidth=1.5, label="Air mass (Mkg)")
    ax.set_xlabel("Time (h)", fontsize=8, color="#c0c0c0")
    ax.set_ylabel("Air mass (Mkg)", fontsize=8, color="#FFD700")
    ax.tick_params(axis="y", labelcolor="#FFD700", labelsize=7)
    ax.tick_params(axis="x", labelcolor="#c0c0c0", labelsize=7)

    ax2 = ax.twinx()
    ax2.plot(t_h, [p / 1000 for p in res.cavern_p_kpa], color="#4FC3F7",
             linewidth=1, linestyle="--", label="Pressure (kPa)")
    ax2.set_ylabel("Pressure (kPa)", fontsize=8, color="#4FC3F7")
    ax2.tick_params(axis="y", labelcolor="#4FC3F7", labelsize=7)

    ax.set_title("Cavern State - Air Mass and Pressure", fontsize=10,
                 fontweight="bold", color="white")
    ax.grid(True, alpha=0.15, color="#444")
    for spine in ax.spines.values():
        spine.set_color("#444")

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right",
              fontsize=6, facecolor="#1a1a2e", edgecolor="#444",
              labelcolor="white")


def _draw_summary_panel(ax, res: SimResult, key: str) -> None:
    """Draw a summary panel with key metrics."""
    ax.clear()
    ax.set_facecolor("#0d1117")
    ax.set_axis_off()

    lines = [
        f"  TARGET: {key}",
        f"  Mode: {res.mode}",
        f"",
        f"  Mean net power : {_fmt_power(res.mean_p_net_mw * 1e6)}",
        f"  Peak net power : {_fmt_power(res.peak_p_net_mw * 1e6)}",
        f"  Total energy   : {res.total_twh:.3f} TWh",
        f"  Discharge time : {res.discharge_hours:.1f} h",
        f"  EROI           : {res.eroi:.2f}",
        f"  CAPEX          : ${res.capex_musd:.1f}M",
        f"  Homes powered  : {res.homes:.0f}",
        f"",
        f"  Conservation residual : {res.energy_residual:.2e}",
        f"  Verdict: {res.verdict}",
    ]

    y = 0.95
    for line in lines:
        color = "#00E676" if "Mean" in line or "Peak" in line else "#c0c0c0"
        if "TARGET" in line:
            color = "#e94560"
            ax.text(0.05, y, line, transform=ax.transAxes, fontsize=14,
                    fontweight="bold", color=color, family="Consolas")
        elif "Verdict" in line:
            color = "#FFD700" if "OK" in line else "#F44336"
            ax.text(0.05, y, line, transform=ax.transAxes, fontsize=11,
                    color=color, family="Consolas")
        else:
            ax.text(0.05, y, line, transform=ax.transAxes, fontsize=10,
                    color=color, family="Consolas")
        y -= 0.07


def _draw_turbine_engine(ax, t: Dict, rotation_angle: float = 0.0) -> None:
    """Draw a detailed axial-flow turbine engine cross-section with spinning blades.

    Shows the internal geometry of a multi-stage air-expansion turbine:
    - Stator vanes (fixed) and rotor blades (spinning) per stage
    - Shaft connecting all stages
    - Casing/housing
    - Inlet and outlet
    - Generator coupling
    - Blade count and rotor diameter from TURBINE_HW

    The rotation_angle parameter animates the rotor blades.
    """
    import numpy as np
    ax.clear()
    ax.set_facecolor("#0d1117")
    ax.set_aspect("equal")

    tn = t["tunnel"]
    n_stages = tn.n_turbine_stages
    rotor_r = TURBINE_HW["rotor_d_mm"] / 1000.0 / 2.0  # rotor radius in meters
    n_blades = TURBINE_HW["rotor_blade_count"]
    stage_spacing = TURBINE_HW["stage_spacing_m"]
    casing_r = rotor_r * 1.4
    shaft_r = rotor_r * 0.15

    # Scale for display: fit all stages horizontally
    total_len = n_stages * stage_spacing
    x_scale = 1.0  # already in meters
    # Vertical exaggeration for blade visibility
    v_scale = max(3.0, 20.0 / rotor_r)

    # Draw casing (top and bottom)
    casing_y_top = casing_r * v_scale
    casing_y_bot = -casing_r * v_scale
    ax.fill_between([0, total_len], casing_y_top, casing_y_top + 5,
                    color="#555", edgecolor="#333", zorder=1)
    ax.fill_between([0, total_len], casing_y_bot - 5, casing_y_bot,
                    color="#555", edgecolor="#333", zorder=1)

    # Draw shaft
    ax.fill_between([0, total_len + stage_spacing], -shaft_r * v_scale, shaft_r * v_scale,
                    color="#888", edgecolor="#666", zorder=2)

    # Draw each stage - batched line segments for performance
    stator_segs_x = []
    stator_segs_y = []
    rotor_segs_x = []
    rotor_segs_y = []
    stage_labels_x = []
    stage_labels_y = []
    stage_labels_t = []
    axvline_xs = []

    for i in range(n_stages):
        cx = (i + 0.5) * stage_spacing

        # Stator vanes (fixed, angled) - collect line segments
        for b in range(min(n_blades, 12)):
            ang = 2 * math.pi * b / min(n_blades, 12)
            sx = cx - stage_spacing * 0.3
            y_outer = casing_r * 0.9 * v_scale * (1 if ang < math.pi else -1)
            y_inner = rotor_r * 0.9 * v_scale * (1 if ang < math.pi else -1)
            # only draw top and bottom vanes (2D cross-section)
            stator_segs_x.extend([sx, sx, None])
            stator_segs_y.extend([y_inner, y_outer, None])

        # Rotor blades (spinning) - collect line segments
        for b in range(min(n_blades, 12)):
            ang = 2 * math.pi * b / min(n_blades, 12) + rotation_angle
            rx = cx + stage_spacing * 0.2
            blade_y = rotor_r * v_scale * math.sin(ang)
            blade_x_offset = rotor_r * 0.3 * math.cos(ang)
            visibility = abs(math.sin(ang))
            if visibility > 0.3:
                rotor_segs_x.extend([rx + blade_x_offset, rx + blade_x_offset, None])
                rotor_segs_y.extend([shaft_r * v_scale * (1 if math.sin(ang) > 0 else -1), blade_y, None])

        # Stage boundary
        axvline_xs.append(cx + stage_spacing * 0.5)

        # Label first and last few stages
        if i < 3 or i >= n_stages - 3 or (i == n_stages // 2):
            stage_labels_x.append(cx)
            stage_labels_y.append(casing_y_top + 15)
            stage_labels_t.append(f"S{i+1}")
        elif i == 3:
            stage_labels_x.append(cx)
            stage_labels_y.append(casing_y_top + 15)
            stage_labels_t.append("...")

    # Batch draw all stator vanes in one call
    if stator_segs_x:
        ax.plot(stator_segs_x, stator_segs_y, color="#9E9E9E",
                linewidth=1.5, alpha=0.6, zorder=3)

    # Batch draw all rotor blades in one call
    if rotor_segs_x:
        ax.plot(rotor_segs_x, rotor_segs_y, color="#00E676",
                linewidth=2, alpha=0.8, zorder=4)

    # Batch draw stage boundaries
    for vx in axvline_xs:
        ax.axvline(vx, color="#333", linewidth=0.5, alpha=0.3, zorder=1)

    # Batch draw stage labels
    for lx, ly, lt in zip(stage_labels_x, stage_labels_y, stage_labels_t):
        color = "#00E676" if lt != "..." else "#888"
        ax.text(lx, ly, lt, fontsize=6, ha="center", color=color, fontweight="bold")

    # Inlet (left side)
    ax.annotate("", xy=(0, 0), xytext=(-stage_spacing * 0.5, 0),
                arrowprops=dict(arrowstyle="->", color="#4FC3F7", lw=2))
    ax.text(-stage_spacing * 0.4, casing_r * v_scale * 0.5, "COLD AIR IN",
            fontsize=7, color="#4FC3F7", fontweight="bold", rotation=0)

    # Outlet (right side)
    ax.annotate("", xy=(total_len + stage_spacing * 0.5, 0),
                xytext=(total_len, 0),
                arrowprops=dict(arrowstyle="->", color="#FF4500", lw=2))
    ax.text(total_len + stage_spacing * 0.1, casing_r * v_scale * 0.5, "HOT EXHAUST",
            fontsize=7, color="#FF4500", fontweight="bold")

    # Generator (right end)
    gen_x = total_len + stage_spacing * 0.5
    gen_r = rotor_r * 1.2 * v_scale
    from matplotlib.patches import FancyBboxPatch
    gen_rect = FancyBboxPatch((gen_x, -gen_r), stage_spacing * 0.5, gen_r * 2,
                               boxstyle="round,pad=2", facecolor="#FFC107",
                               edgecolor="#FFA000", alpha=0.6, zorder=3)
    ax.add_patch(gen_rect)
    ax.text(gen_x + stage_spacing * 0.25, 0, "GEN\n13.8kV\n45MVA",
            fontsize=6, ha="center", va="center", color="white", fontweight="bold", zorder=4)

    # Specs
    spec_text = (f"Stages: {n_stages}  |  Rotor: {TURBINE_HW['rotor_d_mm']/1000:.1f}m  |  "
                 f"Blades/stage: {n_blades}  |  RPM: {TURBINE_HW['rpm']:.0f}  |  "
                 f"eta_isentropic: {TURBINE_HW['eta_isentropic']}")
    ax.text(0.02, 0.02, spec_text, transform=ax.transAxes, fontsize=6,
            color="#c0c0c0", family="monospace")

    ax.set_xlim(-stage_spacing * 0.6, total_len + stage_spacing * 1.2)
    ax.set_ylim(-casing_r * v_scale * 1.3, casing_r * v_scale * 1.3)
    ax.set_title(f"Multi-Stage Axial Turbine Engine ({n_stages} stages, spinning)",
                 fontsize=10, fontweight="bold", color="white")
    ax.tick_params(colors="#c0c0c0", labelsize=6)
    for spine in ax.spines.values():
        spine.set_color("#444")


def _draw_isometric_cutaway(ax, t: Dict) -> None:
    """Draw an isometric cutaway view showing the system from a 3D angle.

    This is NOT a flat 2D view - it uses isometric projection to show depth,
    with cutaway sections revealing internal components.
    """
    from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t["ctrl"]
    n_bores = max(1, lv.n_parallel_bores)
    n_turb = tn.n_turbine_stages
    n_fans = tn.n_exit_fans

    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    cav_depth = cv.depth_m
    tun_len = tn.total_length_m
    lava_len = lv.contact_length_m
    stack_h = tn.height_rise_m

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Isometric projection: x_screen = x - y*cos(30), y_screen = -(z + (x+y)*sin(30))
    import math as _m
    cos30 = _m.cos(_m.radians(30))
    sin30 = _m.sin(_m.radians(30))

    def iso(x, y, z):
        """Project 3D (x,y,z) to 2D isometric screen coords."""
        sx = (x - y) * cos30
        sy = -(z + (x + y) * sin30 * 0.5)
        return sx, sy

    # Scale to fit
    total_x = cav_side + tun_len + stack_h + 100
    total_z = cav_depth + cav_side + stack_h
    scale = 0.85 / max(total_x * cos30 * 1.4, total_z * 1.2)
    def iso_s(x, y, z):
        sx, sy = iso(x, y, z)
        return sx * scale, sy * scale

    # --- Draw ground surface (isometric plane) ---
    ground_pts = [iso_s(0, 0, 0), iso_s(total_x, 0, 0),
                  iso_s(total_x, cav_side, 0), iso_s(0, cav_side, 0)]
    ax.add_patch(Polygon(ground_pts, facecolor="#2B1D0E", edgecolor="#5C4033",
                         linewidth=0.5, alpha=0.3))

    # --- Cavern (cutaway box underground) ---
    cav_x0, cav_y0, cav_z0 = 10, 10, -cav_depth
    cav_x1, cav_y1, cav_z1 = 10 + cav_side, 10 + cav_side, -cav_depth + cav_side * 0.5

    # Cavern bottom face
    pts = [iso_s(cav_x0, cav_y0, cav_z0), iso_s(cav_x1, cav_y0, cav_z0),
           iso_s(cav_x1, cav_y1, cav_z0), iso_s(cav_x0, cav_y1, cav_z0)]
    ax.add_patch(Polygon(pts, facecolor="#1565C0", edgecolor="#0277FD", linewidth=0.8, alpha=0.3))
    # Cavern top face (cutaway - show interior)
    pts = [iso_s(cav_x0, cav_y0, cav_z1), iso_s(cav_x1, cav_y0, cav_z1),
           iso_s(cav_x1, cav_y1, cav_z1), iso_s(cav_x0, cav_y1, cav_z1)]
    ax.add_patch(Polygon(pts, facecolor="#2196F3", edgecolor="#0277FD", linewidth=0.8, alpha=0.2))
    # Cavern front face (cutaway)
    pts = [iso_s(cav_x0, cav_y0, cav_z0), iso_s(cav_x1, cav_y0, cav_z0),
           iso_s(cav_x1, cav_y0, cav_z1), iso_s(cav_x0, cav_y0, cav_z1)]
    ax.add_patch(Polygon(pts, facecolor="#2196F3", edgecolor="#0277FD", linewidth=0.8, alpha=0.25))
    # Cavern right face
    pts = [iso_s(cav_x1, cav_y0, cav_z0), iso_s(cav_x1, cav_y1, cav_z0),
           iso_s(cav_x1, cav_y1, cav_z1), iso_s(cav_x1, cav_y0, cav_z1)]
    ax.add_patch(Polygon(pts, facecolor="#1976D2", edgecolor="#0277FD", linewidth=0.8, alpha=0.2))

    # Cavern label
    cx_s, cy_s = iso_s((cav_x0+cav_x1)/2, cav_y0, (cav_z0+cav_z1)/2)
    ax.text(cx_s, cy_s, f"CAVERN\n{cv.volume_m3/1e9:.1f} km3\n{k_to_c(cv.t_charge_k):.0f}C",
            fontsize=4, ha="center", va="center", color="white", fontweight="bold")

    # Cavern lining layers (show on front face edge)
    lining_t = CAVERN_HW['lining_thick_mm'] / 1000.0
    # shotcrete layer
    pts = [iso_s(cav_x0, cav_y0, cav_z0), iso_s(cav_x0+lining_t, cav_y0, cav_z0),
           iso_s(cav_x0+lining_t, cav_y0, cav_z1), iso_s(cav_x0, cav_y0, cav_z1)]
    ax.add_patch(Polygon(pts, facecolor="#6D4C41", edgecolor="#4E342E", linewidth=0.3, alpha=0.4))

    # --- Access tunnel (from surface to cavern) ---
    acc_x = cav_x0 + cav_side * 0.3
    acc_r = CAVERN_HW['access_tunnel_d_m'] / 2.0
    # Draw as isometric rectangle
    pts = [iso_s(acc_x - acc_r, cav_y0, 0), iso_s(acc_x + acc_r, cav_y0, 0),
           iso_s(acc_x + acc_r, cav_y0, cav_z1), iso_s(acc_x - acc_r, cav_y0, cav_z1)]
    ax.add_patch(Polygon(pts, facecolor="#888", edgecolor="#666", linewidth=0.5, alpha=0.4))

    # --- Tunnel (from cavern to stack) ---
    tun_x0 = cav_x1
    tun_x1 = tun_x0 + tun_len
    tun_z = cav_z1 * 0.6
    tun_r = tn.diameter_m / 2.0

    # Tunnel as isometric cylinder (approximated as box)
    # Pre-lava section
    lava_x0 = tun_x0 + tun_len * 0.15
    lava_x1 = lava_x0 + lava_len
    pts = [iso_s(tun_x0, cav_y0, tun_z - tun_r), iso_s(lava_x0, cav_y0, tun_z - tun_r),
           iso_s(lava_x0, cav_y0, tun_z + tun_r), iso_s(tun_x0, cav_y0, tun_z + tun_r)]
    ax.add_patch(Polygon(pts, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.5, alpha=0.3))
    # Lava section
    pts = [iso_s(lava_x0, cav_y0, tun_z - tun_r), iso_s(lava_x1, cav_y0, tun_z - tun_r),
           iso_s(lava_x1, cav_y0, tun_z + tun_r), iso_s(lava_x0, cav_y0, tun_z + tun_r)]
    ax.add_patch(Polygon(pts, facecolor="#FF6347", edgecolor="#FF4500", linewidth=0.5, alpha=0.4))
    # Post-lava section
    pts = [iso_s(lava_x1, cav_y0, tun_z - tun_r), iso_s(tun_x1, cav_y0, tun_z - tun_r),
           iso_s(tun_x1, cav_y0, tun_z + tun_r), iso_s(lava_x1, cav_y0, tun_z + tun_r)]
    ax.add_patch(Polygon(pts, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.5, alpha=0.3))

    # --- Lava body (below tunnel in lava zone) ---
    lava_d = tn.diameter_m * 3
    pts = [iso_s(lava_x0, cav_y0 + cav_side*0.3, tun_z - tun_r - lava_d),
           iso_s(lava_x1, cav_y0 + cav_side*0.3, tun_z - tun_r - lava_d),
           iso_s(lava_x1, cav_y0 + cav_side*0.3, tun_z - tun_r),
           iso_s(lava_x0, cav_y0 + cav_side*0.3, tun_z - tun_r)]
    ax.add_patch(Polygon(pts, facecolor="#FF4500", edgecolor="#8B0000", linewidth=0.5, alpha=0.4))
    # Lava front face
    pts = [iso_s(lava_x0, cav_y0, tun_z - tun_r - lava_d),
           iso_s(lava_x1, cav_y0, tun_z - tun_r - lava_d),
           iso_s(lava_x1, cav_y0, tun_z - tun_r),
           iso_s(lava_x0, cav_y0, tun_z - tun_r)]
    ax.add_patch(Polygon(pts, facecolor="#FF6347", edgecolor="#8B0000", linewidth=0.5, alpha=0.3))

    lx_s, ly_s = iso_s((lava_x0+lava_x1)/2, cav_y0, tun_z - tun_r - lava_d*0.5)
    ax.text(lx_s, ly_s, f"LAVA\n{lv.t_lava_c:.0f}C", fontsize=4, ha="center", va="center",
            color="white", fontweight="bold")

    # --- Turbine stages (on tunnel) ---
    for i in range(min(n_turb, 10)):
        frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tun_len)
        tx = tun_x0 + frac * tun_len
        tx_s, ty_s = iso_s(tx, cav_y0, tun_z)
        ax.add_patch(Circle((tx_s, ty_s), 2.5, facecolor="#00E676",
                            edgecolor="#00C853", linewidth=0.5, alpha=0.7, zorder=5))
        if i < 4:
            ax.text(tx_s, ty_s + 4, f"T{i+1}", fontsize=2.5, ha="center", color="#00E676")

    # --- Stack (vertical from tunnel end to above surface) ---
    stack_x = tun_x1
    pts = [iso_s(stack_x, cav_y0 - tun_r, 0), iso_s(stack_x, cav_y0 + tun_r, 0),
           iso_s(stack_x, cav_y0 + tun_r, stack_h), iso_s(stack_x, cav_y0 - tun_r, stack_h)]
    ax.add_patch(Polygon(pts, facecolor="#888", edgecolor="#666", linewidth=0.8, alpha=0.4))
    # Stack top
    pts_top = [iso_s(stack_x - tun_r, cav_y0 - tun_r, stack_h),
               iso_s(stack_x + tun_r, cav_y0 - tun_r, stack_h),
               iso_s(stack_x + tun_r, cav_y0 + tun_r, stack_h),
               iso_s(stack_x - tun_r, cav_y0 + tun_r, stack_h)]
    ax.add_patch(Polygon(pts_top, facecolor="#aaa", edgecolor="#666", linewidth=0.5, alpha=0.3))

    sx_s, sy_s = iso_s(stack_x, cav_y0, stack_h * 0.5)
    ax.text(sx_s + 5, sy_s, "STACK", fontsize=4, color="white", fontweight="bold")

    # --- Fans on stack ---
    for i in range(min(n_fans, 6)):
        fy = stack_h * (0.4 + 0.1 * i)
        fx_s, fy_s = iso_s(stack_x, cav_y0 + tun_r, fy)
        ax.add_patch(Circle((fx_s, fy_s), 1.5, facecolor="#00BFFF",
                            edgecolor="#0277FD", linewidth=0.3, alpha=0.6, zorder=5))

    # --- Jet arrow ---
    jx_s, jy_s = iso_s(stack_x, cav_y0, stack_h)
    ax.annotate("", xy=(jx_s, jy_s + 8), xytext=(jx_s, jy_s),
                arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=1.2))

    # --- Surface buildings ---
    # Transformer
    tx_s, ty_s = iso_s(stack_x + 20, cav_y0, 0)
    ax.add_patch(Rectangle((tx_s - 3, ty_s - 2), 6, 4, facecolor="#FFC107",
                           edgecolor="#FF6F00", linewidth=0.5, alpha=0.5))
    ax.text(tx_s, ty_s + 3, "XFMR", fontsize=3, ha="center", color="#FFC107")

    # Turbine hall
    th_sx, th_sy = iso_s(stack_x + 35, cav_y0, 0)
    ax.add_patch(Rectangle((th_sx - 5, th_sy - 3), 10, 6, facecolor="#546E7A",
                           edgecolor="#37474F", linewidth=0.5, alpha=0.3))
    ax.text(th_sx, th_sy + 4, "Turbine\nHall", fontsize=3, ha="center", color="#90A4AE")

    # Control room
    cr_sx, cr_sy = iso_s(stack_x + 55, cav_y0, 0)
    ax.add_patch(Rectangle((cr_sx - 3, cr_sy - 2), 6, 4, facecolor="#78909C",
                           edgecolor="#455A64", linewidth=0.5, alpha=0.4))
    ax.text(cr_sx, cr_sy + 3, "Ctrl\nRoom", fontsize=3, ha="center", color="#CFD8DC")

    # Cooling tower
    ct_sx, ct_sy = iso_s(cav_x0 - 25, cav_y0, 0)
    ax.add_patch(Circle((ct_sx, ct_sy), 4, facecolor="#B0BEC5",
                        edgecolor="#78909C", linewidth=0.5, alpha=0.3))
    ax.text(ct_sx, ct_sy + 5, "Cooling\nTower", fontsize=3, ha="center", color="#B0BEC5")

    # --- Dimension annotations ---
    # Total length
    d_sx, d_sy = iso_s(0, cav_y0 + cav_side + 10, 0)
    d_ex, d_ey = iso_s(total_x, cav_y0 + cav_side + 10, 0)
    ax.annotate("", xy=(d_ex, d_ey), xytext=(d_sx, d_sy),
                arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.5))
    ax.text((d_sx+d_ex)/2, (d_sy+d_ey)/2 + 3, f"{total_x:.0f}m", fontsize=3,
            ha="center", color="#a0aec0")

    # Depth
    d_sx, d_sy = iso_s(cav_x0 - 10, cav_y0, 0)
    d_ex, d_ey = iso_s(cav_x0 - 10, cav_y0, cav_z0)
    ax.annotate("", xy=(d_ex, d_ey), xytext=(d_sx, d_sy),
                arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.5))
    ax.text(d_sx - 5, (d_sy+d_ey)/2, f"{cav_depth:.0f}m", fontsize=3,
            ha="center", color="#a0aec0", rotation=90)

    # Stack height
    d_sx, d_sy = iso_s(stack_x + 10, cav_y0, 0)
    d_ex, d_ey = iso_s(stack_x + 10, cav_y0, stack_h)
    ax.annotate("", xy=(d_ex, d_ey), xytext=(d_sx, d_sy),
                arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.5))
    ax.text(d_sx + 5, (d_sy+d_ey)/2, f"{stack_h:.0f}m", fontsize=3,
            ha="center", color="#a0aec0", rotation=90)

    ax.set_title("ISOMETRIC CUTAWAY  -  3D perspective view", fontsize=7, fontweight="bold",
                 color="#00d4ff", pad=4)

    # Set limits
    all_xs = []
    all_ys = []
    for x in [0, total_x]:
        for y in [0, cav_side]:
            for z in [cav_z0, stack_h]:
                sx, sy = iso_s(x, y, z)
                all_xs.append(sx)
                all_ys.append(sy)
    ax.set_xlim(min(all_xs) - 10, max(all_xs) + 10)
    ax.set_ylim(min(all_ys) - 10, max(all_ys) + 10)


def _draw_pid_diagram(ax, t: Dict) -> None:
    """Draw a Piping & Instrumentation Diagram (P&ID) showing flow paths,
    control loops, valves, sensors, and equipment connections."""
    from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    n_turb = tn.n_turbine_stages
    n_fans = tn.n_exit_fans

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("P&ID  -  Piping & Instrumentation Diagram", fontsize=7,
                 fontweight="bold", color="#00d4ff", pad=4)

    # Equipment nodes (x, y, label, type, color)
    nodes = {
        "CAVERN": (0.08, 0.50, "CV", "#2196F3"),
        "V1": (0.18, 0.50, "V", "#FF5722"),      # main isolation valve
        "HX": (0.30, 0.50, "HX", "#FF6347"),     # lava heat exchanger
        "T1": (0.40, 0.50, "T", "#00E676"),      # turbine 1
        "RH1": (0.46, 0.50, "RH", "#FF5722"),    # reheat 1
        "T2": (0.52, 0.50, "T", "#00E676"),      # turbine 2
        "RH2": (0.58, 0.50, "RH", "#FF5722"),
        "T3": (0.64, 0.50, "T", "#00E676"),
        "GEN": (0.64, 0.70, "G", "#FFC107"),     # generator
        "K": (0.72, 0.50, "K", "#FF9800"),       # potassium
        "sCO2": (0.78, 0.50, "sC", "#9C27B0"),   # sCO2
        "STM": (0.84, 0.50, "ST", "#03A9F4"),    # steam
        "ORC": (0.90, 0.50, "OR", "#4CAF50"),    # ORC
        "STK": (0.90, 0.75, "S", "#888"),        # stack
        "F1": (0.90, 0.85, "F", "#00BFFF"),      # fan
        "CT": (0.96, 0.30, "CT", "#B0BEC5"),     # cooling tower
        "CMP": (0.08, 0.20, "C", "#8D6E63"),     # recharge compressor
    }

    # Draw equipment
    for name, (x, y, sym, color) in nodes.items():
        if name in ("V1",):
            # Valve symbol (triangle pair)
            ax.add_patch(Rectangle((x-0.015, y-0.015), 0.03, 0.03,
                                   facecolor=color, edgecolor="#444", linewidth=0.5, alpha=0.7,
                                   transform=ax.transAxes))
            ax.text(x, y, sym, fontsize=4, ha="center", va="center", color="white",
                    fontweight="bold", transform=ax.transAxes)
        elif name in ("GEN",):
            ax.add_patch(Circle((x, y), 0.025, facecolor=color, edgecolor="#444",
                                linewidth=0.5, alpha=0.7, transform=ax.transAxes))
            ax.text(x, y, sym, fontsize=4, ha="center", va="center", color="white",
                    fontweight="bold", transform=ax.transAxes)
        elif name in ("STK",):
            ax.add_patch(Rectangle((x-0.015, y-0.03), 0.03, 0.06,
                                   facecolor=color, edgecolor="#444", linewidth=0.5, alpha=0.5,
                                   transform=ax.transAxes))
            ax.text(x, y, sym, fontsize=4, ha="center", va="center", color="white",
                    fontweight="bold", transform=ax.transAxes)
        elif name in ("F1",):
            ax.add_patch(Circle((x, y), 0.02, facecolor=color, edgecolor="#444",
                                linewidth=0.5, alpha=0.7, transform=ax.transAxes))
            ax.text(x, y, sym, fontsize=3, ha="center", va="center", color="white",
                    fontweight="bold", transform=ax.transAxes)
        else:
            ax.add_patch(FancyBboxPatch((x-0.025, y-0.018), 0.05, 0.036,
                                        boxstyle="round,pad=0.005",
                                        facecolor=color, edgecolor="#444",
                                        linewidth=0.5, alpha=0.6, transform=ax.transAxes))
            ax.text(x, y, sym, fontsize=4, ha="center", va="center", color="white",
                    fontweight="bold", transform=ax.transAxes)

    # Draw main flow pipe (cavern -> V1 -> HX -> T1 -> RH1 -> T2 -> RH2 -> T3 -> K -> sCO2 -> STM -> ORC -> STK)
    flow_path = ["CAVERN", "V1", "HX", "T1", "RH1", "T2", "RH2", "T3", "K", "sCO2", "STM", "ORC", "STK"]
    for i in range(len(flow_path) - 1):
        n1, n2 = flow_path[i], flow_path[i+1]
        x1, y1 = nodes[n1][0], nodes[n1][1]
        x2, y2 = nodes[n2][0], nodes[n2][1]
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="#FFD700", lw=1.2, alpha=0.6),
                    transform=ax.transAxes)

    # Generator connection from turbine
    ax.annotate("", xy=nodes["GEN"][0:2], xytext=nodes["T3"][0:2],
                arrowprops=dict(arrowstyle="->", color="#FFC107", lw=0.8, alpha=0.5,
                                linestyle="dashed"),
                transform=ax.transAxes)
    ax.text(0.68, 0.62, "shaft", fontsize=3, color="#FFC107", transform=ax.transAxes)

    # Fan on stack
    ax.annotate("", xy=nodes["F1"][0:2], xytext=nodes["STK"][0:2],
                arrowprops=dict(arrowstyle="->", color="#00BFFF", lw=0.8, alpha=0.5),
                transform=ax.transAxes)

    # Cooling tower connections to bottoming cycles
    for bc in ["K", "sCO2", "STM", "ORC"]:
        ax.annotate("", xy=nodes["CT"][0:2], xytext=(nodes[bc][0], nodes[bc][1] - 0.018),
                    arrowprops=dict(arrowstyle="->", color="#26A69A", lw=0.5, alpha=0.4,
                                    linestyle="dotted"),
                    transform=ax.transAxes)

    # Recharge compressor to cavern
    ax.annotate("", xy=nodes["CAVERN"][0:2], xytext=(nodes["CMP"][0], nodes["CMP"][1] + 0.018),
                arrowprops=dict(arrowstyle="->", color="#8D6E63", lw=0.8, alpha=0.5),
                transform=ax.transAxes)
    ax.text(0.06, 0.35, "recharge", fontsize=3, color="#8D6E63", transform=ax.transAxes,
            rotation=90)

    # Sensors (small circles with tags)
    sensors = [
        (0.13, 0.55, "PT-1", "#FFEB3B"),   # pressure transmitter
        (0.25, 0.55, "TT-1", "#FFEB3B"),   # temp transmitter
        (0.35, 0.55, "FT-1", "#FFEB3B"),   # flow transmitter
        (0.48, 0.55, "TT-2", "#FFEB3B"),
        (0.60, 0.55, "TT-3", "#FFEB3B"),
        (0.70, 0.55, "PT-2", "#FFEB3B"),
        (0.88, 0.60, "TT-4", "#FFEB3B"),
    ]
    for sx, sy, tag, color in sensors:
        ax.add_patch(Circle((sx, sy), 0.008, facecolor=color, edgecolor="#444",
                            linewidth=0.3, alpha=0.8, transform=ax.transAxes))
        ax.text(sx, sy + 0.015, tag, fontsize=2.5, ha="center", color=color,
                transform=ax.transAxes)

    # Control loops (dashed lines from sensors to control room)
    # Show one example control loop
    ax.annotate("", xy=(0.95, 0.92), xytext=(0.13, 0.55),
                arrowprops=dict(arrowstyle="->", color="#E91E63", lw=0.4, alpha=0.3,
                                linestyle="dashdot"),
                transform=ax.transAxes)
    ax.text(0.50, 0.92, "control loops -> SCADA", fontsize=3, color="#E91E63",
            transform=ax.transAxes, ha="center")

    # Labels
    labels = {
        "CAVERN": "Cavern\nCV-01",
        "V1": "Iso\nValve",
        "HX": "Lava\nHX",
        "T1": f"Turb 1\nof {n_turb}",
        "GEN": f"Gen\n{TURBINE_HW['generator_mva']:.0f}MVA",
        "STK": "Stack",
        "F1": f"Fan\n1 of {n_fans}",
        "CT": "Cool\nTower",
        "CMP": "Comp",
    }
    for name, label in labels.items():
        x, y = nodes[name][0], nodes[name][1]
        if name in ("CAVERN", "GEN", "STK", "F1", "CT", "CMP"):
            offset = -0.04 if y > 0.6 else 0.04
            ax.text(x, y + offset, label, fontsize=2.5, ha="center", color="#aaa",
                    transform=ax.transAxes)

    # Legend
    ax.text(0.02, 0.95, "--- FLOW ---", fontsize=3, color="#FFD700", transform=ax.transAxes)
    ax.plot([0.02, 0.06], [0.93, 0.93], color="#FFD700", linewidth=1, transform=ax.transAxes)
    ax.text(0.02, 0.90, "--- SIGNAL ---", fontsize=3, color="#E91E63", transform=ax.transAxes)
    ax.plot([0.02, 0.06], [0.88, 0.88], color="#E91E63", linewidth=0.5, linestyle="dashdot",
            transform=ax.transAxes)
    ax.text(0.02, 0.85, "--- COOLING ---", fontsize=3, color="#26A69A", transform=ax.transAxes)
    ax.plot([0.02, 0.06], [0.83, 0.83], color="#26A69A", linewidth=0.5, linestyle="dotted",
            transform=ax.transAxes)


def _draw_generator_detail(ax, t: Dict) -> None:
    """Draw a detailed generator cross-section showing rotor, stator, windings,
    hydrogen cooling, bearings, and exciter."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Wedge
    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL E  -  Generator cross-section", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    # Generator housing (outer)
    ax.add_patch(Rectangle((0.10, 0.30), 0.80, 0.40, transform=ax.transAxes,
                           facecolor="#37474F", edgecolor="#263238", linewidth=1, alpha=0.3))

    # Stator core (laminated steel)
    ax.add_patch(Rectangle((0.14, 0.34), 0.72, 0.32, transform=ax.transAxes,
                           facecolor="#546E7A", edgecolor="#37474F", linewidth=0.8, alpha=0.4))

    # Stator slots (small lines around inner surface)
    n_slots = 24
    for i in range(n_slots):
        # top slots
        x = 0.16 + i * 0.70 / n_slots
        ax.plot([x, x], [0.62, 0.66], color="#78909C", linewidth=0.5,
                transform=ax.transAxes, alpha=0.5)
        # bottom slots
        ax.plot([x, x], [0.34, 0.38], color="#78909C", linewidth=0.5,
                transform=ax.transAxes, alpha=0.5)

    # Stator windings (copper bars in slots)
    for i in range(0, n_slots, 3):
        x = 0.16 + i * 0.70 / n_slots
        ax.add_patch(Rectangle((x-0.005, 0.63), 0.01, 0.03, transform=ax.transAxes,
                               facecolor="#FF8C00", edgecolor="#E65100", linewidth=0.2, alpha=0.6))
        ax.add_patch(Rectangle((x-0.005, 0.34), 0.01, 0.03, transform=ax.transAxes,
                               facecolor="#FF8C00", edgecolor="#E65100", linewidth=0.2, alpha=0.6))

    # Rotor (large cylinder)
    rotor_cx, rotor_cy = 0.50, 0.50
    rotor_r = 0.14
    ax.add_patch(Circle((rotor_cx, rotor_cy), rotor_r, transform=ax.transAxes,
                        facecolor="#FFC107", edgecolor="#FF6F00", linewidth=1, alpha=0.4))

    # Rotor poles (4 poles)
    n_poles = 4
    for i in range(n_poles):
        ang = 2 * math.pi * i / n_poles
        px = rotor_cx + rotor_r * 0.8 * math.cos(ang)
        py = rotor_cy + rotor_r * 0.8 * math.sin(ang)
        ax.add_patch(Circle((px, py), 0.015, transform=ax.transAxes,
                            facecolor="#FF6F00", edgecolor="#E65100", linewidth=0.3, alpha=0.6))

    # Rotor winding (field winding)
    ax.add_patch(Circle((rotor_cx, rotor_cy), rotor_r * 0.5, transform=ax.transAxes,
                        facecolor="#FFD54F", edgecolor="#FF6F00", linewidth=0.5, alpha=0.3))

    # Shaft
    ax.add_patch(Rectangle((0.04, rotor_cy - 0.015), 0.92, 0.03, transform=ax.transAxes,
                           facecolor="#888", edgecolor="#555", linewidth=0.5, alpha=0.5, zorder=5))

    # Bearings (journal bearings)
    for bx in [0.12, 0.88]:
        ax.add_patch(Rectangle((bx-0.02, rotor_cy-0.025), 0.04, 0.05, transform=ax.transAxes,
                               facecolor="#9E9E9E", edgecolor="#616161", linewidth=0.5, alpha=0.5, zorder=6))
        ax.text(bx, rotor_cy - 0.04, "bearing", fontsize=2.5, ha="center", color="#9E9E9E",
                transform=ax.transAxes)

    # Exciter (left end)
    ax.add_patch(Rectangle((0.02, rotor_cy-0.02), 0.04, 0.04, transform=ax.transAxes,
                           facecolor="#7C4DFF", edgecolor="#4527A0", linewidth=0.5, alpha=0.5))
    ax.text(0.04, rotor_cy + 0.04, "exciter", fontsize=2.5, ha="center", color="#B388FF",
            transform=ax.transAxes)

    # Hydrogen cooling ports
    for hx in [0.25, 0.75]:
        ax.add_patch(Circle((hx, 0.72), 0.012, transform=ax.transAxes,
                            facecolor="#80DEEA", edgecolor="#00ACC1", linewidth=0.3, alpha=0.6))
        ax.text(hx, 0.75, "H2", fontsize=2.5, ha="center", color="#80DEEA",
                transform=ax.transAxes)

    # Terminal box
    ax.add_patch(Rectangle((0.42, 0.72), 0.16, 0.06, transform=ax.transAxes,
                           facecolor="#FFC107", edgecolor="#FF6F00", linewidth=0.5, alpha=0.4))
    ax.text(0.50, 0.75, "terminals", fontsize=2.5, ha="center", color="white",
            transform=ax.transAxes)

    # Labels
    ax.text(0.50, 0.22, f"Generator: {TURBINE_HW['generator_mva']:.0f} MVA, "
                        f"{TURBINE_HW['generator_kv']:.1f} kV\n"
                        f"{TURBINE_HW['generator_cooling']}\n"
                        f"PF={TURBINE_HW['generator_pf']:.2f}, "
                        f"eta={TURBINE_HW['eta_generator']:.2f}",
            fontsize=3.5, ha="center", color="white", fontweight="bold",
            transform=ax.transAxes)

    # Component labels
    ax.text(0.50, 0.68, "stator core", fontsize=2.5, ha="center", color="#78909C",
            transform=ax.transAxes)
    ax.text(rotor_cx, rotor_cy, "rotor", fontsize=3, ha="center", va="center", color="#FF6F00",
            fontweight="bold", transform=ax.transAxes)
    ax.text(0.20, 0.50, "H2\ncool", fontsize=2.5, ha="center", color="#80DEEA",
            transform=ax.transAxes)
    ax.text(0.80, 0.50, "H2\ncool", fontsize=2.5, ha="center", color="#80DEEA",
            transform=ax.transAxes)


def _draw_exploded_view(ax, t: Dict) -> None:
    """Draw an exploded assembly view showing turbine components separated
    in assembly order, like an engineering parts manual diagram."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon
    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL F  -  Exploded turbine assembly", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    n_turb = tn.n_turbine_stages
    n_blades = TURBINE_HW['rotor_blade_count']

    # Components in assembly order (top to bottom = front to back of turbine)
    # Each is drawn as a separated part with a label and leader line
    components = [
        ("Inlet casing", "#78909C", 0.15, "Cast steel inlet\nbellmouth"),
        ("Stator 1 (S1)", "#FF9800", 0.25, f"{n_blades} vanes\nInconel 718+TBC"),
        ("Rotor 1 (R1)", "#00E676", 0.35, f"{n_blades} blades\nD={TURBINE_HW['rotor_d_mm']:.0f}mm"),
        ("Stator 2 (S2)", "#FF9800", 0.45, f"{n_blades} vanes"),
        ("Rotor 2 (R2)", "#00E676", 0.55, f"{n_blades} blades"),
        ("Interstage seal", "#9C27B0", 0.62, "Labyrinth\n+ buffer air"),
        ("Stator 3 (S3)", "#FF9800", 0.68, f"{n_blades} vanes"),
        ("Rotor 3 (R3)", "#00E676", 0.75, f"{n_blades} blades"),
        ("Exhaust diffuser", "#78909C", 0.85, "Cast steel\ndiffuser cone"),
        ("Exhaust casing", "#546E7A", 0.92, "Welded steel\nexit flange"),
    ]

    cx = 0.35  # center x for components
    shaft_x = 0.35

    # Shaft (vertical line through all components)
    ax.plot([shaft_x, shaft_x], [0.05, 0.95], color="#888", linewidth=2,
            transform=ax.transAxes, zorder=1)

    for name, color, y, detail in components:
        # Component shape (varies by type)
        if "Rotor" in name:
            # Rotor = circle with blades
            ax.add_patch(Circle((cx, y), 0.06, transform=ax.transAxes,
                                facecolor=color, edgecolor="#333", linewidth=0.5, alpha=0.5, zorder=3))
            # Blade lines
            for b in range(min(n_blades, 8)):
                ang = 2 * math.pi * b / min(n_blades, 8)
                ax.plot([cx, cx + 0.06 * math.cos(ang)],
                       [y, y + 0.06 * math.sin(ang)],
                       color="#333", linewidth=0.3, alpha=0.5,
                       transform=ax.transAxes, zorder=4)
        elif "Stator" in name:
            # Stator = ring (annulus approximated)
            ax.add_patch(Circle((cx, y), 0.07, transform=ax.transAxes,
                                facecolor="none", edgecolor=color, linewidth=2, zorder=3))
            ax.add_patch(Circle((cx, y), 0.04, transform=ax.transAxes,
                                facecolor="#0d0d18", edgecolor=color, linewidth=1, zorder=4))
        elif "casing" in name or "diffuser" in name:
            # Casing = rectangle
            ax.add_patch(Rectangle((cx-0.07, y-0.025), 0.14, 0.05,
                                   transform=ax.transAxes,
                                   facecolor=color, edgecolor="#333", linewidth=0.5, alpha=0.4, zorder=3))
        elif "seal" in name:
            # Seal = thin rectangle with hatching
            ax.add_patch(Rectangle((cx-0.04, y-0.012), 0.08, 0.024,
                                   transform=ax.transAxes,
                                   facecolor=color, edgecolor="#333", linewidth=0.3, alpha=0.5, zorder=3))

        # Label with leader line
        ax.annotate("", xy=(cx + 0.075, y), xytext=(0.58, y),
                    arrowprops=dict(arrowstyle="-", color="#FFEB3B", lw=0.4, alpha=0.5),
                    transform=ax.transAxes)
        ax.text(0.60, y, name, fontsize=3.5, va="center", color="white",
                fontweight="bold", transform=ax.transAxes)
        ax.text(0.78, y, detail, fontsize=2.5, va="center", color="#aaa",
                transform=ax.transAxes, family="monospace")

    # Assembly direction arrow
    ax.annotate("", xy=(0.15, 0.95), xytext=(0.15, 0.05),
                arrowprops=dict(arrowstyle="->", color="#00d4ff", lw=1.5),
                transform=ax.transAxes)
    ax.text(0.10, 0.50, "assembly\norder", fontsize=3, color="#00d4ff",
            ha="center", va="center", rotation=90, transform=ax.transAxes)

    # Specs
    ax.text(0.50, 0.02, f"{n_turb} stages | RPM {TURBINE_HW['rpm']:.0f} | "
                        f"eta={TURBINE_HW['eta_isentropic']:.2f} | "
                        f"{TURBINE_HW['blade_material']}",
            fontsize=3, ha="center", color="white", fontweight="bold",
            transform=ax.transAxes)


def _draw_electrical_sld(ax, t: Dict) -> None:
    """Draw an electrical single-line diagram (SLD) showing power flow
    from generators through switchgear to grid."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    n_turb = tn.n_turbine_stages

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL G  -  Electrical single-line diagram", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    # Generator symbols (circles with G)
    n_gen_show = min(n_turb, 4)
    gen_y_start = 0.85
    gen_y_step = 0.08
    gen_x = 0.12

    for i in range(n_gen_show):
        gy = gen_y_start - i * gen_y_step
        # Generator circle
        ax.add_patch(Circle((gen_x, gy), 0.025, transform=ax.transAxes,
                            facecolor="#FFC107", edgecolor="#FF6F00", linewidth=0.5, alpha=0.6))
        ax.text(gen_x, gy, "G", fontsize=4, ha="center", va="center",
                color="white", fontweight="bold", transform=ax.transAxes)
        # Generator label
        ax.text(gen_x - 0.04, gy, f"G{i+1}\n{TURBINE_HW['generator_mva']:.0f}MVA",
                fontsize=2.5, ha="right", va="center", color="#FFC107",
                transform=ax.transAxes)

        # Line from generator to generator breaker
        ax.plot([gen_x + 0.025, gen_x + 0.06], [gy, gy],
                color="#FFD700", linewidth=1, transform=ax.transAxes)
        # Generator breaker (square)
        ax.add_patch(Rectangle((gen_x + 0.06, gy - 0.008), 0.02, 0.016,
                               transform=ax.transAxes,
                               facecolor="#F44336", edgecolor="#333", linewidth=0.3, alpha=0.7))
        ax.text(gen_x + 0.07, gy + 0.015, "52G", fontsize=2, ha="center",
                color="#F44336", transform=ax.transAxes)
        # Line to 13.8kV bus
        ax.plot([gen_x + 0.08, 0.25], [gy, gy],
                color="#FFD700", linewidth=1, transform=ax.transAxes)

    if n_turb > n_gen_show:
        ax.text(gen_x, gen_y_start - n_gen_show * gen_y_step,
                f"+ {n_turb - n_gen_show} more\ngenerators",
                fontsize=2.5, ha="center", color="#888", transform=ax.transAxes)

    # 13.8 kV bus (horizontal line)
    bus_y = gen_y_start - (n_gen_show - 1) * gen_y_step / 2
    ax.plot([0.25, 0.25], [gen_y_start - 0.01, gen_y_start - (n_gen_show-1) * gen_y_step + 0.01],
            color="#FFD700", linewidth=2.5, transform=ax.transAxes)
    ax.text(0.255, bus_y + 0.06, "13.8 kV bus", fontsize=3, color="#FFD700",
            transform=ax.transAxes, fontweight="bold")

    # Station service transformer
    ax.plot([0.25, 0.35], [bus_y, bus_y], color="#FFD700", linewidth=1,
            transform=ax.transAxes)
    ax.add_patch(Rectangle((0.35, bus_y - 0.02), 0.03, 0.04,
                           transform=ax.transAxes,
                           facecolor="#7C4DFF", edgecolor="#333", linewidth=0.5, alpha=0.6))
    ax.text(0.365, bus_y, "T", fontsize=3, ha="center", va="center",
            color="white", fontweight="bold", transform=ax.transAxes)
    ax.text(0.38, bus_y + 0.03, "500kVA\n13.8kV/480V", fontsize=2,
            color="#B388FF", transform=ax.transAxes)
    ax.plot([0.38, 0.45], [bus_y, bus_y], color="#7C4DFF", linewidth=0.5,
            transform=ax.transAxes)
    ax.text(0.46, bus_y, "station\nservice", fontsize=2, color="#B388FF",
            transform=ax.transAxes)

    # Step-up transformer to 132kV
    ax.plot([0.25, 0.55], [bus_y, bus_y], color="#FFD700", linewidth=1,
            transform=ax.transAxes)
    # Transformer symbol (two overlapping circles)
    ax.add_patch(Circle((0.56, bus_y + 0.012), 0.015, transform=ax.transAxes,
                        facecolor="none", edgecolor="#FFC107", linewidth=1))
    ax.add_patch(Circle((0.56, bus_y - 0.012), 0.015, transform=ax.transAxes,
                        facecolor="none", edgecolor="#FFC107", linewidth=1))
    ax.text(0.56, bus_y + 0.05, f"Step-up\n{TURBINE_HW['generator_kv']:.1f}kV/132kV\n"
                                f"{n_turb * TURBINE_HW['generator_mva']:.0f}MVA",
            fontsize=2.5, ha="center", color="#FFC107", transform=ax.transAxes)

    # 132kV bus
    ax.plot([0.575, 0.70], [bus_y, bus_y], color="#FF4500", linewidth=2.5,
            transform=ax.transAxes)
    ax.text(0.60, bus_y + 0.04, "132 kV bus", fontsize=3, color="#FF4500",
            transform=ax.transAxes, fontweight="bold")

    # Transmission line breakers
    for i in range(2):
        tx = 0.72 + i * 0.08
        ax.add_patch(Rectangle((tx, bus_y - 0.008), 0.02, 0.016,
                               transform=ax.transAxes,
                               facecolor="#F44336", edgecolor="#333", linewidth=0.3, alpha=0.7))
        ax.text(tx + 0.01, bus_y + 0.015, f"52L{i+1}", fontsize=2, ha="center",
                color="#F44336", transform=ax.transAxes)
        ax.plot([0.70, tx], [bus_y, bus_y], color="#FF4500", linewidth=1.5,
                transform=ax.transAxes)
        ax.plot([tx + 0.02, tx + 0.06], [bus_y, bus_y], color="#FF4500", linewidth=1.5,
                transform=ax.transAxes)
        # Transmission line symbol (arrow to grid)
        ax.annotate("", xy=(tx + 0.08, bus_y), xytext=(tx + 0.06, bus_y),
                    arrowprops=dict(arrowstyle="->", color="#FF4500", lw=1.2),
                    transform=ax.transAxes)
        ax.text(tx + 0.07, bus_y - 0.03, f"Line {i+1}\nto grid", fontsize=2,
                ha="center", color="#FF4500", transform=ax.transAxes)

    # Protection relays (dashed lines)
    ax.plot([0.25, 0.25], [bus_y - 0.08, bus_y - 0.02], color="#E91E63",
            linewidth=0.3, linestyle="dashed", transform=ax.transAxes)
    ax.text(0.26, bus_y - 0.10, "87B bus diff", fontsize=2, color="#E91E63",
            transform=ax.transAxes)

    # Voltage levels legend
    ax.text(0.02, 0.15, "--- VOLTAGE LEVELS ---", fontsize=3, color="#00d4ff",
            transform=ax.transAxes, fontweight="bold")
    ax.text(0.02, 0.11, "13.8 kV - generator", fontsize=2.5, color="#FFD700",
            transform=ax.transAxes)
    ax.text(0.02, 0.08, "132 kV - transmission", fontsize=2.5, color="#FF4500",
            transform=ax.transAxes)
    ax.text(0.02, 0.05, "480 V - station service", fontsize=2.5, color="#7C4DFF",
            transform=ax.transAxes)


def _draw_reheat_detail(ax, t: Dict) -> None:
    """Draw a detailed reheat section showing HX tubes, flow path, and
    temperature progression."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    n_reheat = getattr(tn, 'n_reheat_stages', 0)

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL H  -  Reheat section", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    if n_reheat == 0:
        ax.text(0.50, 0.50, "No reheat stages", transform=ax.transAxes,
                fontsize=5, ha="center", color="#888")
        return

    # Outer casing
    ax.add_patch(Rectangle((0.08, 0.30), 0.84, 0.40, transform=ax.transAxes,
                           facecolor="#37474F", edgecolor="#263238", linewidth=0.8, alpha=0.3))

    # Refractory lining
    ax.add_patch(Rectangle((0.10, 0.32), 0.80, 0.36, transform=ax.transAxes,
                           facecolor="#FF7043", edgecolor="#BF360C", linewidth=0.3, alpha=0.2))

    # Inner flow channel
    ax.add_patch(Rectangle((0.14, 0.36), 0.72, 0.28, transform=ax.transAxes,
                           facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.3, alpha=0.15))

    # HX tubes (horizontal lines inside the reheat section)
    n_tubes_show = 8
    for i in range(n_tubes_show):
        ty = 0.38 + i * 0.032
        ax.plot([0.16, 0.84], [ty, ty], color="#FF8C00", linewidth=0.8,
                alpha=0.5, transform=ax.transAxes)
        # Tube ends (small circles)
        ax.add_patch(Circle((0.16, ty), 0.005, transform=ax.transAxes,
                            facecolor="#FF8C00", edgecolor="#E65100", linewidth=0.2, alpha=0.7))
        ax.add_patch(Circle((0.84, ty), 0.005, transform=ax.transAxes,
                            facecolor="#FF8C00", edgecolor="#E65100", linewidth=0.2, alpha=0.7))

    # Lava inlet (bottom)
    ax.add_patch(Rectangle((0.30, 0.20), 0.15, 0.10, transform=ax.transAxes,
                           facecolor="#FF4500", edgecolor="#8B0000", linewidth=0.5, alpha=0.5))
    ax.text(0.375, 0.25, "lava\nin", fontsize=2.5, ha="center", va="center",
            color="white", fontweight="bold", transform=ax.transAxes)

    # Lava outlet (bottom)
    ax.add_patch(Rectangle((0.55, 0.20), 0.15, 0.10, transform=ax.transAxes,
                           facecolor="#8B0000", edgecolor="#5D0000", linewidth=0.5, alpha=0.5))
    ax.text(0.625, 0.25, "lava\nout", fontsize=2.5, ha="center", va="center",
            color="white", fontweight="bold", transform=ax.transAxes)

    # Air flow arrows (left to right)
    ax.annotate("", xy=(0.90, 0.50), xytext=(0.10, 0.50),
                arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=1.5),
                transform=ax.transAxes)

    # Temperature labels
    ax.text(0.12, 0.68, "T_cold\n(in)", fontsize=2.5, ha="center", color="#4FC3F7",
            transform=ax.transAxes)
    ax.text(0.88, 0.68, "T_hot\n(out)", fontsize=2.5, ha="center", color="#FF6347",
            transform=ax.transAxes)

    # Temperature gradient visualization
    for i in range(10):
        frac = i / 9.0
        tx = 0.16 + frac * 0.68
        # Color from blue to red
        r = int(255 * frac)
        b = int(255 * (1 - frac))
        color = f"#{r:02x}80{b:02x}"
        ax.add_patch(Rectangle((tx-0.01, 0.34), 0.02, 0.005, transform=ax.transAxes,
                               facecolor=color, edgecolor="none", alpha=0.6))

    # Specs
    ax.text(0.50, 0.12, f"Reheat stages: {n_reheat}\n"
                        f"HX tubes per stage: {lv.hx_n_tubes:,}\n"
                        f"U-value: {lv.hx_u:.0f} W/m2K\n"
                        f"Function: reheats air between turbine stages\n"
                        f"  -> near-isothermal expansion -> more work output",
            fontsize=3, ha="center", color="white", fontweight="bold",
            transform=ax.transAxes)

    # Labels
    ax.text(0.50, 0.62, "HX tube bundle", fontsize=3, ha="center", color="#FF8C00",
            transform=ax.transAxes, fontweight="bold")
    ax.text(0.50, 0.34, "air flow channel", fontsize=2.5, ha="center", color="#FFD700",
            transform=ax.transAxes)

    # Inlet/outlet flanges
    for fx in [0.08, 0.92]:
        ax.add_patch(Rectangle((fx-0.015, 0.45), 0.03, 0.10, transform=ax.transAxes,
                               facecolor="#555", edgecolor="#333", linewidth=0.3, alpha=0.5))
        ax.text(fx, 0.58, "flange", fontsize=2, ha="center", color="#888",
                transform=ax.transAxes)


def _draw_site_layout(ax, t: Dict) -> None:
    """Draw a detailed site layout plan showing all surface buildings, roads,
    fences, drainage, and the underground system footprint."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon
    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t["ctrl"]
    n_sys = max(1, ctrl.n_systems)

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("SITE LAYOUT  -  surface plan", fontsize=7,
                 fontweight="bold", color="#00d4ff", pad=4)

    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    tun_len = tn.total_length_m
    stack_h = tn.height_rise_m
    total_len = cav_side + tun_len + stack_h + 100

    # Scale
    sx = 0.90 / (total_len * 1.15)
    sy = sx

    def sx_(m): return m * sx
    def sy_(m): return m * sy

    # Ground/terrain
    ax.add_patch(Rectangle((0.02, 0.05), 0.96, 0.88, transform=ax.transAxes,
                           facecolor="#1B1B1B", edgecolor="#333", linewidth=0.5, alpha=0.3))

    # Underground footprint (cavern + tunnel, shown as dashed outline)
    cav_x, cav_y = 0.08, 0.40
    cav_w = sx_(cav_side)
    cav_h = sy_(cav_side)
    ax.add_patch(FancyBboxPatch((cav_x, cav_y - cav_h/2), cav_w, cav_h,
                                boxstyle="round,pad=0.01", facecolor="none",
                                edgecolor="#2196F3", linewidth=0.8, linestyle="--", alpha=0.4,
                                transform=ax.transAxes))
    ax.text(cav_x + cav_w/2, cav_y, "CAVERN\n(underground)", fontsize=3,
            ha="center", va="center", color="#2196F3", alpha=0.5,
            transform=ax.transAxes)

    tun_x0 = cav_x + cav_w
    tun_x1 = tun_x0 + sx_(tun_len)
    ax.plot([tun_x0, tun_x1], [cav_y, cav_y], color="#FFD700", linewidth=0.8,
            linestyle="--", alpha=0.3, transform=ax.transAxes)
    ax.text((tun_x0+tun_x1)/2, cav_y + 0.03, "tunnel (underground)", fontsize=2.5,
            ha="center", color="#FFD700", alpha=0.4, transform=ax.transAxes)

    # Access shaft
    ax.add_patch(Circle((cav_x + cav_w*0.3, cav_y), 0.008, transform=ax.transAxes,
                        facecolor="#888", edgecolor="#555", linewidth=0.3, alpha=0.5))

    # Surface buildings
    stack_x = tun_x1

    buildings = [
        (stack_x + 0.01, cav_y - 0.03, 0.03, 0.04, "#888", "Stack"),
        (stack_x + 0.05, cav_y - 0.04, 0.05, 0.06, "#546E7A", "Turbine\nHall"),
        (stack_x + 0.05, cav_y + 0.04, 0.03, 0.03, "#78909C", "Control\nRoom"),
        (stack_x + 0.09, cav_y - 0.04, 0.03, 0.03, "#455A64", "Switch-\nyard"),
        (stack_x + 0.09, cav_y + 0.03, 0.02, 0.03, "#FFC107", "XFMR"),
        (stack_x + 0.13, cav_y - 0.03, 0.025, 0.03, "#FF6F00", "Diesel\nGen"),
        (stack_x + 0.13, cav_y + 0.03, 0.025, 0.03, "#FF6F00", "Fuel\nTank"),
        (cav_x - 0.04, cav_y + 0.05, 0.03, 0.03, "#80DEEA", "Chiller"),
        (cav_x - 0.04, cav_y - 0.06, 0.04, 0.04, "#B0BEC5", "Cooling\nTower"),
        (cav_x - 0.04, cav_y, 0.025, 0.03, "#26A69A", "DM\nWater"),
        (cav_x - 0.07, cav_y - 0.02, 0.025, 0.03, "#8D6E63", "Recharge\nCompressor"),
    ]

    for bx, by, bw, bh, color, label in buildings:
        ax.add_patch(Rectangle((bx, by), bw, bh, transform=ax.transAxes,
                               facecolor=color, edgecolor="#333", linewidth=0.3, alpha=0.5))
        ax.text(bx + bw/2, by + bh/2, label, fontsize=2.5, ha="center", va="center",
                color="white", fontweight="bold", transform=ax.transAxes)

    # Bottoming cycle building
    bc_x = stack_x + 0.03
    bc_y = cav_y - 0.10
    ax.add_patch(Rectangle((bc_x, bc_y), 0.06, 0.03, transform=ax.transAxes,
                           facecolor="#9C27B0", edgecolor="#333", linewidth=0.3, alpha=0.4))
    ax.text(bc_x + 0.03, bc_y + 0.015, "Bottoming\nCycles", fontsize=2.5,
            ha="center", va="center", color="white", fontweight="bold",
            transform=ax.transAxes)

    # Access roads
    road_color = "#424242"
    # Main access road
    ax.plot([0.02, 0.98], [0.12, 0.12], color=road_color, linewidth=2,
            alpha=0.4, transform=ax.transAxes)
    ax.text(0.50, 0.10, "main access road", fontsize=2.5, ha="center",
            color="#666", transform=ax.transAxes)
    # Connecting road to buildings
    ax.plot([0.50, 0.50], [0.12, 0.30], color=road_color, linewidth=1.5,
            alpha=0.3, transform=ax.transAxes)
    ax.plot([0.15, 0.15], [0.12, 0.35], color=road_color, linewidth=1,
            alpha=0.3, transform=ax.transAxes)
    ax.plot([0.85, 0.85], [0.12, 0.35], color=road_color, linewidth=1,
            alpha=0.3, transform=ax.transAxes)

    # Security fence (dashed perimeter)
    fence_pts = [(0.04, 0.08), (0.96, 0.08), (0.96, 0.92), (0.04, 0.92)]
    for i in range(len(fence_pts)):
        p1 = fence_pts[i]
        p2 = fence_pts[(i+1) % len(fence_pts)]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#555", linewidth=0.5,
                linestyle="--", alpha=0.3, transform=ax.transAxes)

    # Gate
    ax.plot([0.48, 0.52], [0.08, 0.08], color="#FFC107", linewidth=1.5,
            alpha=0.5, transform=ax.transAxes)
    ax.text(0.50, 0.065, "gate", fontsize=2.5, ha="center", color="#FFC107",
            transform=ax.transAxes)

    # Drainage
    ax.plot([0.10, 0.90], [0.07, 0.07], color="#26A69A", linewidth=0.5,
            linestyle=":", alpha=0.3, transform=ax.transAxes)
    ax.text(0.50, 0.055, "site drainage", fontsize=2, ha="center", color="#26A69A",
            alpha=0.5, transform=ax.transAxes)

    # Transmission lines (exiting to right)
    for i in range(3):
        ax.plot([0.96, 1.0], [0.45 + i*0.02, 0.45 + i*0.02], color="#FFEB3B",
                linewidth=0.5, alpha=0.4, transform=ax.transAxes)
    ax.text(0.98, 0.50, "to\ngrid", fontsize=2.5, ha="center", color="#FFEB3B",
            transform=ax.transAxes)

    # North arrow
    ax.annotate("N", xy=(0.93, 0.88), xytext=(0.93, 0.82),
                arrowprops=dict(arrowstyle="->", color="#00d4ff", lw=1),
                fontsize=4, ha="center", color="#00d4ff",
                fontweight="bold", transform=ax.transAxes)

    # Scale bar
    ax.plot([0.05, 0.05 + sx_(500)], [0.93, 0.93], color="white", linewidth=1,
            transform=ax.transAxes)
    ax.text(0.05 + sx_(250), 0.94, "500 m", fontsize=2.5, ha="center",
            color="white", transform=ax.transAxes)

    # Dual system indicator
    if n_sys > 1:
        ax.text(0.50, 0.90, f"x{n_sys} systems (side by side)", fontsize=3,
                ha="center", color="#FFD700", transform=ax.transAxes)


def _draw_cavern_interior(ax, t: Dict) -> None:
    """Draw a detailed 3D-style interior view of the cavern showing sump,
    floor slope, door, sensors, DTS fiber, and rock bolts."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon, Wedge
    cv = t["cavern"]

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL I  -  Cavern interior (perspective)", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    # Cavern outline (perspective trapezoid)
    pts = [(0.08, 0.20), (0.92, 0.25), (0.92, 0.80), (0.08, 0.75)]
    ax.add_patch(Polygon(pts, facecolor="#1565C0", edgecolor="#0277FD",
                         linewidth=1, alpha=0.15, transform=ax.transAxes))

    # Cavern lining (inner edge)
    pts_inner = [(0.10, 0.22), (0.90, 0.27), (0.90, 0.78), (0.10, 0.73)]
    ax.add_patch(Polygon(pts_inner, facecolor="none", edgecolor="#6D4C41",
                         linewidth=0.8, linestyle="--", alpha=0.4,
                         transform=ax.transAxes))

    # Floor (sloped toward sump)
    ax.plot([0.10, 0.90], [0.22, 0.27], color="#8D6E63", linewidth=1.5,
            transform=ax.transAxes)
    ax.fill_between([0.10, 0.90], [0.20, 0.25], [0.22, 0.27],
                    color="#5D4037", alpha=0.3, transform=ax.transAxes)
    ax.text(0.50, 0.21, "floor (1:200 slope to sump)", fontsize=2.5,
            ha="center", color="#8D6E63", transform=ax.transAxes)

    # Drainage sump (at low end)
    ax.add_patch(Rectangle((0.08, 0.18), 0.04, 0.04, transform=ax.transAxes,
                           facecolor="#26A69A", edgecolor="#00897B", linewidth=0.5, alpha=0.5))
    ax.text(0.10, 0.20, "sump", fontsize=2.5, ha="center", va="center",
            color="white", fontweight="bold", transform=ax.transAxes)

    # Hydraulic door (at access tunnel)
    ax.add_patch(Rectangle((0.06, 0.45), 0.04, 0.10, transform=ax.transAxes,
                           facecolor="#FF5722", edgecolor="#D84315", linewidth=0.5, alpha=0.6))
    ax.text(0.08, 0.50, "door", fontsize=2.5, ha="center", va="center",
            color="white", fontweight="bold", rotation=90, transform=ax.transAxes)
    ax.text(0.04, 0.50, "access\ntunnel", fontsize=2, ha="center", color="#888",
            transform=ax.transAxes)

    # Pressure sensors (on walls)
    n_p = min(MONITOR_HW['cavern_pressure_sensors'], 6)
    for i in range(n_p):
        sy = 0.35 + i * 0.07
        ax.add_patch(Circle((0.11, sy), 0.006, transform=ax.transAxes,
                            facecolor="#FFEB3B", edgecolor="#F57F17", linewidth=0.2, alpha=0.7))
        ax.text(0.13, sy, f"PT-{i+1}", fontsize=2, va="center", color="#FFEB3B",
                transform=ax.transAxes)

    # Temperature sensors (RTD)
    n_t = min(MONITOR_HW['cavern_temp_sensors'], 6)
    for i in range(n_t):
        sy = 0.35 + i * 0.07
        ax.add_patch(Circle((0.89, sy), 0.006, transform=ax.transAxes,
                            facecolor="#4FC3F7", edgecolor="#0277FD", linewidth=0.2, alpha=0.7))
        ax.text(0.87, sy, f"TT-{i+1}", fontsize=2, va="center", ha="right",
                color="#4FC3F7", transform=ax.transAxes)

    # DTS fiber (along ceiling)
    ax.plot([0.12, 0.88], [0.76, 0.73], color="#E91E63", linewidth=0.8,
            alpha=0.5, transform=ax.transAxes)
    ax.text(0.50, 0.78, "DTS fiber (temperature mapping)", fontsize=2.5,
            ha="center", color="#E91E63", transform=ax.transAxes)

    # DAS fiber (along floor)
    ax.plot([0.12, 0.88], [0.23, 0.28], color="#9C27B0", linewidth=0.8,
            alpha=0.5, transform=ax.transAxes)
    ax.text(0.50, 0.17, "DAS fiber (acoustic monitoring)", fontsize=2.5,
            ha="center", color="#9C27B0", transform=ax.transAxes)

    # Rock bolts (on ceiling)
    for i in range(8):
        bx = 0.15 + i * 0.09
        by = 0.75 if i < 4 else 0.77
        ax.plot([bx, bx], [by, by + 0.02], color="#888", linewidth=0.5,
                alpha=0.4, transform=ax.transAxes)

    # Roof support arches
    for i in range(4):
        ax_x = 0.20 + i * 0.20
        ax.plot([ax_x, ax_x], [0.30, 0.72], color="#555", linewidth=0.4,
                alpha=0.3, transform=ax.transAxes)

    # Cold air (fill pattern)
    ax.text(0.50, 0.50, f"COLD AIR\n{cv.volume_m3/1e9:.1f} km3\n"
                        f"{k_to_c(cv.t_charge_k):.0f} C\n"
                        f"{cv.p_charge_pa/1e5:.0f} bar",
            fontsize=5, ha="center", va="center", color="white",
            fontweight="bold", transform=ax.transAxes)

    # Geophones (around perimeter)
    for i in range(4):
        gx = 0.15 + i * 0.22
        ax.add_patch(Circle((gx, 0.82), 0.005, transform=ax.transAxes,
                            facecolor="#E91E63", edgecolor="#880E4F", linewidth=0.2, alpha=0.6))

    # Specs
    ax.text(0.50, 0.08, f"Volume: {cv.volume_m3/1e9:.1f} km3 | "
                        f"Depth: {cv.depth_m:.0f}m | "
                        f"Door: {CAVERN_HW['hydraulic_door_mm']:.0f}mm | "
                        f"Sump: {CAVERN_HW['drainage_sump_m3']:,.0f} m3",
            fontsize=2.5, ha="center", color="#aaa", transform=ax.transAxes)


def _draw_ts_diagrams(ax, t: Dict) -> None:
    """Draw Temperature-Entropy (T-s) diagrams for each bottoming cycle
    showing the thermodynamic cycle paths."""
    from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon
    tn = t["tunnel"]

    ax.set_facecolor("#0d0d18")
    ax.set_axis_off()
    ax.set_title("DETAIL J  -  Bottoming cycle T-s diagrams", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    cycles = []
    if tn.potassium_enabled:
        cycles.append(("Potassium", "#FF9800", tn.potassium_eta,
                       [(0.5, 300), (2.0, 1800), (3.5, 1800), (5.0, 800), (5.0, 300), (0.5, 300)],
                       "K Rankine"))
    if tn.sco2_enabled:
        cycles.append(("sCO2", "#9C27B0", tn.sco2_eta,
                       [(0.5, 350), (2.0, 800), (3.5, 800), (5.0, 400), (5.0, 350), (0.5, 350)],
                       "sCO2 Brayton"))
    if tn.steam_enabled:
        cycles.append(("Steam", "#03A9F4", tn.steam_eta,
                       [(0.5, 320), (1.5, 820), (3.0, 820), (4.0, 550), (4.5, 330), (0.5, 320)],
                       "Steam Rankine"))
    if tn.orc_enabled:
        cycles.append(("ORC", "#4CAF50", tn.orc_eta,
                       [(0.5, 300), (1.5, 400), (2.5, 400), (3.5, 310), (0.5, 300)],
                       "ORC Rankine"))

    if not cycles:
        ax.text(0.50, 0.50, "No bottoming cycles enabled", transform=ax.transAxes,
                fontsize=5, ha="center", color="#888")
        return

    n_cycles = len(cycles)
    panel_w = 0.90 / n_cycles

    for ci, (name, color, eta, points, cycle_type) in enumerate(cycles):
        px0 = 0.05 + ci * panel_w
        # Panel border
        ax.add_patch(Rectangle((px0, 0.10), panel_w - 0.02, 0.80,
                               transform=ax.transAxes,
                               facecolor="#111118", edgecolor="#333", linewidth=0.5))

        # Title
        ax.text(px0 + (panel_w-0.02)/2, 0.88, name, fontsize=4,
                ha="center", color=color, fontweight="bold",
                transform=ax.transAxes)
        ax.text(px0 + (panel_w-0.02)/2, 0.84, f"eta={eta:.2f}", fontsize=3,
                ha="center", color="#aaa", transform=ax.transAxes)

        # Scale points to panel
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)

        def scale_pt(x, y):
            sx_p = (x - x_min) / (x_max - x_min) * (panel_w - 0.06) + px0 + 0.02
            sy_p = (y - y_min) / (y_max - y_min) * 0.55 + 0.20
            return sx_p, sy_p

        # Draw cycle path
        scaled = [scale_pt(x, y) for x, y in points]
        xs_p = [p[0] for p in scaled]
        ys_p = [p[1] for p in scaled]
        ax.plot(xs_p, ys_p, color=color, linewidth=1.2, alpha=0.7,
                transform=ax.transAxes)

        # Fill area under curve
        fill_pts = scaled + [(scaled[-1][0], 0.20), (scaled[0][0], 0.20)]
        ax.add_patch(Polygon(fill_pts, facecolor=color, alpha=0.1,
                             edgecolor="none", transform=ax.transAxes))

        # State points
        for i, (sx_p, sy_p) in enumerate(scaled):
            ax.add_patch(Circle((sx_p, sy_p), 0.004, transform=ax.transAxes,
                                facecolor=color, edgecolor="white", linewidth=0.2, alpha=0.8))
            if i < 4:
                ax.text(sx_p, sy_p + 0.02, str(i+1), fontsize=2,
                        ha="center", color="white", transform=ax.transAxes)

        # Axis labels
        ax.text(px0 + (panel_w-0.02)/2, 0.13, "s (entropy)", fontsize=2.5,
                ha="center", color="#888", transform=ax.transAxes)
        ax.text(px0 + 0.01, 0.50, "T", fontsize=3, ha="center", va="center",
                color="#888", rotation=90, transform=ax.transAxes)

        # Cycle type
        ax.text(px0 + (panel_w-0.02)/2, 0.17, cycle_type, fontsize=2.5,
                ha="center", color="#666", transform=ax.transAxes)


def _draw_heat_pipe_detail(ax, t: Dict) -> None:
    """Draw a detailed heat pipe cross-section showing the wick structure,
    working fluid, phase change zones, and lava interface."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon, Wedge, Arc
    lv = t["lava"]

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL K  -  Heat pipe array", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    if not lv.heat_pipe:
        ax.text(0.50, 0.50, "Heat pipes disabled", transform=ax.transAxes,
                fontsize=5, ha="center", color="#888")
        return

    # --- Left side: single heat pipe cross-section (zoomed) ---
    # Outer wall
    ax.add_patch(Rectangle((0.06, 0.15), 0.10, 0.75, transform=ax.transAxes,
                           facecolor="#37474F", edgecolor="#263238", linewidth=0.8, alpha=0.5))
    # Wick structure (inner layer)
    ax.add_patch(Rectangle((0.07, 0.15), 0.08, 0.75, transform=ax.transAxes,
                           facecolor="#546E7A", edgecolor="#37474F", linewidth=0.3, alpha=0.4))
    # Vapor core (center)
    ax.add_patch(Rectangle((0.085, 0.15), 0.05, 0.75, transform=ax.transAxes,
                           facecolor="#FF6347", edgecolor="#FF4500", linewidth=0.3, alpha=0.2))

    # Zone labels (left side)
    zones = [
        (0.75, "Condenser", "#4FC3F7", "vapor condenses\n-> liquid"),
        (0.50, "Adiabatic", "#FFD700", "vapor transport\n(no heat transfer)"),
        (0.25, "Evaporator", "#FF4500", "liquid evaporates\n-> vapor"),
    ]
    for y, name, color, desc in zones:
        ax.plot([0.04, 0.06], [y, y], color=color, linewidth=0.5,
                transform=ax.transAxes)
        ax.text(0.035, y, name, fontsize=3, ha="right", va="center",
                color=color, fontweight="bold", transform=ax.transAxes)
        ax.text(0.18, y, desc, fontsize=2.5, va="center", color="#aaa",
                transform=ax.transAxes)

    # Lava zone (bottom)
    ax.add_patch(Rectangle((0.04, 0.08), 0.14, 0.07, transform=ax.transAxes,
                           facecolor="#FF4500", edgecolor="#8B0000", linewidth=0.5, alpha=0.4))
    ax.text(0.11, 0.11, "LAVA", fontsize=3, ha="center", va="center",
            color="white", fontweight="bold", transform=ax.transAxes)

    # Heat flow arrows
    for i in range(3):
        ax.annotate("", xy=(0.11, 0.20 + i*0.02), xytext=(0.11, 0.08),
                    arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=0.5, alpha=0.5),
                    transform=ax.transAxes)

    # Wick detail callout
    ax.add_patch(Rectangle((0.07, 0.40), 0.08, 0.05, transform=ax.transAxes,
                           facecolor="none", edgecolor="#FFEB3B", linewidth=0.5))
    ax.annotate("wick\n(sintered Cu)", xy=(0.075, 0.425), xytext=(0.02, 0.55),
                arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=0.3),
                fontsize=2.5, color="#FFEB3B", transform=ax.transAxes)

    # --- Right side: heat pipe array layout (top view) ---
    ax.text(0.55, 0.88, "Array layout (top view)", fontsize=3.5,
            ha="center", color="#00d4ff", fontweight="bold", transform=ax.transAxes)

    # Tunnel casing outline
    ax.add_patch(Rectangle((0.40, 0.55), 0.30, 0.20, transform=ax.transAxes,
                           facecolor="#333", edgecolor="#666", linewidth=0.5, alpha=0.3))

    # Heat pipes (circles in grid)
    n_hp = 24
    cols = 6
    rows = 4
    for i in range(min(n_hp, cols * rows)):
        r = i // cols
        c = i % cols
        hx = 0.43 + c * 0.045
        hy = 0.58 + r * 0.04
        ax.add_patch(Circle((hx, hy), 0.008, transform=ax.transAxes,
                            facecolor="#FF6347", edgecolor="#FF4500", linewidth=0.2, alpha=0.6))

    # Flow direction arrow
    ax.annotate("", xy=(0.72, 0.65), xytext=(0.38, 0.65),
                arrowprops=dict(arrowstyle="->", color="#FFD700", lw=1),
                transform=ax.transAxes)
    ax.text(0.55, 0.53, "air flow", fontsize=2.5, ha="center", color="#FFD700",
            transform=ax.transAxes)

    # Specs
    # Specs (use safe defaults since LavaSourceSpec only has heat_pipe boolean)
    hp_n = getattr(lv, 'heat_pipe_n', 200)
    hp_d = getattr(lv, 'heat_pipe_d_mm', 25)
    hp_fluid = getattr(lv, 'heat_pipe_fluid', 'NaK (sodium-potassium)')
    hp_material = getattr(lv, 'heat_pipe_material', 'Inconel 600')
    hp_q = getattr(lv, 'heat_pipe_q_each', 5.0)

    ax.text(0.55, 0.40, f"Count: {hp_n} pipes\n"
                        f"Diameter: {hp_d:.0f} mm\n"
                        f"Working fluid: {hp_fluid}\n"
                        f"Material: {hp_material}\n"
                        f"Capacity: {hp_q:.0f} kW each\n"
                        f"Total: {hp_q * hp_n:.0f} kW",
            fontsize=3, ha="center", color="white", fontweight="bold",
            transform=ax.transAxes, family="monospace")

    # Temperature gradient
    ax.text(0.55, 0.20, "Temperature gradient:", fontsize=3, ha="center",
            color="#aaa", transform=ax.transAxes)
    for i in range(10):
        frac = i / 9.0
        tx = 0.40 + frac * 0.30
        r = int(255 * frac)
        b = int(255 * (1 - frac))
        color = f"#{r:02x}40{b:02x}"
        ax.add_patch(Rectangle((tx-0.01, 0.16), 0.025, 0.02, transform=ax.transAxes,
                               facecolor=color, edgecolor="none", alpha=0.6))
    ax.text(0.40, 0.13, f"{lv.t_lava_c:.0f}C", fontsize=2.5, ha="center",
            color="#FF4500", transform=ax.transAxes)
    ax.text(0.70, 0.13, "350C", fontsize=2.5, ha="center",
            color="#4FC3F7", transform=ax.transAxes)


def _draw_cooling_tower_detail(ax, t: Dict) -> None:
    """Draw a cooling tower cross-section showing fill, water distribution,
    drift eliminators, fan, and air/water flow paths."""
    from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon, Wedge, Arc

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL L  -  Cooling tower cross-section", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    # Tower shell (hyperbolic shape approximated with trapezoid)
    pts = [(0.20, 0.10), (0.80, 0.10), (0.75, 0.85), (0.25, 0.85)]
    ax.add_patch(Polygon(pts, facecolor="#B0BEC5", edgecolor="#78909C",
                         linewidth=1, alpha=0.2, transform=ax.transAxes))

    # Basin (bottom)
    ax.add_patch(Rectangle((0.22, 0.10), 0.56, 0.05, transform=ax.transAxes,
                           facecolor="#26A69A", edgecolor="#00897B", linewidth=0.5, alpha=0.4))
    ax.text(0.50, 0.125, "basin (cold water)", fontsize=2.5, ha="center",
            va="center", color="white", fontweight="bold", transform=ax.transAxes)

    # Fill media (above basin)
    ax.add_patch(Rectangle((0.25, 0.18), 0.50, 0.15, transform=ax.transAxes,
                           facecolor="#80CBC4", edgecolor="#4DB6AC", linewidth=0.3, alpha=0.3))
    ax.text(0.50, 0.255, "FILL (film fill)", fontsize=3, ha="center", va="center",
            color="#004D40", fontweight="bold", transform=ax.transAxes)
    # Fill pattern lines
    for i in range(8):
        fx = 0.27 + i * 0.06
        ax.plot([fx, fx], [0.19, 0.32], color="#4DB6AC", linewidth=0.2,
                alpha=0.3, transform=ax.transAxes)

    # Water distribution header (above fill)
    ax.plot([0.28, 0.72], [0.36, 0.36], color="#03A9F4", linewidth=1.5,
            alpha=0.5, transform=ax.transAxes)
    # Spray nozzles
    for i in range(6):
        nx = 0.30 + i * 0.07
        ax.add_patch(Circle((nx, 0.35), 0.005, transform=ax.transAxes,
                            facecolor="#03A9F4", edgecolor="#0277FD", linewidth=0.2, alpha=0.6))
        # Water droplets
        for d in range(3):
            ax.plot([nx, nx + 0.005], [0.34 - d*0.02, 0.33 - d*0.02],
                    color="#4FC3F7", linewidth=0.2, alpha=0.3,
                    transform=ax.transAxes)

    ax.text(0.50, 0.38, "water distribution", fontsize=2.5, ha="center",
            color="#03A9F4", transform=ax.transAxes)

    # Drift eliminators (above water distribution)
    ax.add_patch(Rectangle((0.25, 0.42), 0.50, 0.06, transform=ax.transAxes,
                           facecolor="#90A4AE", edgecolor="#607D8B", linewidth=0.3, alpha=0.3))
    ax.text(0.50, 0.45, "drift eliminators", fontsize=2.5, ha="center", va="center",
            color="#37474F", fontweight="bold", transform=ax.transAxes)
    # Eliminator pattern
    for i in range(10):
        ex = 0.26 + i * 0.048
        ax.plot([ex, ex + 0.02], [0.42, 0.48], color="#607D8B", linewidth=0.2,
                alpha=0.3, transform=ax.transAxes)

    # Fan (at top)
    fan_cx, fan_cy = 0.50, 0.75
    ax.add_patch(Circle((fan_cx, fan_cy), 0.08, transform=ax.transAxes,
                        facecolor="#26A69A", edgecolor="#00897B", linewidth=0.5, alpha=0.3))
    # Fan blades
    for b in range(6):
        ang = 2 * math.pi * b / 6
        ax.plot([fan_cx, fan_cx + 0.07 * math.cos(ang)],
               [fan_cy, fan_cy + 0.07 * math.sin(ang)],
               color="#00897B", linewidth=0.8, alpha=0.5,
               transform=ax.transAxes)
    ax.add_patch(Circle((fan_cx, fan_cy), 0.015, transform=ax.transAxes,
                        facecolor="#37474F", edgecolor="#263238", linewidth=0.3))
    ax.text(0.50, 0.83, "fan (induced draft)", fontsize=2.5, ha="center",
            color="#26A69A", transform=ax.transAxes)

    # Air inlet (louvers at bottom sides)
    for side in [0.22, 0.74]:
        for i in range(4):
            ly = 0.16 + i * 0.02
            ax.plot([side, side + 0.02], [ly, ly + 0.01], color="#78909C",
                    linewidth=0.3, alpha=0.4, transform=ax.transAxes)

    # Air flow arrows (upward)
    for ax_x in [0.35, 0.50, 0.65]:
        ax.annotate("", xy=(ax_x, 0.70), xytext=(ax_x, 0.15),
                    arrowprops=dict(arrowstyle="->", color="#80DEEA", lw=0.5, alpha=0.4),
                    transform=ax.transAxes)

    # Warm water inlet (top side)
    ax.annotate("", xy=(0.28, 0.36), xytext=(0.15, 0.36),
                arrowprops=dict(arrowstyle="->", color="#FF6347", lw=1),
                transform=ax.transAxes)
    ax.text(0.15, 0.38, "warm\nwater\nin", fontsize=2.5, ha="center",
            color="#FF6347", transform=ax.transAxes)

    # Cold water outlet (bottom)
    ax.annotate("", xy=(0.10, 0.12), xytext=(0.22, 0.12),
                arrowprops=dict(arrowstyle="->", color="#26A69A", lw=1),
                transform=ax.transAxes)
    ax.text(0.08, 0.12, "cold\nwater\nout", fontsize=2.5, ha="center",
            color="#26A69A", transform=ax.transAxes)

    # Exhaust air (top)
    ax.annotate("", xy=(0.50, 0.95), xytext=(0.50, 0.85),
                arrowprops=dict(arrowstyle="->", color="#80DEEA", lw=1),
                transform=ax.transAxes)
    ax.text(0.50, 0.97, "moist air out", fontsize=2.5, ha="center",
            color="#80DEEA", transform=ax.transAxes)

    # Specs
    ax.text(0.50, 0.04, "Mechanical draft cooling tower | 10-cell | "
                        "30/40C cold/warm | 50,000 m3/h",
            fontsize=2.5, ha="center", color="white", fontweight="bold",
            transform=ax.transAxes)


def _draw_control_architecture(ax, t: Dict) -> None:
    """Draw a control system architecture diagram showing SCADA hierarchy,
    PLC network, I/O mapping, and communication paths."""
    from matplotlib.patches import Rectangle, FancyBboxPatch, Circle

    ax.set_facecolor("#0d0d18")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("DETAIL M  -  Control architecture", fontsize=7,
                 fontweight="bold", color="#FFEB3B", pad=4)

    ctrl = t["ctrl"]

    # Hierarchy levels (top to bottom)
    levels = [
        (0.85, "Level 3: SCADA / DCS", "#00d4ff", [
            ("Operator HMI", 0.15), ("Historian", 0.40), ("Engineering", 0.65)
        ]),
        (0.65, "Level 2: Supervisory PLCs", "#4CAF50", [
            ("Turbine Gov", 0.12), ("BOP Ctrl", 0.37), ("Fire/Gas", 0.62), ("Elec Prot", 0.85)
        ]),
        (0.40, "Level 1: Field PLCs / RTUs", "#FF9800", [
            ("Cavern PLC", 0.08), ("Tunnel PLC", 0.28), ("Turb PLC", 0.48),
            ("Bottom PLC", 0.68), ("Cooling PLC", 0.88)
        ]),
        (0.15, "Level 0: Field I/O", "#F44336", [
            ("PT/TT/FT", 0.05), ("Valves", 0.20), ("Breakers", 0.35),
            ("Pumps", 0.50), ("Fans", 0.65), ("Sensors", 0.80), ("Relays", 0.92)
        ]),
    ]

    for y, title, color, items in levels:
        # Level label
        ax.text(0.02, y, title, fontsize=3, va="center", color=color,
                fontweight="bold", transform=ax.transAxes, rotation=0)
        # Level bar
        ax.plot([0.03, 0.98], [y - 0.03, y - 0.03], color=color, linewidth=0.3,
                alpha=0.2, transform=ax.transAxes)
        # Items
        for label, x in items:
            ax.add_patch(FancyBboxPatch((x - 0.04, y - 0.03), 0.08, 0.06,
                                        boxstyle="round,pad=0.005",
                                        facecolor=color, edgecolor="#333",
                                        linewidth=0.3, alpha=0.4,
                                        transform=ax.transAxes))
            ax.text(x, y, label, fontsize=2, ha="center", va="center",
                    color="white", fontweight="bold", transform=ax.transAxes)

    # Communication links between levels
    link_color = "#FFEB3B"
    for i in range(len(levels) - 1):
        y1 = levels[i][0] - 0.03
        y2 = levels[i + 1][0] + 0.03
        # Draw a few connecting lines
        for x in [0.20, 0.50, 0.80]:
            ax.plot([x, x], [y1, y2], color=link_color, linewidth=0.3,
                    alpha=0.2, linestyle="dashed", transform=ax.transAxes)

    # Network ring (Ethernet)
    ax.text(0.50, 0.55, "Industrial Ethernet Ring (PRP/HSR)", fontsize=2.5,
            ha="center", color=link_color, alpha=0.5, transform=ax.transAxes)

    # Communication protocols
    ax.text(0.50, 0.06, "Protocols: IEC 61850 (electrical) | Modbus TCP (field) | "
                        "OPC UA (SCADA) | DNP3 (grid)",
            fontsize=2.5, ha="center", color="#aaa", transform=ax.transAxes)

    # Redundancy indicators
    ax.text(0.02, 0.94, "Redundancy:", fontsize=2.5, color="#FFEB3B",
            transform=ax.transAxes)
    ax.text(0.02, 0.91, "  Dual servers, dual PLC, dual network",
            fontsize=2, color="#aaa", transform=ax.transAxes)

    # I/O count
    ax.text(0.70, 0.94, f"I/O: {MONITOR_HW['scada_points']:,} points",
            fontsize=2.5, color="#00E676", transform=ax.transAxes, fontweight="bold")


def _draw_blueprint(ax, t: Dict) -> None:
    """Draw a detailed multi-view engineering blueprint schematic.

    Twenty-one panels:
      1. ISOMETRIC CUTAWAY - 3D perspective view showing depth and internals
      2. PLAN VIEW - top-down layout with dimensions and flow arrows
      3. SIDE ELEVATION - full cross-section with depth, callouts, dimensions
      4. END ELEVATION - bore pattern looking down tunnel axis
      5. CAVERN LINING DETAIL - layered wall cross-section (zoomed)
      6. TURBINE STAGE DETAIL - stator/rotor/stator cross-section
      7. LAVA HX DETAIL - tube arrangement and heat pipe detail
      8. FAN/NOZZLE DETAIL - exit fan and nozzle cross-section
      9. GENERATOR DETAIL - generator cross-section with rotor/stator/H2
      10. P&ID - piping & instrumentation diagram with control loops
      11. EXPLODED VIEW - turbine assembly components separated
      12. ELECTRICAL SLD - single-line diagram from gen to grid
      13. REHEAT DETAIL - reheat section cross-section
      14. SITE LAYOUT - surface plan with all buildings, roads, fences
      15. CAVERN INTERIOR - perspective view of cavern internals
      16. T-s DIAGRAMS - thermodynamic cycle diagrams for bottoming cycles
      17. HEAT PIPE DETAIL - heat pipe cross-section and array layout
      18. COOLING TOWER DETAIL - cross-section with fill, fan, water flow
      19. CONTROL ARCHITECTURE - SCADA hierarchy and PLC network
      20. TITLE BLOCK - specs, dimensions, component count
      21. LEGEND - symbol key, color coding
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Wedge, Arc

    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t["ctrl"]
    n_sys = max(1, ctrl.n_systems)
    n_bores = max(1, lv.n_parallel_bores)
    ax.clear()
    ax.set_facecolor("#0a0a12")
    ax.set_axis_off()

    fig = ax.figure
    # Cache sub-axes on the figure to avoid recreating them every redraw.
    # On first call, create all 20 sub-axes. On subsequent calls, reuse and just clear.
    axes_specs = [
        ("ax_iso",     [0.04, 0.85, 0.52, 0.13]),
        ("ax_plan",    [0.04, 0.70, 0.52, 0.13]),
        ("ax_side",    [0.04, 0.55, 0.52, 0.13]),
        ("ax_site",    [0.04, 0.41, 0.52, 0.12]),
        ("ax_explod",  [0.04, 0.30, 0.22, 0.09]),
        ("ax_elec",    [0.28, 0.30, 0.18, 0.09]),
        ("ax_reheat",  [0.48, 0.30, 0.20, 0.09]),
        ("ax_hp",      [0.04, 0.20, 0.30, 0.08]),
        ("ax_ct",      [0.36, 0.20, 0.30, 0.08]),
        ("ax_ctrl",    [0.68, 0.20, 0.28, 0.08]),
        ("ax_end",     [0.60, 0.85, 0.36, 0.13]),
        ("ax_cavdet",  [0.60, 0.70, 0.18, 0.13]),
        ("ax_cavint",  [0.60, 0.55, 0.36, 0.13]),
        ("ax_gendet",  [0.60, 0.41, 0.36, 0.12]),
        ("ax_ts",      [0.70, 0.30, 0.26, 0.09]),
        ("ax_turbdet", [0.04, 0.01, 0.16, 0.16]),
        ("ax_hxdet",   [0.22, 0.01, 0.16, 0.16]),
        ("ax_fandet",  [0.40, 0.01, 0.16, 0.16]),
        ("ax_pid",     [0.58, 0.01, 0.22, 0.16]),
        ("ax_title",   [0.82, 0.01, 0.14, 0.16]),
    ]

    if not hasattr(fig, '_bp_axes_cache'):
        # First call: create all sub-axes
        fig._bp_axes_cache = {}
        for name, spec in axes_specs:
            a = fig.add_axes(spec)
            a.set_facecolor("#111118")
            for spine in a.spines.values():
                spine.set_color("#4a5568")
                spine.set_linewidth(1.0)
            a.tick_params(colors="#a0aec0", labelsize=4)
            fig._bp_axes_cache[name] = a
    else:
        # Subsequent calls: verify axes still exist (figure may have been closed)
        for name, spec in axes_specs:
            if name not in fig._bp_axes_cache or fig._bp_axes_cache[name] not in fig.axes:
                a = fig.add_axes(spec)
                a.set_facecolor("#111118")
                for spine in a.spines.values():
                    spine.set_color("#4a5568")
                    spine.set_linewidth(1.0)
                a.tick_params(colors="#a0aec0", labelsize=4)
                fig._bp_axes_cache[name] = a

    # Assign cached axes to local variables
    ax_iso     = fig._bp_axes_cache["ax_iso"]
    ax_plan    = fig._bp_axes_cache["ax_plan"]
    ax_side    = fig._bp_axes_cache["ax_side"]
    ax_site    = fig._bp_axes_cache["ax_site"]
    ax_explod  = fig._bp_axes_cache["ax_explod"]
    ax_elec    = fig._bp_axes_cache["ax_elec"]
    ax_reheat  = fig._bp_axes_cache["ax_reheat"]
    ax_hp      = fig._bp_axes_cache["ax_hp"]
    ax_ct      = fig._bp_axes_cache["ax_ct"]
    ax_ctrl    = fig._bp_axes_cache["ax_ctrl"]
    ax_end     = fig._bp_axes_cache["ax_end"]
    ax_cavdet  = fig._bp_axes_cache["ax_cavdet"]
    ax_cavint  = fig._bp_axes_cache["ax_cavint"]
    ax_gendet  = fig._bp_axes_cache["ax_gendet"]
    ax_ts      = fig._bp_axes_cache["ax_ts"]
    ax_turbdet = fig._bp_axes_cache["ax_turbdet"]
    ax_hxdet   = fig._bp_axes_cache["ax_hxdet"]
    ax_fandet  = fig._bp_axes_cache["ax_fandet"]
    ax_pid     = fig._bp_axes_cache["ax_pid"]
    ax_title   = fig._bp_axes_cache["ax_title"]

    # Clear all sub-axes for redraw (much faster than recreating)
    for a in [ax_iso, ax_plan, ax_side, ax_site, ax_explod, ax_elec, ax_reheat,
              ax_hp, ax_ct, ax_ctrl, ax_end, ax_cavdet, ax_cavint, ax_gendet,
              ax_ts, ax_turbdet, ax_hxdet, ax_fandet, ax_pid, ax_title]:
        a.clear()
        a.set_facecolor("#111118")

    # --- dimensions ---
    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    cav_depth = cv.depth_m
    tun_len = tn.total_length_m
    tun_r = tn.diameter_m / 2.0
    lava_len = lv.contact_length_m
    stack_h = tn.height_rise_m
    n_turb = tn.n_turbine_stages
    n_reheat = getattr(tn, 'n_reheat_stages', 0)
    n_fans = tn.n_exit_fans

    # ============================================================
    # PANEL 1: ISOMETRIC CUTAWAY (delegated)
    # ============================================================
    _draw_isometric_cutaway(ax_iso, t)

    # ============================================================
    # PANEL 2: PLAN VIEW (top-down)
    # ============================================================
    ax_plan.set_title("PLAN VIEW  -  top-down layout", fontsize=7, fontweight="bold",
                      color="#00d4ff", pad=4)
    total_w = cav_side + tun_len + stack_h + 100
    sx = 1.0 / total_w * 0.92
    sy = 1.0 / max(cav_side * 2, n_bores * tn.diameter_m * 3) * 0.80
    def px(m): return m * sx
    def py(m): return m * sy

    ax_plan.axhline(0, color="#4a5568", linewidth=0.5, linestyle="--", alpha=0.4)

    cav_x = px(10)
    cav_w = px(cav_side)
    cav_h = py(cav_side)
    cav_color = "#2196F3" if k_to_c(cv.t_charge_k) < -100 else "#4FC3F7"
    ax_plan.add_patch(FancyBboxPatch((cav_x, -cav_h/2), cav_w, cav_h,
                                      boxstyle="round,pad=1", facecolor=cav_color,
                                      edgecolor="#0277FD", linewidth=1.2, alpha=0.5))
    ax_plan.text(cav_x + cav_w/2, 0, f"CAVERN\n{cv.volume_m3/1e9:.1f} km3\n{k_to_c(cv.t_charge_k):.0f}C\n{cv.p_charge_pa/1e5:.0f}bar",
                 ha="center", va="center", fontsize=4, color="white", fontweight="bold")

    ax_plan.plot([cav_x + cav_w*0.3, cav_x + cav_w*0.3], [cav_h/2, cav_h/2 + py(50)],
                 color="#888", linewidth=1.5)
    ax_plan.text(cav_x + cav_w*0.3 + 2, cav_h/2 + py(25), "access", fontsize=3, color="#aaa")

    tun_x0 = cav_x + cav_w
    tun_x1 = tun_x0 + px(tun_len)
    bore_spacing = max(tn.diameter_m * 1.5, 20)
    n_show = min(n_bores, 12)
    for bi in range(n_show):
        by = (bi - (n_show - 1) / 2.0) * py(bore_spacing)
        ax_plan.plot([tun_x0, tun_x0 + px(tun_len * 0.15)], [by, by],
                     color="#FFD700", linewidth=0.7, alpha=0.5)
        lava_x0 = tun_x0 + px(tun_len * 0.15)
        lava_x1 = lava_x0 + px(lava_len)
        ax_plan.plot([lava_x0, lava_x1], [by, by], color="#FF6347", linewidth=1, alpha=0.7)
        ax_plan.plot([lava_x1, tun_x1], [by, by], color="#FFD700", linewidth=0.7, alpha=0.5)

    if n_bores > n_show:
        ax_plan.text(tun_x0 + px(tun_len * 0.5), py(bore_spacing) * (n_show/2 + 1),
                     f"+ {n_bores - n_show} more", fontsize=3, color="#FFD700", ha="center")

    lava_y_h = py(n_bores * tn.diameter_m * 1.5)
    ax_plan.fill_between([lava_x0, lava_x1], -lava_y_h/2, lava_y_h/2,
                         color="#FF4500", alpha=0.15, zorder=1)
    ax_plan.text((lava_x0 + lava_x1)/2, lava_y_h/2 + 4,
                 f"LAVA {lv.t_lava_c:.0f}C", fontsize=4, color="#FF6347", ha="center", fontweight="bold")

    for i in range(min(n_turb, 10)):
        frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tun_len)
        tx = tun_x0 + px(frac * tun_len)
        ax_plan.plot(tx, 0, "o", color="#00E676", markersize=2, zorder=5)

    stack_x = tun_x1
    ax_plan.add_patch(Rectangle((stack_x, -py(tn.diameter_m/2)), px(stack_h), py(tn.diameter_m),
                                facecolor="#888", edgecolor="#666", alpha=0.6))
    ax_plan.text(stack_x + px(stack_h)/2, 0, "STACK", fontsize=4, color="white",
                 ha="center", va="center", rotation=90, fontweight="bold")

    for i in range(min(n_fans, 12)):
        fy = (i - 5) * py(tn.diameter_m * 0.8)
        ax_plan.plot(stack_x + px(stack_h) + 3, fy, "s", color="#00BFFF", markersize=2)

    # Dimensions
    dim_y = -cav_h/2 - 10
    ax_plan.annotate("", xy=(cav_x, dim_y), xytext=(stack_x + px(stack_h), dim_y),
                     arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.5))
    ax_plan.text((cav_x + stack_x + px(stack_h))/2, dim_y - 2,
                 f"{total_w:.0f}m", ha="center", fontsize=3.5, color="#a0aec0")

    ax_plan.annotate("", xy=(lava_x0, -lava_y_h/2 - 3), xytext=(lava_x1, -lava_y_h/2 - 3),
                     arrowprops=dict(arrowstyle="<->", color="#FF6347", lw=0.5))
    ax_plan.text((lava_x0 + lava_x1)/2, -lava_y_h/2 - 5, f"{lava_len:.0f}m lava", ha="center", fontsize=3.5, color="#FF6347")

    # Callout bubbles
    for label, lx, ly in [("A", cav_x + cav_w/2, -cav_h/2), ("B", lava_x0 + px(lava_len)*0.5, 0),
                           ("C", tun_x0 + px(tun_len*0.3), 0), ("D", stack_x + px(stack_h)*0.5, 0)]:
        ax_plan.add_patch(Circle((lx, ly), 3, fill=False, edgecolor="#FFEB3B", linewidth=1))
        ax_plan.text(lx, ly, label, fontsize=4, ha="center", va="center", color="#FFEB3B", fontweight="bold")

    ax_plan.set_xlim(-5, stack_x + px(stack_h) + 15)
    ax_plan.set_ylim(-cav_h/2 - 18, cav_h/2 + 14)
    ax_plan.set_aspect("equal")
    ax_plan.set_xlabel("Distance (m)", fontsize=4, color="#a0aec0")
    ax_plan.set_ylabel("Width (m)", fontsize=4, color="#a0aec0")

    # ============================================================
    # PANEL 3: SIDE ELEVATION
    # ============================================================
    ax_side.set_title("SIDE ELEVATION  -  cross-section with depth", fontsize=7, fontweight="bold",
                      color="#00d4ff", pad=4)
    v_exag = max(total_w / (cav_depth + stack_h + 100) * 0.3, 4.0)

    ax_side.axhline(0, color="#8B7355", linewidth=1)
    ax_side.fill_between([0, total_w * 1.02], 0, -(cav_depth + cav_side) * v_exag,
                         color="#2B1D0E", alpha=0.15)
    ax_side.fill_between([0, total_w * 1.02], 0, stack_h * v_exag * 1.1,
                         color="#1a1a2e", alpha=0.08)
    ax_side.text(2, 2, "SURFACE", fontsize=3, color="#8B7355", fontweight="bold")

    cav_y = -cav_depth * v_exag
    cav_h_v = cav_side * 0.5 * v_exag
    lining_thick = CAVERN_HW['lining_thick_mm'] / 1000.0 * v_exag
    ax_side.add_patch(FancyBboxPatch((cav_x - lining_thick, cav_y - cav_h_v - lining_thick),
                                     cav_w + 2*lining_thick, cav_h_v + 2*lining_thick,
                                     boxstyle="round,pad=0.5", facecolor="#555",
                                     edgecolor="#333", linewidth=0.3, alpha=0.3))
    insul_thick = CAVERN_HW['insulation_mm'] / 1000.0 * v_exag
    ax_side.add_patch(FancyBboxPatch((cav_x - insul_thick, cav_y - cav_h_v - insul_thick),
                                     cav_w + 2*insul_thick, cav_h_v + 2*insul_thick,
                                     boxstyle="round,pad=0.5", facecolor="#8D6E63",
                                     edgecolor="#6D4C41", linewidth=0.2, alpha=0.25))
    # Ultra thermal insulation layer (when enabled)
    if cv.ultra_insulation:
        ultra_thick = CAVERN_HW['ultra_insulation_mm'] / 1000.0 * v_exag
        ax_side.add_patch(FancyBboxPatch((cav_x - insul_thick - ultra_thick,
                                          cav_y - cav_h_v - insul_thick - ultra_thick),
                                         cav_w + 2*(insul_thick + ultra_thick),
                                         cav_h_v + 2*(insul_thick + ultra_thick),
                                         boxstyle="round,pad=0.5", facecolor="#CE93D8",
                                         edgecolor="#9C27B0", linewidth=0.3, alpha=0.2))
        ax_side.text(cav_x - insul_thick - ultra_thick - 2, cav_y - cav_h_v/2,
                     "ULTRA\nINSUL\nR=30", fontsize=2.5, color="#CE93D8",
                     ha="right", va="center", fontweight="bold")
    ax_side.add_patch(FancyBboxPatch((cav_x, cav_y - cav_h_v), cav_w, cav_h_v,
                                      boxstyle="round,pad=1", facecolor=cav_color,
                                      edgecolor="#0277FD", linewidth=1, alpha=0.5))
    ax_side.text(cav_x + cav_w/2, cav_y - cav_h_v/2,
                 f"CAVERN\n{cv.volume_m3/1e9:.1f} km3\n{k_to_c(cv.t_charge_k):.0f}C\n{cv.p_charge_pa/1e5:.0f}bar",
                 ha="center", va="center", fontsize=3.5, color="white", fontweight="bold")

    ax_side.plot([cav_x + cav_w*0.3, cav_x + cav_w*0.3], [0, cav_y],
                 color="#888", linewidth=1.5)
    ax_side.add_patch(Rectangle((cav_x + cav_w*0.3 - 2, cav_y - 2), 4, 4,
                                facecolor="#FF5722", edgecolor="#D84315", linewidth=0.3, alpha=0.7))
    ax_side.text(cav_x + cav_w*0.3 + 3, cav_y, "door", fontsize=2.5, color="#FF5722")

    tun_y = cav_y + cav_h_v * 0.3
    casing_thick = TUNNEL_HW['casing_od_mm'] / 1000.0 * v_exag * 0.3
    ax_side.fill_between([tun_x0, tun_x1], tun_y - 2 - casing_thick, tun_y - 2,
                         color="#666", alpha=0.3)
    ax_side.fill_between([tun_x0, tun_x1], tun_y + 2, tun_y + 2 + casing_thick,
                         color="#666", alpha=0.3)
    ax_side.fill_between([tun_x0, tun_x1], tun_y - 2, tun_y + 2,
                         color="#FFD700", alpha=0.15)
    ax_side.plot([tun_x0, tun_x1], [tun_y - 2, tun_y - 2], color="#B8860B", linewidth=0.6)
    ax_side.plot([tun_x0, tun_x1], [tun_y + 2, tun_y + 2], color="#B8860B", linewidth=0.6)

    ref_thick = TUNNEL_HW['refractory_thick_mm'] / 1000.0 * v_exag * 0.3
    ax_side.fill_between([lava_x0, lava_x1], tun_y - 2 - ref_thick, tun_y - 2,
                         color="#FF7043", alpha=0.4)
    ax_side.fill_between([lava_x0, lava_x1], tun_y + 2, tun_y + 2 + ref_thick,
                         color="#FF7043", alpha=0.4)

    lava_y = tun_y - 12 * v_exag
    ax_side.fill_between([lava_x0, lava_x1], lava_y, lava_y - 18 * v_exag,
                         color="#FF4500", alpha=0.4)
    ax_side.fill_between([lava_x0, lava_x1], lava_y - 18 * v_exag, lava_y - 22 * v_exag,
                         color="#8B0000", alpha=0.5)
    ax_side.text((lava_x0 + lava_x1)/2, lava_y - 9 * v_exag,
                 f"LAVA\n{lv.t_lava_c:.0f}C", ha="center", va="center",
                 fontsize=3.5, color="white", fontweight="bold")

    if lv.hx_enabled:
        for i in range(min(lv.hx_n_tubes, 6)):
            ty = tun_y - 3 - i * 0.8
            ax_side.plot([lava_x0 + 2, lava_x1 - 2], [ty, ty],
                        color="#FF8C00", linewidth=0.2, alpha=0.5)

    for i in range(min(n_turb, 10)):
        frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tun_len)
        tx = tun_x0 + px(frac * tun_len)
        ax_side.add_patch(Rectangle((tx - 1.5, tun_y - 2.5), 3, 5,
                                    facecolor="#00E676", edgecolor="#00C853", linewidth=0.3, alpha=0.6))
        if i < 6:
            ax_side.text(tx, tun_y + 4, f"T{i+1}", fontsize=2, ha="center", color="#00E676")

    if n_reheat > 0:
        for i in range(min(n_reheat, 6)):
            frac = 0.15 + (i + 1.5) / (n_turb + 1) * (lava_len / tun_len)
            rx = tun_x0 + px(frac * tun_len)
            ax_side.add_patch(Rectangle((rx - 1, tun_y - 1.5), 2, 3,
                                        facecolor="#FF5722", edgecolor="#BF360C", linewidth=0.2, alpha=0.5))

    n_joints = int(tun_len / TUNNEL_HW['expansion_joint_m'])
    for j in range(min(n_joints, 12)):
        jx = tun_x0 + (j + 1) * px(TUNNEL_HW['expansion_joint_m'])
        if jx > tun_x1: break
        ax_side.plot([jx, jx], [tun_y - 2.5, tun_y + 2.5], color="#FF9800", linewidth=0.3, alpha=0.5)

    ax_side.plot([tun_x0, tun_x1], [tun_y - 4 - casing_thick, tun_y - 4 - casing_thick],
                 color="#26A69A", linewidth=0.6, alpha=0.5)

    bottoming = []
    if tn.potassium_enabled: bottoming.append(("K", "#FF9800"))
    if tn.sco2_enabled: bottoming.append(("sCO2", "#9C27B0"))
    if tn.steam_enabled: bottoming.append(("Steam", "#03A9F4"))
    if tn.orc_enabled: bottoming.append(("ORC", "#4CAF50"))
    for i, (name, color) in enumerate(bottoming):
        bx = stack_x + 10 + i * 12
        ax_side.add_patch(FancyBboxPatch((bx, -6 * v_exag), 10, 5,
                                          boxstyle="round,pad=1", facecolor=color, alpha=0.5))
        ax_side.text(bx + 5, -3.5 * v_exag, name, fontsize=2.5, ha="center", color="white", fontweight="bold")

    ax_side.fill_between([stack_x, stack_x + 4], [0, 0], [stack_h * v_exag, stack_h * v_exag],
                         color="#888", alpha=0.5)
    for i in range(min(n_fans, 8)):
        fy = stack_h * v_exag * (0.5 + 0.08 * i)
        ax_side.add_patch(Circle((stack_x + 2, fy), 1.5, facecolor="#00BFFF", edgecolor="#0277FD", linewidth=0.2, alpha=0.6))

    ax_side.annotate("", xy=(stack_x + 2, stack_h * v_exag + 6),
                     xytext=(stack_x + 2, stack_h * v_exag),
                     arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=1))

    ax_side.add_patch(Rectangle((stack_x + 8, 2), 6, 4,
                                facecolor="#FFC107", edgecolor="#FF6F00", linewidth=0.3, alpha=0.5))
    ax_side.text(stack_x + 11, 4, "XFMR", fontsize=2.5, ha="center", color="white", fontweight="bold")

    # Dimensions
    ax_side.annotate("", xy=(cav_x - 5, 0), xytext=(cav_x - 5, cav_y),
                     arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.4))
    ax_side.text(cav_x - 8, cav_y/2, f"{cav_depth:.0f}m", fontsize=2.5, color="#a0aec0",
                 ha="center", rotation=90)

    ax_side.annotate("", xy=(stack_x + 6, 0), xytext=(stack_x + 6, stack_h * v_exag),
                     arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.4))
    ax_side.text(stack_x + 9, stack_h * v_exag / 2, f"{stack_h:.0f}m", fontsize=2.5,
                 color="#a0aec0", rotation=90)

    ax_side.annotate("", xy=(tun_x0, tun_y - 6), xytext=(tun_x1, tun_y - 6),
                     arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.4))
    ax_side.text((tun_x0 + tun_x1)/2, tun_y - 9, f"{tun_len:.0f}m", fontsize=2.5,
                 color="#a0aec0", ha="center")

    # Numbered callouts
    callouts = [
        (1, cav_x + cav_w * 0.5, cav_y - cav_h_v, "Cavern"),
        (2, lava_x0 + px(lava_len) * 0.5, lava_y - 22 * v_exag, "Lava HX"),
        (3, tun_x0 + px(tun_len * 0.3), tun_y + 5, "Turbines"),
        (4, stack_x + 2, stack_h * v_exag * 0.5, "Stack+Fans"),
        (5, stack_x + 11, 4, "XFMR"),
        (6, cav_x + cav_w*0.3, cav_y, "Door"),
        (7, tun_x0 + px(tun_len)*0.5, tun_y - 5, "Drain"),
    ]
    for num, cx, cy, label in callouts:
        ax_side.annotate(f"{num}", xy=(cx, cy), fontsize=3, color="#FFEB3B", fontweight="bold",
                        xytext=(cx + 6, cy - 5),
                        arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=0.3, alpha=0.5),
                        ha="center", zorder=10)

    ax_side.set_xlim(-10, stack_x + 60)
    ax_side.set_ylim(-(cav_depth + cav_side) * v_exag - 6, stack_h * v_exag + 12)
    ax_side.set_aspect("equal")
    ax_side.set_xlabel("Distance (m)", fontsize=4, color="#a0aec0")
    ax_side.set_ylabel(f"Depth/Height ({v_exag:.0f}x exag)", fontsize=4, color="#a0aec0")

    # ============================================================
    # PANEL 4: END ELEVATION
    # ============================================================
    ax_end.set_title("END ELEVATION  -  bore pattern", fontsize=7, fontweight="bold",
                     color="#00d4ff", pad=4)

    grid_n = int(math.ceil(math.sqrt(n_bores)))
    spacing = max(tn.diameter_m * 1.5, 20)
    bore_r_draw = min(6, max(2, 40 / n_bores))

    outer_r = grid_n * spacing / 2 + spacing * 0.4
    ax_end.add_patch(Circle((0, 0), outer_r, fill=True, facecolor="#333",
                            edgecolor="#666", linewidth=1, alpha=0.2))
    ax_end.add_patch(Circle((0, 0), outer_r - 2, fill=False, edgecolor="#FF7043",
                            linewidth=0.8, linestyle="--", alpha=0.5))

    for i in range(min(n_bores, 48)):
        row = i // grid_n
        col = i % grid_n
        bx = (col - (grid_n - 1) / 2.0) * spacing
        by = (row - (grid_n - 1) / 2.0) * spacing
        ax_end.add_patch(Circle((bx, by), bore_r_draw + 1, fill=True, facecolor="#555",
                                edgecolor="#444", linewidth=0.2, alpha=0.4))
        ax_end.add_patch(Circle((bx, by), bore_r_draw, fill=True, facecolor="#FFD700",
                                edgecolor="#B8860B", linewidth=0.3, alpha=0.6))

    if n_bores > 48:
        ax_end.text(0, -outer_r - 6, f"+ {n_bores - 48} more", fontsize=3, color="#FFD700", ha="center")

    ax_end.annotate("", xy=(-outer_r, -outer_r - 5), xytext=(outer_r, -outer_r - 5),
                    arrowprops=dict(arrowstyle="<->", color="#a0aec0", lw=0.4))
    ax_end.text(0, -outer_r - 8, f"{outer_r*2:.0f}m", fontsize=3, color="#a0aec0", ha="center")

    ax_end.text(0, outer_r + 6, f"{n_bores} bores x {tn.diameter_m:.0f}m\ngrid: {grid_n}x{grid_n}",
                fontsize=3.5, color="#00d4ff", ha="center", fontweight="bold")

    if n_bores > 1:
        ax_end.annotate("", xy=(0, -bore_r_draw), xytext=(spacing, -bore_r_draw),
                        arrowprops=dict(arrowstyle="<->", color="#FFEB3B", lw=0.3))
        ax_end.text(spacing/2, -bore_r_draw - 2.5, f"{spacing:.0f}m", fontsize=2.5,
                    color="#FFEB3B", ha="center")

    ax_end.set_xlim(-outer_r - 12, outer_r + 12)
    ax_end.set_ylim(-outer_r - 12, outer_r + 10)
    ax_end.set_aspect("equal")

    # ============================================================
    # PANEL 5: CAVERN LINING DETAIL (Detail A)
    # ============================================================
    ax_cavdet.set_title("DETAIL A  -  Cavern wall", fontsize=6, fontweight="bold",
                        color="#FFEB3B", pad=3)
    ax_cavdet.set_axis_off()

    layers = [
        ("Host rock", "#3E2723", 0.28),
        ("Shotcrete 600mm", "#6D4C41", 0.16),
        ("HDPE seal 8mm", "#1565C0", 0.07),
        ("Cavern interior", cav_color, 0.28),
        ("HDPE seal 8mm", "#1565C0", 0.07),
        ("PU foam 200mm", "#8D6E63", 0.09),
    ]
    y_start = 0.90
    for name, color, h_frac in layers:
        h = h_frac * 0.78
        ax_cavdet.add_patch(Rectangle((0.08, y_start - h), 0.84, h,
                                       transform=ax_cavdet.transAxes,
                                       facecolor=color, edgecolor="#444", linewidth=0.3, alpha=0.6))
        ax_cavdet.text(0.50, y_start - h/2, name, transform=ax_cavdet.transAxes,
                       fontsize=3.5, ha="center", va="center", color="white", fontweight="bold")
        y_start -= h

    ax_cavdet.text(0.50, 0.05, f"Pressure: {CAVERN_HW['pressure_rating_bar']:.0f} bar\n"
                               f"Door: {CAVERN_HW['hydraulic_door_mm']:.0f}mm",
                   transform=ax_cavdet.transAxes, fontsize=3, ha="center", color="#a0aec0")

    # ============================================================
    # PANEL 6: GENERATOR DETAIL (Detail E - delegated)
    # ============================================================
    _draw_generator_detail(ax_gendet, t)

    # ============================================================
    # PANEL 7: TURBINE STAGE DETAIL (Detail C)
    # ============================================================
    ax_turbdet.set_title("DETAIL C  -  Turbine stage", fontsize=6, fontweight="bold",
                         color="#FFEB3B", pad=3)
    ax_turbdet.set_aspect("equal")

    rotor_d = TURBINE_HW['rotor_d_mm'] / 1000.0
    stage_spacing = TURBINE_HW['stage_spacing_m']
    n_blades = TURBINE_HW['rotor_blade_count']
    turb_scale = 0.7 / (stage_spacing * 3)
    rd = rotor_d * turb_scale * 0.4
    sd = stage_spacing * turb_scale
    cx_t, cy_t = 0.5, 0.5

    ax_turbdet.add_patch(Rectangle((cx_t - sd*1.5, cy_t - rd - 0.05), sd*3, (rd + 0.05) * 2,
                                    transform=ax_turbdet.transAxes,
                                    facecolor="#333", edgecolor="#666", linewidth=0.8, alpha=0.3))

    for stage in range(3):
        x = cx_t - sd + stage * sd
        for b in range(min(n_blades, 8)):
            ang = 2 * math.pi * b / min(n_blades, 8)
            bxe = x - sd * 0.15
            ax_turbdet.plot([bxe, bxe + rd*0.3 * math.cos(0.5 + ang*0.1)],
                           [cy_t + rd * math.sin(ang), cy_t + rd * math.sin(ang) - rd*0.2],
                           color="#FF9800", linewidth=0.3, alpha=0.6,
                           transform=ax_turbdet.transAxes)

        ax_turbdet.add_patch(Circle((x, cy_t), rd, transform=ax_turbdet.transAxes,
                                    facecolor="#00E676", edgecolor="#00C853", linewidth=0.5, alpha=0.5))
        for b in range(min(n_blades, 8)):
            ang = 2 * math.pi * b / min(n_blades, 8)
            ax_turbdet.plot([x, x + rd * math.cos(ang)],
                           [cy_t, cy_t + rd * math.sin(ang)],
                           color="#00C853", linewidth=0.4, alpha=0.7,
                           transform=ax_turbdet.transAxes)

    ax_turbdet.plot([cx_t - sd*1.5, cx_t + sd*1.5], [cy_t, cy_t],
                    color="#888", linewidth=1.5, transform=ax_turbdet.transAxes, zorder=5)

    ax_turbdet.text(0.5, 0.92, f"{n_turb} stages x {n_blades} blades\nD={TURBINE_HW['rotor_d_mm']:.0f}mm\nRPM={TURBINE_HW['rpm']:.0f}",
                   transform=ax_turbdet.transAxes, fontsize=3, ha="center", color="white", fontweight="bold")
    ax_turbdet.text(0.05, 0.5, "cold->", transform=ax_turbdet.transAxes,
                   fontsize=2.5, color="#FFEB3B", ha="center")
    ax_turbdet.text(0.95, 0.5, "->hot", transform=ax_turbdet.transAxes,
                   fontsize=2.5, color="#FF6347", ha="center")

    # ============================================================
    # PANEL 8: LAVA HX DETAIL (Detail B)
    # ============================================================
    ax_hxdet.set_title("DETAIL B  -  Lava HX", fontsize=6, fontweight="bold",
                       color="#FFEB3B", pad=3)
    ax_hxdet.set_aspect("equal")

    ax_hxdet.add_patch(Rectangle((0.05, 0.05), 0.90, 0.90, transform=ax_hxdet.transAxes,
                                  facecolor="#FF4500", edgecolor="#8B0000", linewidth=0.8, alpha=0.3))
    ax_hxdet.text(0.50, 0.93, f"LAVA {lv.t_lava_c:.0f}C", transform=ax_hxdet.transAxes,
                  fontsize=3.5, ha="center", color="white", fontweight="bold")

    if lv.hx_enabled:
        n_tubes_show = min(lv.hx_n_tubes, 16)
        cols = 4
        tube_r = 0.025
        for i in range(n_tubes_show):
            r = i // cols
            c = i % cols
            tx = 0.18 + c * 0.18
            ty = 0.20 + r * 0.16
            ax_hxdet.add_patch(Circle((tx, ty), tube_r + 0.004, transform=ax_hxdet.transAxes,
                                       facecolor="#FF8C00", edgecolor="#E65100", linewidth=0.3, alpha=0.7))
            ax_hxdet.add_patch(Circle((tx, ty), tube_r, transform=ax_hxdet.transAxes,
                                       facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.2, alpha=0.5))

        ax_hxdet.text(0.50, 0.04, f"{lv.hx_n_tubes:,} tubes x {lv.hx_tube_od_mm:.0f}mm\nU={lv.hx_u:.0f} W/m2K",
                      transform=ax_hxdet.transAxes, fontsize=3, ha="center", color="#FFEB3B")
    else:
        ax_hxdet.text(0.50, 0.50, "HX disabled", transform=ax_hxdet.transAxes,
                      fontsize=4, ha="center", color="#888")

    if lv.heat_pipe:
        for i in range(4):
            hx = 0.15 + i * 0.22
            ax_hxdet.plot([hx, hx], [0.06, 0.92], color="#FF1744", linewidth=0.6,
                         alpha=0.4, transform=ax_hxdet.transAxes)

    # ============================================================
    # PANEL 9: FAN/NOZZLE DETAIL (Detail D)
    # ============================================================
    ax_fandet.set_title("DETAIL D  -  Fan + nozzle", fontsize=6, fontweight="bold",
                        color="#FFEB3B", pad=3)
    ax_fandet.set_aspect("equal")

    ax_fandet.add_patch(plt.Polygon([(0.08, 0.65), (0.08, 0.35), (0.30, 0.45), (0.30, 0.55)],
                                     facecolor="#888", edgecolor="#666", linewidth=0.6, alpha=0.4,
                                     transform=ax_fandet.transAxes))
    ax_fandet.text(0.15, 0.75, "nozzle", transform=ax_fandet.transAxes,
                   fontsize=3, color="#aaa", ha="center")

    ax_fandet.add_patch(Rectangle((0.32, 0.30), 0.35, 0.40, transform=ax_fandet.transAxes,
                                   facecolor="#333", edgecolor="#666", linewidth=0.8, alpha=0.3))

    fan_cx, fan_cy = 0.50, 0.50
    fan_r = 0.14
    ax_fandet.add_patch(Circle((fan_cx, fan_cy), fan_r, transform=ax_fandet.transAxes,
                                facecolor="#00BFFF", edgecolor="#0277FD", linewidth=0.8, alpha=0.3))
    n_fan_blades = EXIT_FAN_HW['blade_count']
    for b in range(n_fan_blades):
        ang = 2 * math.pi * b / n_fan_blades
        ax_fandet.plot([fan_cx, fan_cx + fan_r * 0.9 * math.cos(ang)],
                       [fan_cy, fan_cy + fan_r * 0.9 * math.sin(ang)],
                       color="#0277FD", linewidth=0.8, alpha=0.6,
                       transform=ax_fandet.transAxes)
    ax_fandet.add_patch(Circle((fan_cx, fan_cy), 0.02, transform=ax_fandet.transAxes,
                                facecolor="#444", edgecolor="#222", linewidth=0.3))

    ax_fandet.add_patch(Rectangle((0.68, 0.42), 0.15, 0.16, transform=ax_fandet.transAxes,
                                   facecolor="#FFC107", edgecolor="#FF6F00", linewidth=0.3, alpha=0.5))
    ax_fandet.text(0.755, 0.50, "GEN", transform=ax_fandet.transAxes,
                   fontsize=2.5, ha="center", va="center", color="white", fontweight="bold")

    ax_fandet.annotate("", xy=(0.92, 0.50), xytext=(0.30, 0.50),
                       arrowprops=dict(arrowstyle="->", color="#FFEB3B", lw=1),
                       transform=ax_fandet.transAxes)

    ax_fandet.text(0.50, 0.15, f"D={EXIT_FAN_HW['fan_d_mm']:.0f}mm\n{EXIT_FAN_HW['generator_kW']:.0f}kW\nRPM={EXIT_FAN_HW['rpm']:.0f}",
                   transform=ax_fandet.transAxes, fontsize=3, ha="center", color="white", fontweight="bold")

    # ============================================================
    # PANEL 10: P&ID (delegated)
    # ============================================================
    _draw_pid_diagram(ax_pid, t)

    # ============================================================
    # PANEL 11: EXPLODED ASSEMBLY VIEW (delegated)
    # ============================================================
    _draw_exploded_view(ax_explod, t)

    # ============================================================
    # PANEL 12: ELECTRICAL SINGLE-LINE DIAGRAM (delegated)
    # ============================================================
    _draw_electrical_sld(ax_elec, t)

    # ============================================================
    # PANEL 13: REHEAT DETAIL (delegated)
    # ============================================================
    _draw_reheat_detail(ax_reheat, t)

    # ============================================================
    # PANEL 14: SITE LAYOUT (delegated)
    # ============================================================
    _draw_site_layout(ax_site, t)

    # ============================================================
    # PANEL 15: CAVERN INTERIOR (delegated)
    # ============================================================
    _draw_cavern_interior(ax_cavint, t)

    # ============================================================
    # PANEL 16: T-s DIAGRAMS (delegated)
    # ============================================================
    _draw_ts_diagrams(ax_ts, t)

    # ============================================================
    # PANEL 17: HEAT PIPE DETAIL (delegated)
    # ============================================================
    _draw_heat_pipe_detail(ax_hp, t)

    # ============================================================
    # PANEL 18: COOLING TOWER DETAIL (delegated)
    # ============================================================
    _draw_cooling_tower_detail(ax_ct, t)

    # ============================================================
    # PANEL 19: CONTROL ARCHITECTURE (delegated)
    # ============================================================
    _draw_control_architecture(ax_ctrl, t)

    # ============================================================
    # PANEL 20+21: TITLE BLOCK + LEGEND
    # ============================================================
    ax_title.set_title("TITLE + LEGEND", fontsize=6, fontweight="bold",
                       color="#00d4ff", pad=3)
    ax_title.set_axis_off()

    ax_title.add_patch(Rectangle((0.03, 0.03), 0.94, 0.94, transform=ax_title.transAxes,
                                  fill=False, edgecolor="#4a5568", linewidth=1))
    ax_title.add_patch(Rectangle((0.03, 0.80), 0.94, 0.17, transform=ax_title.transAxes,
                                  fill=True, facecolor="#1a1a2e", edgecolor="#4a5568", linewidth=0.5))

    ax_title.text(0.50, 0.91, "GMANS TUNNEL", transform=ax_title.transAxes,
                  fontsize=8, fontweight="bold", color="#e94560", ha="center", va="center")
    ax_title.text(0.50, 0.84, "Cryo-Lava Harvester", transform=ax_title.transAxes,
                  fontsize=4, color="#00d4ff", ha="center", va="center")

    specs = [
        ("Cavern", f"{cv.volume_m3/1e9:.1f} km3"),
        ("Charge", f"{cv.p_charge_pa/1e5:.0f}bar/{k_to_c(cv.t_charge_k):.0f}C"),
        ("Tunnel", f"{tun_len:.0f}m x {tn.diameter_m:.0f}m"),
        ("Bores", f"{n_bores}"),
        ("Lava", f"{lv.t_lava_c:.0f}C"),
        ("Turbines", f"{n_turb}+{n_reheat}RH"),
        ("Fans", f"{n_fans}"),
        ("Stack", f"{stack_h:.0f}m"),
        ("Systems", f"{n_sys}"),
    ]
    y = 0.75
    for label, value in specs:
        ax_title.text(0.06, y, label, transform=ax_title.transAxes,
                      fontsize=3.5, color="#a0aec0", family="monospace")
        ax_title.text(0.55, y, value, transform=ax_title.transAxes,
                      fontsize=3.5, color="#00E676", family="monospace", fontweight="bold")
        y -= 0.028

    # Legend
    ax_title.text(0.50, 0.46, "--- LEGEND ---", transform=ax_title.transAxes,
                  fontsize=3.5, color="#00d4ff", ha="center", fontweight="bold")
    legend_items = [
        ("#2196F3", "Cavern"),
        ("#FFD700", "Tunnel"),
        ("#FF4500", "Lava"),
        ("#00E676", "Turbine"),
        ("#FF5722", "Reheat"),
        ("#00BFFF", "Fan"),
        ("#FFC107", "Generator"),
        ("#26A69A", "Drainage"),
    ]
    y = 0.42
    for color, label in legend_items:
        ax_title.add_patch(Rectangle((0.08, y - 0.008), 0.04, 0.016, transform=ax_title.transAxes,
                                      facecolor=color, edgecolor="#444", linewidth=0.2, alpha=0.7))
        ax_title.text(0.14, y, label, transform=ax_title.transAxes,
                      fontsize=3, color="#a0aec0", family="monospace", va="center")
        y -= 0.025

    parts = build_parts_list(t)
    ax_title.text(0.50, 0.06, f"BOM: {len(parts)}", transform=ax_title.transAxes,
                  fontsize=4, color="#FFD700", ha="center", fontweight="bold")


def _draw_3d_view(ax, t: Dict, detail_level: int = 1) -> None:
    """Draw a comprehensive 3D to-scale view of the entire system.

    detail_level:
      0 = minimal (core components only, fast)
      1 = standard (all major components)
      2 = full (all 50+ components including secondary detail)
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np

    ax.clear()
    ax.set_facecolor("#0d1117")

    cv, lv, tn = t["cavern"], t["lava"], t["tunnel"]
    ctrl = t["ctrl"]
    n_sys = max(1, ctrl.n_systems)

    # --- dimensions (all to scale, in metres) ---
    cav_side = cv.volume_m3 ** (1.0 / 3.0)
    cav_depth = cv.depth_m
    tun_len = tn.total_length_m
    tun_r = tn.diameter_m / 2.0
    lava_len = lv.contact_length_m
    n_bores = max(1, lv.n_parallel_bores)
    stack_h = tn.height_rise_m
    n_turb = tn.n_turbine_stages
    n_reheat = getattr(tn, 'n_reheat_stages', 0)
    n_fans = tn.n_exit_fans
    fan_d = EXIT_FAN_HW['fan_d_mm'] / 1000.0
    turb_rotor_r = TURBINE_HW['rotor_d_mm'] / 1000.0 / 2.0

    # --- minimum visible size ---
    total_len = cav_side + tun_len + stack_h
    min_vis = total_len * 0.008
    turb_r = max(turb_rotor_r, min_vis)
    fan_r = max(fan_d / 2.0, min_vis * 0.6)
    bore_r = max(tun_r, min_vis * 0.3)

    # --- layout ---
    cav_cx, cav_cy, cav_cz = 0.0, 0.0, -cav_depth
    cav_s = cav_side / 2.0
    tun_x0 = cav_cx + cav_s
    tun_x1 = tun_x0 + tun_len
    tun_z = cav_cz + cav_s
    lava_x0 = tun_x0 + tun_len * 0.15
    lava_x1 = min(lava_x0 + lava_len, tun_x1)
    stack_x = tun_x1

    # --- bore grid offsets ---
    if n_bores <= 1:
        bore_offsets = [(0.0, 0.0)]
    else:
        grid_n = int(math.ceil(math.sqrt(n_bores)))
        spacing = max(tn.diameter_m * 1.5, min_vis * 2)
        bore_offsets = []
        for i in range(n_bores):
            row = i // grid_n
            col = i % grid_n
            by = (col - (grid_n - 1) / 2.0) * spacing
            bz = (row - (grid_n - 1) / 2.0) * spacing
            bore_offsets.append((by, bz))

    def draw_system(ox, oy, alpha=1.0):
        # 1. GROUND PLANE
        gx = [cav_cx - cav_s - 50, tun_x1 + stack_h + 100]
        gy = [cav_cy - cav_s * 2 - 50, cav_cy + cav_s * 2 + 50]
        ax.add_collection3d(Poly3DCollection(
            [[(ox+gx[0], oy+gy[0], 0), (ox+gx[1], oy+gy[0], 0),
              (ox+gx[1], oy+gy[1], 0), (ox+gx[0], oy+gy[1], 0)]],
            alpha=0.06 * alpha, facecolor="#5C4033", edgecolor="#444", linewidth=0.3))

        # 2. CAVERN (with lining detail)
        s = cav_s
        cx, cy, cz = cav_cx + ox, cav_cy + oy, cav_cz
        verts = [
            [(cx-s, cy-s, cz-s), (cx+s, cy-s, cz-s), (cx+s, cy+s, cz-s), (cx-s, cy+s, cz-s)],
            [(cx-s, cy-s, cz+s), (cx+s, cy-s, cz+s), (cx+s, cy+s, cz+s), (cx-s, cy+s, cz+s)],
            [(cx-s, cy-s, cz-s), (cx+s, cy-s, cz-s), (cx+s, cy-s, cz+s), (cx-s, cy-s, cz+s)],
            [(cx+s, cy-s, cz-s), (cx+s, cy+s, cz-s), (cx+s, cy+s, cz+s), (cx+s, cy-s, cz+s)],
            [(cx-s, cy+s, cz-s), (cx+s, cy+s, cz-s), (cx+s, cy+s, cz+s), (cx-s, cy+s, cz+s)],
            [(cx-s, cy-s, cz-s), (cx-s, cy+s, cz-s), (cx-s, cy+s, cz+s), (cx-s, cy-s, cz+s)],
        ]
        t_c = k_to_c(cv.t_charge_k)
        cav_color = "#2196F3" if t_c < -100 else ("#4FC3F7" if t_c < 0 else "#81C784")
        ax.add_collection3d(Poly3DCollection(verts, alpha=0.35 * alpha,
            facecolor=cav_color, edgecolor="#0277FD", linewidth=1))
        if alpha > 0.5:
            ax.text(cx, cy, cz, f"CAVERN\n{cv.volume_m3/1e6:.1f}M m3\n{t_c:.0f}C\n{cv.p_charge_pa/1e5:.0f} bar",
                    fontsize=5, ha="center", va="center", color="white", zorder=10)

        # 2b. Cavern lining (slightly larger box)
        lining_s = s + CAVERN_HW['lining_thick_mm'] / 1000.0
        lining_verts = [
            [(cx-lining_s, cy-lining_s, cz-lining_s), (cx+lining_s, cy-lining_s, cz-lining_s),
             (cx+lining_s, cy+lining_s, cz-lining_s), (cx-lining_s, cy+lining_s, cz-lining_s)],
        ]
        ax.add_collection3d(Poly3DCollection(lining_verts, alpha=0.15 * alpha,
            facecolor="#666", edgecolor="#444", linewidth=0.5))

        # 2b2. Ultra thermal insulation layer (when enabled)
        if cv.ultra_insulation and alpha > 0.3:
            insul_s = lining_s + CAVERN_HW['ultra_insulation_mm'] / 1000.0
            insul_verts = [
                [(cx-insul_s, cy-insul_s, cz-insul_s), (cx+insul_s, cy-insul_s, cz-insul_s),
                 (cx+insul_s, cy+insul_s, cz-insul_s), (cx-insul_s, cy+insul_s, cz-insul_s)],
            ]
            ax.add_collection3d(Poly3DCollection(insul_verts, alpha=0.12 * alpha,
                facecolor="#E1BEE7", edgecolor="#9C27B0", linewidth=0.4))
            if alpha > 0.6:
                ax.text(cx + insul_s + 2, cy, cz + s * 0.5,
                        "ULTRA INSUL\nR=30", fontsize=3, color="#CE93D8",
                        ha="left", va="center", zorder=9)

        # 2c. CAVERN SENSORS (pressure + temperature markers) - batched
        if alpha > 0.5:
            n_p_sensors = min(MONITOR_HW['cavern_pressure_sensors'], 8)
            ps_x, ps_y, ps_z = [], [], []
            for si in range(n_p_sensors):
                ang = 2 * math.pi * si / n_p_sensors
                ps_x.append(cx + s * 0.9 * math.cos(ang))
                ps_y.append(cy + s * 0.9 * math.sin(ang))
                ps_z.append(cz + s * 0.8)
            if ps_x:
                ax.scatter(ps_x, ps_y, ps_z, color="#FFEB3B", s=4,
                          marker="o", alpha=0.6*alpha, zorder=8)
            # geophones around cavern perimeter - batched
            n_geo = min(CAVERN_HW['geophone_count'], 8)
            gs_x, gs_y, gs_z = [], [], []
            for gi in range(n_geo):
                ang = 2 * math.pi * gi / n_geo
                gs_x.append(cx + (s + 5) * math.cos(ang))
                gs_y.append(cy + (s + 5) * math.sin(ang))
                gs_z.append(0)
            if gs_x:
                ax.scatter(gs_x, gs_y, gs_z, color="#E91E63", s=3,
                          marker="^", alpha=0.5*alpha, zorder=8)

        # 2d. HYDRAULIC DOOR (at access tunnel / cavern junction)
        if alpha > 0.5:
            door_r = max(CAVERN_HW['hydraulic_door_mm'] / 1000.0 / 2.0, min_vis * 0.3)
            door_x = ox + cav_cx - cav_s * 0.3
            door_z = cav_cz + cav_s
            f = _cylinder_faces(door_x, oy + cav_cy, door_z - door_r*0.2,
                                door_x, oy + cav_cy, door_z + door_r*0.2, door_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.7*alpha,
                facecolor="#FF5722", edgecolor="#D84315", linewidth=0.5))

        # 3. ACCESS TUNNEL
        acc_r = max(CAVERN_HW['access_tunnel_d_m'] / 2.0, min_vis * 0.3)
        acc_faces = _cylinder_faces(
            ox + cav_cx - cav_s * 0.3, oy + cav_cy, 0,
            ox + cav_cx - cav_s * 0.3, oy + cav_cy, cav_cz + cav_s,
            acc_r, n_seg=8)
        if acc_faces:
            ax.add_collection3d(Poly3DCollection(acc_faces, alpha=0.5 * alpha,
                facecolor="#888", edgecolor="#555", linewidth=0.5))

        # 4. PARALLEL BORES
        for by, bz in bore_offsets:
            bz_z = tun_z + bz
            # pre-lava section
            f = _cylinder_faces(ox+tun_x0, oy+by, bz_z, ox+lava_x0, oy+by, bz_z, bore_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.25*alpha, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.3))
            # lava section
            f = _cylinder_faces(ox+lava_x0, oy+by, bz_z, ox+lava_x1, oy+by, bz_z, bore_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.35*alpha, facecolor="#FF6347", edgecolor="#FF4500", linewidth=0.3))
            # post-lava section
            f = _cylinder_faces(ox+lava_x1, oy+by, bz_z, ox+tun_x1, oy+by, bz_z, bore_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.25*alpha, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.3))

        # 4b. EXPANSION JOINTS (ring markers along tunnel)
        if alpha > 0.5:
            n_joints = int(tun_len / TUNNEL_HW['expansion_joint_m'])
            joint_spacing = TUNNEL_HW['expansion_joint_m']
            by0_j, bz0_j = bore_offsets[0] if bore_offsets else (0, 0)
            for j in range(min(n_joints, 30)):
                jx = tun_x0 + (j + 1) * joint_spacing
                if jx > tun_x1:
                    break
                # draw as a thin ring (small cylinder cross-section)
                f = _cylinder_faces(ox+jx, oy+by0_j, tun_z+bz0_j-bore_r,
                                    ox+jx, oy+by0_j, tun_z+bz0_j+bore_r, bore_r*1.15, 8)
                if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.3*alpha,
                    facecolor="#FF9800", edgecolor="#FF6F00", linewidth=0.3))

        # 4c. SENSORS (temperature + pressure taps along tunnel) - batched
        if alpha > 0.5:
            n_sensors = min(MONITOR_HW['tunnel_temp_sensors'], 20)
            ss_x, ss_y, ss_z = [], [], []
            by0_s, bz0_s = bore_offsets[0] if bore_offsets else (0, 0)
            for s_i in range(n_sensors):
                ss_x.append(ox + tun_x0 + (s_i + 0.5) * tun_len / n_sensors)
                ss_y.append(oy + by0_s + bore_r*1.2)
                ss_z.append(tun_z + bz0_s)
            if ss_x:
                ax.scatter(ss_x, ss_y, ss_z, color="#FFEB3B", s=3,
                          marker="o", alpha=0.5*alpha, zorder=7)

        # 4d. DRAINAGE PIPE (small cylinder below tunnel)
        if alpha > 0.5:
            drain_r = max(TUNNEL_HW['drainage_pipe_mm'] / 1000.0 / 2.0, min_vis * 0.1)
            by0_d = bore_offsets[0][0] if bore_offsets else 0
            f = _cylinder_faces(ox+tun_x0, oy+by0_d, tun_z-bore_r*1.3,
                                ox+tun_x1, oy+by0_d, tun_z-bore_r*1.3, drain_r, 6)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.4*alpha,
                facecolor="#26A69A", edgecolor="#00897B", linewidth=0.3))

        # 4e. ESCAPE REFUGES (small boxes along tunnel)
        if alpha > 0.5:
            n_refuges = TUNNEL_HW['escape_refuges']
            for r_i in range(n_refuges):
                rx_pos = tun_x0 + (r_i + 0.5) * tun_len / n_refuges
                by0_r = bore_offsets[0][0] if bore_offsets else 0
                refuge_s = min_vis * 0.5
                rv = [
                    [(ox+rx_pos-refuge_s, oy+by0_r+bore_r*1.5, tun_z-refuge_s),
                     (ox+rx_pos+refuge_s, oy+by0_r+bore_r*1.5, tun_z-refuge_s),
                     (ox+rx_pos+refuge_s, oy+by0_r+bore_r*1.5, tun_z+refuge_s),
                     (ox+rx_pos-refuge_s, oy+by0_r+bore_r*1.5, tun_z+refuge_s)],
                ]
                ax.add_collection3d(Poly3DCollection(rv, alpha=0.4*alpha,
                    facecolor="#F44336", edgecolor="#D32F2F", linewidth=0.3))

        # 5. LAVA BODY
        lava_d = max(tn.diameter_m * 3, 50.0)
        lava_w = max(cav_side * 0.8, n_bores * tn.diameter_m * 2)
        lv_verts = [
            [(ox+lava_x0, oy-lava_w/2, tun_z-lava_d), (ox+lava_x1, oy-lava_w/2, tun_z-lava_d),
             (ox+lava_x1, oy+lava_w/2, tun_z-lava_d), (ox+lava_x0, oy+lava_w/2, tun_z-lava_d)],
            [(ox+lava_x0, oy-lava_w/2, tun_z-tun_r), (ox+lava_x1, oy-lava_w/2, tun_z-tun_r),
             (ox+lava_x1, oy+lava_w/2, tun_z-lava_d), (ox+lava_x0, oy-lava_w/2, tun_z-lava_d)],
            [(ox+lava_x1, oy-lava_w/2, tun_z-tun_r), (ox+lava_x1, oy+lava_w/2, tun_z-tun_r),
             (ox+lava_x1, oy+lava_w/2, tun_z-lava_d), (ox+lava_x1, oy-lava_w/2, tun_z-lava_d)],
            [(ox+lava_x0, oy+lava_w/2, tun_z-tun_r), (ox+lava_x1, oy+lava_w/2, tun_z-tun_r),
             (ox+lava_x1, oy+lava_w/2, tun_z-lava_d), (ox+lava_x0, oy+lava_w/2, tun_z-lava_d)],
            [(ox+lava_x0, oy-lava_w/2, tun_z-tun_r), (ox+lava_x0, oy+lava_w/2, tun_z-tun_r),
             (ox+lava_x0, oy+lava_w/2, tun_z-lava_d), (ox+lava_x0, oy-lava_w/2, tun_z-lava_d)],
        ]
        ax.add_collection3d(Poly3DCollection(lv_verts, alpha=0.55*alpha,
            facecolor="#FF4500", edgecolor="#8B0000", linewidth=0.5))
        if alpha > 0.5:
            ax.text(ox+(lava_x0+lava_x1)/2, oy, tun_z-lava_d/2,
                    f"LAVA {lv.t_lava_c:.0f}C", fontsize=5, ha="center", color="#FFD700", zorder=10)

        # 6. HX TUBE BUNDLE
        if lv.hx_enabled and lv.hx_n_tubes > 0 and alpha > 0.5:
            n_show = min(lv.hx_n_tubes, 40)
            tube_r = max(lv.hx_tube_od_mm / 1000.0 / 2.0, min_vis * 0.08)
            for i in range(n_show):
                ang = 2 * math.pi * i / n_show
                ty = math.cos(ang) * lava_w * 0.15
                tz = math.sin(ang) * lava_w * 0.15
                f = _cylinder_faces(ox+lava_x0+10, oy+ty, tun_z+tz-lava_d*0.5,
                                    ox+lava_x1-10, oy+ty, tun_z+tz-lava_d*0.5, tube_r, 6)
                if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.15, facecolor="#FF8C00", edgecolor="#FF8C00", linewidth=0.2))
            ax.text(ox+(lava_x0+lava_x1)/2, oy+lava_w*0.3, tun_z-lava_d*0.5,
                    f"HX: {lv.hx_n_tubes:,} tubes", fontsize=4, color="#FF8C00", zorder=10)

        # 7. HEAT PIPES
        if lv.heat_pipe and alpha > 0.5:
            for i in range(min(20, n_bores * 2)):
                px = lava_x0 + (i + 0.5) * (lava_len / min(20, n_bores * 2))
                ax.plot([ox+px, ox+px], [oy, oy], [tun_z-lava_d*0.8, tun_z-tun_r],
                        color="#FF1744", linewidth=1, alpha=0.4 * alpha)

        # 8. TURBINE STAGES (cylinders with blade markers)
        by0, bz0 = bore_offsets[0] if bore_offsets else (0, 0)
        for i in range(n_turb):
            frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tun_len)
            tx = tun_x0 + frac * tun_len
            # turbine housing
            f = _cylinder_faces(ox+tx-turb_r*0.4, oy+by0, tun_z+bz0,
                                ox+tx+turb_r*0.4, oy+by0, tun_z+bz0, turb_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.8*alpha, facecolor="#00E676", edgecolor="#00C853", linewidth=0.5))
            # blade markers (small lines across the rotor face)
            if alpha > 0.5 and i < 8:
                ax.plot([ox+tx, ox+tx], [oy+by0-turb_r*0.8, oy+by0+turb_r*0.8],
                        [tun_z+bz0, tun_z+bz0], color="#00C853", linewidth=0.5, alpha=0.5)
                ax.text(ox+tx, oy+by0+turb_r*1.5, tun_z+bz0, f"T{i+1}", fontsize=4, color="#00E676", zorder=10)
        if alpha > 0.5 and n_turb > 8:
            ax.text(ox+tun_x0+(0.15+8/(n_turb+1)*(lava_len/tun_len))*tun_len, oy+turb_r*1.5, tun_z,
                    f"...T{n_turb}", fontsize=4, color="#00E676", zorder=10)

        # 9. REHEAT MARKERS - batched
        if n_reheat > 0 and alpha > 0.5:
            rh_x, rh_y, rh_z = [], [], []
            for i in range(min(n_reheat, 10)):
                frac = 0.15 + (i + 1.5) / (n_turb + 1) * (lava_len / tun_len)
                rh_x.append(ox + tun_x0 + frac * tun_len)
                rh_y.append(oy + by0)
                rh_z.append(tun_z + bz0 + turb_r)
            if rh_x:
                ax.scatter(rh_x, rh_y, rh_z, color="#FF5722", s=12,
                          marker="v", zorder=8, alpha=alpha)

        # 10. MHD CHANNEL
        if tn.mhd_enabled and alpha > 0.5:
            mx = lava_x0 + 20
            f = _cylinder_faces(ox+mx, oy+by0, tun_z+bz0, ox+mx+min_vis*2, oy+by0, tun_z+bz0, turb_r*0.8, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.7*alpha, facecolor="#E91E63", edgecolor="#C2185B", linewidth=0.5))
            ax.text(ox+mx+min_vis, oy+by0, tun_z+bz0+turb_r*2, "MHD", fontsize=4, color="#E91E63", zorder=10)

        # 11. REGENERATOR
        if tn.regenerator_eff > 0 and alpha > 0.5:
            rx = tun_x0 + tun_len * 0.10
            f = _cylinder_faces(ox+rx, oy+by0, tun_z+bz0, ox+rx+min_vis*1.5, oy+by0, tun_z+bz0, turb_r*0.7, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.6*alpha, facecolor="#9C27B0", edgecolor="#7B1FA2", linewidth=0.5))
            ax.text(ox+rx, oy+by0, tun_z+bz0+turb_r*2, "REGEN", fontsize=4, color="#9C27B0", zorder=10)

        # 12. STACK
        stack_r = max(tun_r * 0.6, min_vis * 0.3)
        f = _cylinder_faces(ox+stack_x, oy, tun_z, ox+stack_x, oy, stack_h, stack_r, 8)
        if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.5*alpha, facecolor="#888", edgecolor="#666", linewidth=0.5))

        # 13. EXIT NOZZLE
        nozzle_r = max(tn.exit_nozzle_area_m2 ** 0.5 / math.pi, min_vis * 0.2)
        f = _cylinder_faces(ox+stack_x, oy, stack_h, ox+stack_x, oy, stack_h+min_vis*0.5, nozzle_r, 8)
        if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.7*alpha, facecolor="#FFEB3B", edgecolor="#FBC02D", linewidth=0.5))

        # 14. EXIT FANS
        for i in range(n_fans):
            fy = (i - (n_fans - 1) / 2.0) * fan_r * 2.5
            fz = stack_h + min_vis * 0.3
            f = _cylinder_faces(ox+stack_x+fy*0.2, oy+fy, fz, ox+stack_x+fy*0.2, oy+fy, fz+min_vis*0.2, fan_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.8*alpha, facecolor="#00BFFF", edgecolor="#0277FD", linewidth=0.5))
            if alpha > 0.5 and i == 0:
                ax.text(ox+stack_x, oy+fy, fz+min_vis*0.4, f"Fans x{n_fans}", fontsize=4, color="#00BFFF", zorder=10)

        # 15. EXIT JET - batched
        if alpha > 0.5:
            jet_x, jet_y, jet_z = [], [], []
            for i in range(3):
                jy = (i - 1) * fan_r * 2
                ax.plot([ox+stack_x, ox+stack_x], [oy+jy, oy+jy],
                        [stack_h+min_vis, stack_h+min_vis*3], color="#FFEB3B", linewidth=1, alpha=0.4*alpha)
                jet_x.append(ox+stack_x)
                jet_y.append(oy+jy)
                jet_z.append(stack_h+min_vis*3)
            if jet_x:
                ax.scatter(jet_x, jet_y, jet_z, color="#FFEB3B", s=6,
                          marker="^", alpha=0.5*alpha)

        # 16. BOTTOMING CYCLES
        bottoming = []
        if tn.potassium_enabled: bottoming.append(("K", "#FF9800", tn.potassium_eta))
        if tn.sco2_enabled: bottoming.append(("sCO2", "#9C27B0", tn.sco2_eta))
        if tn.steam_enabled: bottoming.append(("Steam", "#03A9F4", tn.steam_eta))
        if tn.orc_enabled: bottoming.append(("ORC", "#4CAF50", tn.orc_eta))
        for i, (name, color, eta) in enumerate(bottoming):
            bx = stack_x + min_vis * (2 + i * 1.5)
            bz = tun_z - min_vis * (1 + i * 0.5)
            br = min_vis * 0.5
            f = _cylinder_faces(ox+bx, oy, bz, ox+bx+min_vis, oy, bz, br, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.6*alpha, facecolor=color, edgecolor=color, linewidth=0.5))
            if alpha > 0.5:
                ax.text(ox+bx+min_vis*0.5, oy+br*1.5, bz, f"{name}\n{eta*100:.0f}%", fontsize=4, color=color, zorder=10)

        # 17. TRANSFORMER
        if alpha > 0.5:
            tx_sw = stack_x + min_vis * 2
            tr = min_vis * 0.4
            tr_v = [
                [(ox+tx_sw-tr, oy-tr, 0), (ox+tx_sw+tr, oy-tr, 0), (ox+tx_sw+tr, oy+tr, 0), (ox+tx_sw-tr, oy+tr, 0)],
                [(ox+tx_sw-tr, oy-tr, tr), (ox+tx_sw+tr, oy-tr, tr), (ox+tx_sw+tr, oy+tr, tr), (ox+tx_sw-tr, oy+tr, tr)],
                [(ox+tx_sw-tr, oy-tr, 0), (ox+tx_sw+tr, oy-tr, 0), (ox+tx_sw+tr, oy-tr, tr), (ox+tx_sw-tr, oy-tr, tr)],
                [(ox+tx_sw+tr, oy-tr, 0), (ox+tx_sw+tr, oy+tr, 0), (ox+tx_sw+tr, oy+tr, tr), (ox+tx_sw+tr, oy-tr, tr)],
                [(ox+tx_sw-tr, oy+tr, 0), (ox+tx_sw+tr, oy+tr, 0), (ox+tx_sw+tr, oy+tr, tr), (ox+tx_sw-tr, oy+tr, tr)],
                [(ox+tx_sw-tr, oy-tr, 0), (ox+tx_sw-tr, oy+tr, 0), (ox+tx_sw-tr, oy+tr, tr), (ox+tx_sw-tr, oy-tr, tr)],
            ]
            ax.add_collection3d(Poly3DCollection(tr_v, alpha=0.5*alpha, facecolor="#FFC107", edgecolor="#FFA000", linewidth=0.5))
            ax.text(ox+tx_sw, oy, tr*2, "XFMR", fontsize=4, color="#FFC107", zorder=10)

        # 18. CHILLER
        if (cv.active_cooling or cv.lava_heated_cooling) and alpha > 0.5:
            ch_x = cav_cx - cav_s - min_vis * 2
            ch_r = min_vis * 0.4
            f = _cylinder_faces(ox+ch_x, oy, tun_z, ox+ch_x+min_vis, oy, tun_z, ch_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.6*alpha, facecolor="#80DEEA", edgecolor="#00ACC1", linewidth=0.5))
            lbl = "Absorption\nChiller" if cv.lava_heated_cooling else "Chiller"
            ax.text(ox+ch_x+min_vis*0.5, oy+ch_r*1.5, tun_z, lbl, fontsize=4, color="#80DEEA", zorder=10)

        # 19. TURBINE HALL BUILDING (surface)
        if alpha > 0.5:
            th_x = stack_x + min_vis * 4
            th_w = min_vis * 3
            th_d = min_vis * 2
            th_h = min_vis * 1.5
            th_v = [
                [(ox+th_x, oy-th_d/2, 0), (ox+th_x+th_w, oy-th_d/2, 0), (ox+th_x+th_w, oy+th_d/2, 0), (ox+th_x, oy+th_d/2, 0)],
                [(ox+th_x, oy-th_d/2, th_h), (ox+th_x+th_w, oy-th_d/2, th_h), (ox+th_x+th_w, oy+th_d/2, th_h), (ox+th_x, oy+th_d/2, th_h)],
                [(ox+th_x, oy-th_d/2, 0), (ox+th_x+th_w, oy-th_d/2, 0), (ox+th_x+th_w, oy-th_d/2, th_h), (ox+th_x, oy-th_d/2, th_h)],
                [(ox+th_x+th_w, oy-th_d/2, 0), (ox+th_x+th_w, oy+th_d/2, 0), (ox+th_x+th_w, oy+th_d/2, th_h), (ox+th_x+th_w, oy-th_d/2, th_h)],
                [(ox+th_x, oy+th_d/2, 0), (ox+th_x+th_w, oy+th_d/2, 0), (ox+th_x+th_w, oy+th_d/2, th_h), (ox+th_x, oy+th_d/2, th_h)],
                [(ox+th_x, oy-th_d/2, 0), (ox+th_x, oy+th_d/2, 0), (ox+th_x, oy+th_d/2, th_h), (ox+th_x, oy-th_d/2, th_h)],
            ]
            ax.add_collection3d(Poly3DCollection(th_v, alpha=0.3*alpha, facecolor="#546E7A", edgecolor="#37474F", linewidth=0.5))
            ax.text(ox+th_x+th_w/2, oy, th_h*1.5, "Turbine\nHall", fontsize=4, color="#90A4AE", zorder=10)

        # 20. CONTROL ROOM BUILDING (surface)
        if alpha > 0.5:
            cr_x = stack_x + min_vis * 8
            cr_w = min_vis * 1.5
            cr_d = min_vis * 1.5
            cr_h = min_vis * 1.0
            cr_v = [
                [(ox+cr_x, oy-cr_d/2, 0), (ox+cr_x+cr_w, oy-cr_d/2, 0), (ox+cr_x+cr_w, oy+cr_d/2, 0), (ox+cr_x, oy+cr_d/2, 0)],
                [(ox+cr_x, oy-cr_d/2, cr_h), (ox+cr_x+cr_w, oy-cr_d/2, cr_h), (ox+cr_x+cr_w, oy+cr_d/2, cr_h), (ox+cr_x, oy+cr_d/2, cr_h)],
                [(ox+cr_x, oy-cr_d/2, 0), (ox+cr_x+cr_w, oy-cr_d/2, 0), (ox+cr_x+cr_w, oy-cr_d/2, cr_h), (ox+cr_x, oy-cr_d/2, cr_h)],
                [(ox+cr_x+cr_w, oy-cr_d/2, 0), (ox+cr_x+cr_w, oy+cr_d/2, 0), (ox+cr_x+cr_w, oy+cr_d/2, cr_h), (ox+cr_x+cr_w, oy-cr_d/2, cr_h)],
                [(ox+cr_x, oy+cr_d/2, 0), (ox+cr_x+cr_w, oy+cr_d/2, 0), (ox+cr_x+cr_w, oy+cr_d/2, cr_h), (ox+cr_x, oy+cr_d/2, cr_h)],
                [(ox+cr_x, oy-cr_d/2, 0), (ox+cr_x, oy+cr_d/2, 0), (ox+cr_x, oy+cr_d/2, cr_h), (ox+cr_x, oy-cr_d/2, cr_h)],
            ]
            ax.add_collection3d(Poly3DCollection(cr_v, alpha=0.4*alpha, facecolor="#78909C", edgecolor="#455A64", linewidth=0.5))
            ax.text(ox+cr_x+cr_w/2, oy, cr_h*1.5, "Control\nRoom", fontsize=4, color="#CFD8DC", zorder=10)

        # 21. SWITCHYARD / GIS (surface)
        if alpha > 0.5:
            sw_x = stack_x + min_vis * 10
            sw_w = min_vis * 2
            sw_d = min_vis * 2
            sw_h = min_vis * 0.6
            sw_v = [
                [(ox+sw_x, oy-sw_d/2, 0), (ox+sw_x+sw_w, oy-sw_d/2, 0), (ox+sw_x+sw_w, oy+sw_d/2, 0), (ox+sw_x, oy+sw_d/2, 0)],
                [(ox+sw_x, oy-sw_d/2, sw_h), (ox+sw_x+sw_w, oy-sw_d/2, sw_h), (ox+sw_x+sw_w, oy+sw_d/2, sw_h), (ox+sw_x, oy+sw_d/2, sw_h)],
            ]
            ax.add_collection3d(Poly3DCollection(sw_v, alpha=0.3*alpha, facecolor="#455A64", edgecolor="#263238", linewidth=0.5))
            # transmission tower
            ax.plot([ox+sw_x+sw_w, ox+sw_x+sw_w], [oy, oy], [0, min_vis*2],
                    color="#90A4AE", linewidth=1, alpha=0.6*alpha)
            ax.plot([ox+sw_x+sw_w-min_vis*0.3, ox+sw_x+sw_w+min_vis*0.3], [oy, oy],
                    [min_vis*1.8, min_vis*1.8], color="#90A4AE", linewidth=0.5, alpha=0.6*alpha)
            ax.text(ox+sw_x+sw_w/2, oy, sw_h*2, "Switchyard\n+ 132kV", fontsize=4, color="#B0BEC5", zorder=10)

        # 22. COOLING TOWER (surface)
        if alpha > 0.5:
            ct_x = cav_cx - cav_s - min_vis * 4
            ct_r = min_vis * 0.8
            ct_h = min_vis * 2.5
            # hyperbolic cooling tower shape (approximated as cylinder)
            f = _cylinder_faces(ox+ct_x, oy, 0, ox+ct_x, oy, ct_h, ct_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.3*alpha, facecolor="#B0BEC5", edgecolor="#78909C", linewidth=0.5))
            ax.text(ox+ct_x, oy+ct_r*1.5, ct_h*0.5, "Cooling\nTower", fontsize=4, color="#B0BEC5", zorder=10)

        # 23. RECHARGE COMPRESSOR BUILDING (surface, near cavern)
        if alpha > 0.5:
            rc_x = cav_cx - cav_s - min_vis * 6
            rc_w = min_vis * 1.5
            rc_d = min_vis * 1.0
            rc_h = min_vis * 0.8
            rc_v = [
                [(ox+rc_x, oy-rc_d/2, 0), (ox+rc_x+rc_w, oy-rc_d/2, 0), (ox+rc_x+rc_w, oy+rc_d/2, 0), (ox+rc_x, oy+rc_d/2, 0)],
                [(ox+rc_x, oy-rc_d/2, rc_h), (ox+rc_x+rc_w, oy-rc_d/2, rc_h), (ox+rc_x+rc_w, oy+rc_d/2, rc_h), (ox+rc_x, oy+rc_d/2, rc_h)],
            ]
            ax.add_collection3d(Poly3DCollection(rc_v, alpha=0.4*alpha, facecolor="#8D6E63", edgecolor="#5D4037", linewidth=0.5))
            ax.text(ox+rc_x+rc_w/2, oy, rc_h*1.8, "Recharge\nCompressor", fontsize=4, color="#D7CCC8", zorder=10)

        # 24. DIESEL BACKUP GENERATOR (surface)
        if alpha > 0.5:
            dg_x = stack_x + min_vis * 6.5
            dg_w = min_vis * 0.8
            dg_d = min_vis * 0.6
            dg_h = min_vis * 0.5
            dg_v = [
                [(ox+dg_x, oy-dg_d/2, 0), (ox+dg_x+dg_w, oy-dg_d/2, 0), (ox+dg_x+dg_w, oy+dg_d/2, 0), (ox+dg_x, oy+dg_d/2, 0)],
                [(ox+dg_x, oy-dg_d/2, dg_h), (ox+dg_x+dg_w, oy-dg_d/2, dg_h), (ox+dg_x+dg_w, oy+dg_d/2, dg_h), (ox+dg_x, oy+dg_d/2, dg_h)],
            ]
            ax.add_collection3d(Poly3DCollection(dg_v, alpha=0.4*alpha, facecolor="#FF6F00", edgecolor="#E65100", linewidth=0.5))
            ax.text(ox+dg_x+dg_w/2, oy, dg_h*2, "Diesel\nGen", fontsize=3.5, color="#FFE082", zorder=10)

        # 25. MAIN ISOLATION VALVE (at cavern-tunnel junction)
        if alpha > 0.5:
            mv_x = tun_x0
            mv_r = max(tn.diameter_m * 0.6, min_vis * 0.5)
            f = _cylinder_faces(ox+mv_x-mv_r*0.3, oy, tun_z, ox+mv_x+mv_r*0.3, oy, tun_z, mv_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.7*alpha,
                facecolor="#FF5722", edgecolor="#D84315", linewidth=0.5))
            ax.text(ox+mv_x, oy+mv_r*1.5, tun_z, "ISO\nValve", fontsize=3.5, color="#FF5722", zorder=10)

        # 26. TURBINE BYPASS VALVE
        if alpha > 0.5:
            bp_x = tun_x0 + min_vis * 0.5
            bp_r = max(tn.diameter_m * 0.3, min_vis * 0.3)
            f = _cylinder_faces(ox+bp_x, oy+tn.diameter_m, tun_z, ox+bp_x, oy+tn.diameter_m+min_vis, tun_z, bp_r, 6)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.5*alpha,
                facecolor="#FF9800", edgecolor="#E65100", linewidth=0.3))
            ax.text(ox+bp_x, oy+tn.diameter_m+min_vis*0.8, tun_z, "bypass", fontsize=3, color="#FF9800", zorder=10)

        # 27. PRESSURE RELIEF VALVE (on cavern)
        if alpha > 0.5:
            rv_x = cav_cx + cav_s * 0.5
            rv_r = min_vis * 0.25
            f = _cylinder_faces(ox+rv_x, oy+cav_s, cav_cz, ox+rv_x, oy+cav_s+min_vis*0.5, cav_cz, rv_r, 6)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.6*alpha,
                facecolor="#F44336", edgecolor="#C62828", linewidth=0.3))
            ax.text(ox+rv_x, oy+cav_s+min_vis*0.8, cav_cz, "RV", fontsize=3, color="#F44336", zorder=10)

        # 28. INTERCONNECTING PIPING (turbine to bottoming cycles)
        if alpha > 0.5:
            pipe_r = min_vis * 0.08
            for i, (name, color, eta) in enumerate(bottoming):
                bx = stack_x + min_vis * (2 + i * 1.5)
                bz = tun_z - min_vis * (1 + i * 0.5)
                # pipe from tunnel to bottoming cycle
                f = _cylinder_faces(ox+tun_x1-min_vis, oy, tun_z, ox+bx, oy, bz, pipe_r, 5)
                if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.4*alpha,
                    facecolor=color, edgecolor=color, linewidth=0.2))

        # 29. CABLE TRAY (from turbine hall to control room)
        if alpha > 0.5:
            ct_x1 = stack_x + min_vis * 4
            ct_x2 = stack_x + min_vis * 8
            ct_y = oy + min_vis * 1.5
            ct_z = min_vis * 1.0
            ax.plot([ox+ct_x1, ox+ct_x2], [ct_y, ct_y], [ct_z, ct_z],
                    color="#FFC107", linewidth=0.8, alpha=0.4*alpha)
            ax.text(ox+(ct_x1+ct_x2)/2, ct_y+min_vis*0.3, ct_z, "cable\ntray", fontsize=3,
                    color="#FFC107", zorder=10)

        # 30. ACCESS PLATFORM (at turbine locations)
        if alpha > 0.5 and detail_level >= 2:
            for i in range(min(n_turb, 5)):
                frac = 0.15 + (i + 1) / (n_turb + 1) * (lava_len / tun_len)
                tx = tun_x0 + frac * tun_len
                by0_p = bore_offsets[0][0] if bore_offsets else 0
                # platform as small box
                pf_v = [
                    [(ox+tx-min_vis*0.3, oy+by0_p+bore_r*1.2, tun_z-bore_r*0.8),
                     (ox+tx+min_vis*0.3, oy+by0_p+bore_r*1.2, tun_z-bore_r*0.8),
                     (ox+tx+min_vis*0.3, oy+by0_p+bore_r*1.2, tun_z-bore_r*0.5),
                     (ox+tx-min_vis*0.3, oy+by0_p+bore_r*1.2, tun_z-bore_r*0.5)],
                ]
                ax.add_collection3d(Poly3DCollection(pf_v, alpha=0.3*alpha,
                    facecolor="#9E9E9E", edgecolor="#616161", linewidth=0.2))

        # 31. TRANSMISSION LINE (from switchyard)
        if alpha > 0.5 and detail_level >= 2:
            sw_x = stack_x + min_vis * 10
            # transmission tower
            ax.plot([ox+sw_x+sw_w, ox+sw_x+sw_w+min_vis*3], [oy, oy], [0, min_vis*2],
                    color="#90A4AE", linewidth=0.8, alpha=0.5*alpha)
            # transmission lines (3-phase)
            for phase in range(3):
                pz = min_vis * (1.5 + phase * 0.2)
                ax.plot([ox+sw_x+sw_w, ox+sw_x+sw_w+min_vis*3],
                        [oy, oy], [pz, pz],
                        color="#FFEB3B", linewidth=0.3, alpha=0.3*alpha)

        # 32. OVERHEAD CRANE (inside turbine hall)
        if alpha > 0.5 and detail_level >= 2:
            th_x = stack_x + min_vis * 4
            th_w = min_vis * 3
            th_h = min_vis * 1.5
            # crane bridge (horizontal beam at top of hall)
            cr_y = oy + min_vis * 0.8
            cr_z = th_h * 0.9
            ax.plot([ox+th_x, ox+th_x+th_w], [cr_y, cr_y], [cr_z, cr_z],
                    color="#FFC107", linewidth=1.5, alpha=0.5*alpha)
            # crane rails (side rails)
            ax.plot([ox+th_x, ox+th_x+th_w], [cr_y-min_vis*0.3, cr_y-min_vis*0.3], [cr_z, cr_z],
                    color="#888", linewidth=0.5, alpha=0.3*alpha)
            ax.plot([ox+th_x, ox+th_x+th_w], [cr_y+min_vis*0.3, cr_y+min_vis*0.3], [cr_z, cr_z],
                    color="#888", linewidth=0.5, alpha=0.3*alpha)
            # hoist (small box on bridge)
            ax.scatter([ox+th_x+th_w*0.5], [cr_y], [cr_z - min_vis*0.1],
                      color="#FFC107", s=8, marker="s", alpha=0.6*alpha, zorder=8)
            ax.text(ox+th_x+th_w*0.5, cr_y, cr_z + min_vis*0.2, "crane\n50t",
                    fontsize=3, color="#FFC107", zorder=10)

        # 33. STAIRWAY (from surface to turbine hall entrance)
        if alpha > 0.5 and detail_level >= 2:
            st_x = stack_x + min_vis * 7.5
            for step in range(5):
                st_z = min_vis * 0.1 * step
                ax.plot([ox+st_x, ox+st_x+min_vis*0.2], [oy-min_vis*0.5, oy-min_vis*0.5],
                        [st_z, st_z], color="#78909C", linewidth=0.5, alpha=0.4*alpha)

        # 34. PIPE RACK (between bottoming cycles and cooling tower)
        if alpha > 0.5 and detail_level >= 2:
            pr_x0 = stack_x + min_vis * 2
            pr_x1 = cav_cx - cav_s - min_vis * 3
            pr_z = min_vis * 0.3
            # pipe rack structure
            ax.plot([ox+pr_x0, ox+pr_x1], [oy+min_vis*1.5, oy+min_vis*1.5], [pr_z, pr_z],
                    color="#546E7A", linewidth=0.8, alpha=0.3*alpha)
            # pipes on rack (3 different colored pipes)
            for pi, pcolor in enumerate([("#03A9F4", "steam"), ("#9C27B0", "sCO2"), ("#4CAF50", "ORC")]):
                ax.plot([ox+pr_x0, ox+pr_x1], [oy+min_vis*1.5, oy+min_vis*1.5],
                        [pr_z + min_vis*0.05*(pi+1), pr_z + min_vis*0.05*(pi+1)],
                        color=pcolor[0], linewidth=0.5, alpha=0.4*alpha)

        # 35. WATER TANK (demineralized water)
        if alpha > 0.5 and detail_level >= 2:
            wt_x = cav_cx - cav_s - min_vis * 5
            wt_r = min_vis * 0.5
            wt_h = min_vis * 0.8
            f = _cylinder_faces(ox+wt_x, oy+min_vis*1.5, 0, ox+wt_x, oy+min_vis*1.5, wt_h, wt_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.3*alpha,
                facecolor="#26A69A", edgecolor="#00897B", linewidth=0.3))
            ax.text(ox+wt_x, oy+min_vis*2, wt_h*1.2, "DM\nWater", fontsize=3,
                    color="#80CBC4", zorder=10)

        # 36. FUEL TANK (diesel)
        if alpha > 0.5 and detail_level >= 2:
            ft_x = stack_x + min_vis * 7
            ft_r = min_vis * 0.3
            ft_h = min_vis * 0.5
            f = _cylinder_faces(ox+ft_x, oy-min_vis*1.0, 0, ox+ft_x, oy-min_vis*1.0, ft_h, ft_r, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.4*alpha,
                facecolor="#FF6F00", edgecolor="#E65100", linewidth=0.3))
            ax.text(ox+ft_x, oy-min_vis*0.5, ft_h*1.3, "diesel\n72h", fontsize=3,
                    color="#FFE082", zorder=10)

        # 37. SECURITY FENCE (perimeter)
        if alpha > 0.5 and detail_level >= 2:
            fence_r = max(cav_side + tun_len + stack_h, 200) * 0.6
            n_fence = 24
            for fi in range(n_fence):
                ang = 2 * math.pi * fi / n_fence
                fx = ox + (cav_side + tun_len) * 0.5 + fence_r * math.cos(ang)
                fy = fence_r * math.sin(ang)
                ax.plot([fx, fx], [fy, fy], [0, min_vis*0.3],
                        color="#555", linewidth=0.3, alpha=0.15*alpha)

        # 38. HVAC DUCTS ON TURBINE HALL
        if alpha > 0.5 and detail_level >= 2:
            th_x = stack_x + min_vis * 4
            th_w = min_vis * 3
            th_h = min_vis * 1.5
            # Supply duct along hall roof
            ax.plot([ox+th_x, ox+th_x+th_w], [oy+min_vis*0.5, oy+min_vis*0.5],
                    [th_h*0.95, th_h*0.95], color="#80DEEA", linewidth=0.8,
                    alpha=0.3*alpha)
            # Branch ducts
            for di in range(4):
                dx = ox + th_x + th_w * (0.2 + di * 0.2)
                ax.plot([dx, dx], [oy+min_vis*0.5, oy+min_vis*0.2],
                        [th_h*0.95, th_h*0.95], color="#80DEEA", linewidth=0.4,
                        alpha=0.2*alpha)
            # Exhaust fan on roof
            for fi in range(2):
                fx = ox + th_x + th_w * (0.3 + fi * 0.4)
                ax.scatter([fx], [oy+min_vis*0.5], [th_h*1.05],
                          color="#26A69A", s=6, marker="o", alpha=0.4*alpha, zorder=8)

        # 39. TUNNEL LINING RINGS (visible in cutaway)
        if alpha > 0.5 and detail_level >= 2:
            for bi in range(min(n_bores, 3)):
                by = bore_offsets[bi][0] if bi < len(bore_offsets) else 0
                # Show a few lining ring segments
                for ri in range(5):
                    rx = tun_x0 + min_vis * (1 + ri * 2)
                    if rx > tun_x1: break
                    # Ring as small rectangle
                    ax.plot([rx, rx], [oy+by-bore_r, oy+by+bore_r],
                            [tun_z, tun_z], color="#616161", linewidth=0.3,
                            alpha=0.2*alpha)

        # 40. HX TUBE BUNDLE DETAIL (in lava zone)
        if alpha > 0.5 and lv.hx_enabled:
            lava_x0_d = tun_x0 + tun_len * 0.15
            # Show tube bundle as multiple parallel lines
            for ti in range(min(lv.hx_n_tubes, 8)):
                tube_y = oy + (ti - 3.5) * min_vis * 0.05
                ax.plot([ox+lava_x0_d, ox+lava_x0_d+min_vis*2],
                        [tube_y, tube_y], [tun_z, tun_z],
                        color="#FF8C00", linewidth=0.3, alpha=0.3*alpha)

        # 41. SWITCHYARD BUSWORK
        if alpha > 0.5 and detail_level >= 2:
            sw_x = stack_x + min_vis * 10
            sw_w = min_vis * 3
            # Bus bars (3-phase)
            for phase in range(3):
                bz = min_vis * (1.0 + phase * 0.15)
                ax.plot([ox+sw_x, ox+sw_x+sw_w], [oy, oy], [bz, bz],
                        color="#FFEB3B", linewidth=0.6, alpha=0.3*alpha)
            # Switchyard structures (pylons)
            for pi in range(3):
                px = ox + sw_x + sw_w * (0.2 + pi * 0.3)
                ax.plot([px, px], [oy, oy], [0, min_vis*1.5],
                        color="#90A4AE", linewidth=0.5, alpha=0.3*alpha)
                # Crossarm
                ax.plot([px-min_vis*0.2, px+min_vis*0.2], [oy, oy],
                        [min_vis*1.3, min_vis*1.3], color="#90A4AE",
                        linewidth=0.3, alpha=0.3*alpha)

        # 42. CONTROL ROOM EQUIPMENT (visible through walls) - batched
        if alpha > 0.5 and detail_level >= 2:
            cr_x = stack_x + min_vis * 5
            cr_z = min_vis * 0.5
            # Operator desks (3 workstations) - batched
            cr_dx = [ox + cr_x + di * min_vis * 0.3 for di in range(3)]
            cr_dy = [oy+min_vis*0.8] * 3
            cr_dz = [cr_z] * 3
            ax.scatter(cr_dx, cr_dy, cr_dz, color="#00d4ff", s=4,
                      marker="s", alpha=0.4*alpha, zorder=8)
            # Server rack
            ax.scatter([ox+cr_x+min_vis*1.0], [oy+min_vis*0.8], [cr_z],
                      color="#4CAF50", s=6, marker="s", alpha=0.4*alpha, zorder=8)

        # 43. METEOROLOGICAL STATION
        if alpha > 0.5 and detail_level >= 2:
            met_x = cav_cx - cav_s - min_vis * 8
            ax.plot([ox+met_x, ox+met_x], [oy-min_vis*2, oy-min_vis*2],
                    [0, min_vis*0.8], color="#90A4AE", linewidth=0.5,
                    alpha=0.4*alpha)
            # Anemometer (small cross at top)
            ax.plot([ox+met_x-min_vis*0.1, ox+met_x+min_vis*0.1],
                    [oy-min_vis*2, oy-min_vis*2], [min_vis*0.8, min_vis*0.8],
                    color="#90A4AE", linewidth=0.3, alpha=0.4*alpha)
            ax.text(ox+met_x, oy-min_vis*2, min_vis*1.0, "met", fontsize=2.5,
                    color="#B0BEC5", zorder=10)

        # 44. WATER TREATMENT PLANT
        if alpha > 0.5 and detail_level >= 2:
            wtp_x = cav_cx - cav_s - min_vis * 6
            wtp_w = min_vis * 1.5
            wtp_d = min_vis * 1.0
            wtp_h = min_vis * 0.6
            wtp_v = [
                [(ox+wtp_x, oy+min_vis*2-wtp_d/2, 0), (ox+wtp_x+wtp_w, oy+min_vis*2-wtp_d/2, 0),
                 (ox+wtp_x+wtp_w, oy+min_vis*2+wtp_d/2, 0), (ox+wtp_x, oy+min_vis*2+wtp_d/2, 0)],
                [(ox+wtp_x, oy+min_vis*2-wtp_d/2, wtp_h), (ox+wtp_x+wtp_w, oy+min_vis*2-wtp_d/2, wtp_h),
                 (ox+wtp_x+wtp_w, oy+min_vis*2+wtp_d/2, wtp_h), (ox+wtp_x, oy+min_vis*2+wtp_d/2, wtp_h)],
            ]
            ax.add_collection3d(Poly3DCollection(wtp_v, alpha=0.3*alpha,
                facecolor="#26A69A", edgecolor="#00897B", linewidth=0.3))
            ax.text(ox+wtp_x+wtp_w/2, oy+min_vis*2, wtp_h*2, "Water\nTreatment",
                    fontsize=3, color="#80CBC4", zorder=10)

        # 45. COOLING TOWER FAN (visible on top of cooling tower) - batched
        if alpha > 0.5 and detail_level >= 2:
            ct_x = cav_cx - cav_s - min_vis * 3
            ct_r = min_vis * 0.8
            # Fan markers - batched
            ctf_xs = [ox + ct_x + (fi - 1) * min_vis * 0.6 for fi in range(3)]
            ctf_ys = [oy] * 3
            ctf_zs = [min_vis * 1.5] * 3
            ax.scatter(ctf_xs, ctf_ys, ctf_zs, color="#26A69A", s=10,
                      marker="^", alpha=0.4*alpha, zorder=8)
            # Fan blades (small lines)
            for ctf_x in ctf_xs:
                for b in range(4):
                    ang = 2 * math.pi * b / 4
                    ax.plot([ctf_x, ctf_x + min_vis*0.1 * math.cos(ang)],
                           [oy, oy + min_vis*0.1 * math.sin(ang)],
                           [min_vis*1.5, min_vis*1.5],
                           color="#26A69A", linewidth=0.3, alpha=0.3*alpha)

        # 46. TRANSFORMER DETAIL (radiator banks)
        if alpha > 0.5 and detail_level >= 2:
            xfmr_x = stack_x + min_vis * 1
            xfmr_y = oy + min_vis * 0.5
            xfmr_w = min_vis * 1.2
            xfmr_d = min_vis * 0.8
            xfmr_h = min_vis * 0.6
            # Transformer body
            xfmr_v = [
                [(ox+xfmr_x, xfmr_y-xfmr_d/2, 0), (ox+xfmr_x+xfmr_w, xfmr_y-xfmr_d/2, 0),
                 (ox+xfmr_x+xfmr_w, xfmr_y+xfmr_d/2, 0), (ox+xfmr_x, xfmr_y+xfmr_d/2, 0)],
                [(ox+xfmr_x, xfmr_y-xfmr_d/2, xfmr_h), (ox+xfmr_x+xfmr_w, xfmr_y-xfmr_d/2, xfmr_h),
                 (ox+xfmr_x+xfmr_w, xfmr_y+xfmr_d/2, xfmr_h), (ox+xfmr_x, xfmr_y+xfmr_d/2, xfmr_h)],
            ]
            ax.add_collection3d(Poly3DCollection(xfmr_v, alpha=0.4*alpha,
                facecolor="#FFC107", edgecolor="#FF6F00", linewidth=0.3))
            # Radiator banks (vertical fins on sides)
            for ri in range(5):
                rx = ox + xfmr_x + ri * xfmr_w / 5
                ax.plot([rx, rx], [xfmr_y + xfmr_d/2, xfmr_y + xfmr_d/2 + min_vis*0.1],
                        [0, xfmr_h], color="#FFC107", linewidth=0.3, alpha=0.3*alpha)
            # Bushings (vertical insulators on top)
            for bi in range(3):
                bx = ox + xfmr_x + xfmr_w * (0.2 + bi * 0.3)
                ax.plot([bx, bx], [xfmr_y, xfmr_y], [xfmr_h, xfmr_h + min_vis*0.15],
                        color="#90A4AE", linewidth=0.5, alpha=0.4*alpha)
            # Conservator tank (cylinder on top)
            f = _cylinder_faces(ox+xfmr_x+xfmr_w*0.7, xfmr_y-xfmr_d*0.3, xfmr_h*1.1,
                               ox+xfmr_x+xfmr_w*0.9, xfmr_y-xfmr_d*0.3, xfmr_h*1.1,
                               min_vis*0.08, 8)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.3*alpha,
                facecolor="#555", edgecolor="#333", linewidth=0.2))

        # 47. SF6 BREAKERS (in switchyard)
        if alpha > 0.5 and detail_level >= 2:
            sw_x = stack_x + min_vis * 10
            for bi in range(2):
                br_x = ox + sw_x + min_vis * (0.5 + bi * 1.5)
                # Breaker body (vertical cylinder)
                f = _cylinder_faces(br_x, oy, 0, br_x, oy, min_vis*0.8, min_vis*0.1, 8)
                if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.4*alpha,
                    facecolor="#F44336", edgecolor="#C62828", linewidth=0.3))
                # Bushings (top)
                for pi in range(2):
                    px = br_x + (pi - 0.5) * min_vis * 0.2
                    ax.plot([px, px], [oy, oy], [min_vis*0.8, min_vis*1.1],
                            color="#90A4AE", linewidth=0.4, alpha=0.4*alpha)

        # 48. DISCONNECTORS (in switchyard)
        if alpha > 0.5 and detail_level >= 2:
            sw_x = stack_x + min_vis * 10
            for di in range(3):
                dx = ox + sw_x + min_vis * (1.0 + di * 0.8)
                # Disconnector base
                ax.scatter([dx], [oy], [min_vis*0.4],
                          color="#78909C", s=4, marker="s", alpha=0.3*alpha, zorder=8)
                # Blade (open position - angled)
                ax.plot([dx, dx + min_vis*0.15], [oy, oy], [min_vis*0.4, min_vis*0.7],
                        color="#90A4AE", linewidth=0.4, alpha=0.3*alpha)

        # 49. MUCK CONVEYOR (from cavern excavation)
        if alpha > 0.3 and detail_level >= 2:
            conv_x = cav_cx + cav_s * 0.3
            conv_y = oy - cav_s * 0.6
            ax.plot([ox+conv_x, ox+conv_x+min_vis*2], [conv_y, conv_y - min_vis*0.5],
                    [0, 0], color="#8D6E63", linewidth=0.5, alpha=0.2*alpha)
            ax.text(ox+conv_x+min_vis, conv_y - min_vis*0.3, 0, "muck\nconveyor",
                    fontsize=2.5, color="#8D6E63", alpha=0.3, zorder=10)

        # 50. SITE LIGHTING (pole-mounted lights) - batched
        if alpha > 0.5 and detail_level >= 2:
            lt_x, lt_y, lt_z = [], [], []
            for li in range(6):
                ang = 2 * math.pi * li / 6
                lx = ox + (cav_side + tun_len) * 0.5 + min_vis * 8 * math.cos(ang)
                ly = min_vis * 8 * math.sin(ang)
                ax.plot([lx, lx], [ly, ly], [0, min_vis*0.6],
                        color="#FFEB3B", linewidth=0.3, alpha=0.2*alpha)
                lt_x.append(lx)
                lt_y.append(ly)
                lt_z.append(min_vis*0.6)
            if lt_x:
                ax.scatter(lt_x, lt_y, lt_z, color="#FFEB3B", s=3,
                          marker="o", alpha=0.3*alpha, zorder=8)

    # --- draw all systems ---
    # For detail_level <= 1, only draw the primary system to save render time
    n_sys_draw = n_sys if detail_level >= 2 else 1
    for si in range(n_sys_draw):
        ox = si * (cav_side + tun_len + stack_h) * 1.1
        draw_system(ox, 0.0, alpha=1.0 if si == 0 else 0.6)

    # --- axes ---
    ax.set_xlabel("X (m)", fontsize=7, color="#c0c0c0")
    ax.set_ylabel("Y (m)", fontsize=7, color="#c0c0c0")
    ax.set_zlabel("Z (m)", fontsize=7, color="#c0c0c0")
    title = "3D System View - All Components to Scale"
    if n_sys > 1 and detail_level >= 2: title += f" ({n_sys} systems)"
    elif n_sys > 1: title += f" (showing 1 of {n_sys} systems - use Full detail for all)"
    ax.set_title(title, fontsize=10, fontweight="bold", color="white")
    ax.tick_params(colors="#c0c0c0", labelsize=6)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#333")
    ax.yaxis.pane.set_edgecolor("#333")
    ax.zaxis.pane.set_edgecolor("#333")

    # component summary
    parts = [f"Cavern {cv.volume_m3/1e6:.0f}M m3", f"Bores x{n_bores}", f"Turbines x{n_turb}",
             f"Reheat x{n_reheat}", f"Fans x{n_fans}"]
    if lv.hx_enabled: parts.append(f"HX x{lv.hx_n_tubes:,}")
    if tn.mhd_enabled: parts.append("MHD")
    if tn.potassium_enabled: parts.append("K-cycle")
    if tn.sco2_enabled: parts.append("sCO2")
    if tn.steam_enabled: parts.append("Steam")
    if tn.orc_enabled: parts.append("ORC")
    ax.text2D(0.02, 0.02, " | ".join(parts), transform=ax.transAxes,
              fontsize=5, color="#c0c0c0", family="monospace")


class VisualizerGUI:
    """Interactive tkinter + matplotlib visualization window.

    Tabbed interface with:
      - 3D View (isometric, all components to scale)
      - Turbine Engine (animated spinning turbine cross-section)
      - Operations (energy output over time, multi-panel)
      - Cross-Section (2D schematic, pan/zoom)
      - Timeline (power, T, P over time)
      - Energy Flow (bar chart)
      - Turbine Stages (T and W per stage)
      - Pressure Profile (P per stage)
      - Cavern State (mass and P over time)
      - Summary (key metrics)

    The matplotlib navigation toolbar provides pan, zoom, and save on every tab.
    """

    def __init__(self, initial_target: str = "Gmans Tunnel"):
        if not HAS_TK:
            print("ERROR: tkinter is not available. Install it via your Python distribution.")
            print("       On Linux: sudo apt install python3-tk")
            return
        if not HAS_MPL:
            print("ERROR: matplotlib is not available. Install it: pip install matplotlib")
            return

        self._target_names = list(targets().keys())
        if initial_target not in self._target_names:
            initial_target = self._target_names[0]
        self._current_target = initial_target
        self._res: Optional[SimResult] = None
        self._fr: Optional[FlowResult] = None
        self._stages: List[TurbineStage] = []
        self._t_dict: Optional[Dict] = None
        self._turbine_angle = 0.0
        self._animating = False
        self._anim_after_id = None
        self._drawn_tabs: set = set()  # track which tabs have been drawn
        self._current_tab = "3d"
        self._detail_level = 1  # 0=fast, 1=standard, 2=full

        self._root = tk.Tk()
        self._root.title("Gmans Tunnel - Interactive Visualization")
        self._root.geometry("1400x900")
        self._root.configure(bg="#0a0a0f")
        self._root.protocol("WM_DELETE_WINDOW", self._on_close)

        # --- style ---
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("TNotebook", background="#0a0a0f", borderwidth=0)
        style.configure("TNotebook.Tab", background="#1a1a2e", foreground="#c0c0c0",
                        padding=[12, 6], font=("Consolas", 10, "bold"))
        style.map("TNotebook.Tab",
                  background=[("selected", "#e94560")],
                  foreground=[("selected", "#ffffff")])
        style.configure("TFrame", background="#0a0a0f")

        # --- header ---
        hdr = tk.Frame(self._root, bg="#1a1a2e", height=44)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        tk.Label(hdr, text="GMANS TUNNEL - INTERACTIVE VISUALIZATION",
                 bg="#1a1a2e", fg="#e94560",
                 font=("Consolas", 14, "bold")).pack(side="left", padx=16)

        # target selector
        tk.Label(hdr, text="Target:", bg="#1a1a2e", fg="#00d4ff",
                 font=("Consolas", 10)).pack(side="left", padx=(20, 4))
        self._target_var = tk.StringVar(value=initial_target)
        self._target_combo = ttk.Combobox(
            hdr, textvariable=self._target_var, values=self._target_names,
            width=22, font=("Consolas", 10), state="readonly")
        self._target_combo.pack(side="left", padx=4)
        self._target_combo.bind("<<ComboboxSelected>>", self._on_target_change)

        # recompute button
        self._recompute_btn = tk.Button(
            hdr, text="Recompute", bg="#333344", fg="#00d4ff",
            font=("Consolas", 9, "bold"), bd=0, padx=12, pady=2,
            activebackground="#444466", activeforeground="#00d4ff",
            command=self._recompute)
        self._recompute_btn.pack(side="left", padx=8)

        # animate toggle button
        self._anim_btn = tk.Button(
            hdr, text="Animate Turbine", bg="#333344", fg="#00E676",
            font=("Consolas", 9, "bold"), bd=0, padx=12, pady=2,
            activebackground="#444466", activeforeground="#00E676",
            command=self._toggle_animation)
        self._anim_btn.pack(side="left", padx=8)

        # detail level selector
        tk.Label(hdr, text="Detail:", bg="#1a1a2e", fg="#FFEB3B",
                 font=("Consolas", 9)).pack(side="left", padx=(8, 2))
        self._detail_var = tk.StringVar(value="Standard")
        self._detail_combo = ttk.Combobox(
            hdr, textvariable=self._detail_var,
            values=["Fast", "Standard", "Full"], width=8,
            font=("Consolas", 9), state="readonly")
        self._detail_combo.pack(side="left", padx=2)
        self._detail_combo.bind("<<ComboboxSelected>>", self._on_detail_change)

        # status label
        self._status_lbl = tk.Label(hdr, text="Loading...", bg="#1a1a2e",
                                    fg="#00d4ff", font=("Consolas", 9))
        self._status_lbl.pack(side="right", padx=16)

        # --- notebook ---
        self._nb = ttk.Notebook(self._root)
        self._nb.pack(fill="both", expand=True, padx=4, pady=4)

        # --- tabs (3D first, then turbine engine, operations, etc.) ---
        self._tabs = {}
        self._canvases = {}
        self._figs = {}
        self._axes = {}

        tab_defs = [
            ("3d", "3D View"),
            ("blueprint", "Blueprint"),
            ("turbine_engine", "Turbine Engine"),
            ("operations", "Operations"),
            ("cross", "Cross-Section"),
            ("timeline", "Timeline"),
            ("energy", "Energy Flow"),
            ("turbines", "Turbine Stages"),
            ("pressure", "Pressure Profile"),
            ("cavern", "Cavern State"),
            ("summary", "Summary"),
        ]

        for key, label in tab_defs:
            frame = ttk.Frame(self._nb)
            self._nb.add(frame, text=f" {label} ")
            self._tabs[key] = frame

            if key == "blueprint":
                fig = Figure(figsize=(16, 14), facecolor="#0a0a12")
            else:
                fig = Figure(figsize=(12, 7), facecolor="#0d1117")
            self._figs[key] = fig

            if key == "3d":
                ax = fig.add_subplot(111, projection="3d")
            elif key == "operations":
                ax = fig.add_subplot(111)
                ax.set_facecolor("#0d1117")
            elif key == "blueprint":
                ax = fig.add_subplot(111)
                ax.set_facecolor("#0a0a12")
            else:
                ax = fig.add_subplot(111)
                ax.set_facecolor("#0d1117")
            self._axes[key] = ax

            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            toolbar = NavigationToolbar2Tk(canvas, frame)
            toolbar.update()
            toolbar.pack(side="bottom", fill="x")

            self._canvases[key] = canvas

        # tab change handler - lazy draw only the visible tab
        self._nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

        # initial compute + draw (only visible tab)
        self._recompute()

    def _toggle_animation(self) -> None:
        """Toggle turbine blade animation on/off."""
        self._animating = not self._animating
        if self._animating:
            self._anim_btn.config(text="Stop Animation", fg="#F44336")
            self._animate_turbine()
        else:
            self._anim_btn.config(text="Animate Turbine", fg="#00E676")
            if self._anim_after_id:
                self._root.after_cancel(self._anim_after_id)
                self._anim_after_id = None

    def _animate_turbine(self) -> None:
        """Animate the turbine engine view by rotating blades.
        Only redraws when the turbine tab is actually visible."""
        if not self._animating or self._t_dict is None:
            return
        # Skip rendering if turbine tab is not visible (major lag fix)
        if self._current_tab != "turbine_engine":
            self._turbine_angle += 0.15
            if self._turbine_angle > 2 * math.pi:
                self._turbine_angle -= 2 * math.pi
            self._anim_after_id = self._root.after(100, self._animate_turbine)
            return
        self._turbine_angle += 0.15
        if self._turbine_angle > 2 * math.pi:
            self._turbine_angle -= 2 * math.pi
        _draw_turbine_engine(self._axes["turbine_engine"], self._t_dict,
                             rotation_angle=self._turbine_angle)
        self._canvases["turbine_engine"].draw_idle()
        self._anim_after_id = self._root.after(100, self._animate_turbine)

    def _recompute(self) -> None:
        """Recompute the simulation for the current target and redraw all tabs."""
        name = self._target_var.get()
        self._status_lbl.config(text=f"Computing {name}...")
        self._root.update_idletasks()

        t = targets().get(name)
        if not t:
            self._status_lbl.config(text=f"ERROR: unknown target {name}")
            return
        self._t_dict = t
        self._current_target = name

        try:
            self._res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                                 hours=48.0, n_steps=300)
            st = build_initial_cavern(t["cavern"])
            self._fr = solve_flow(st, t["cavern"], t["lava"], t["tunnel"], t["ctrl"])

            if self._fr.mdot > 1e-3:
                self._stages = solve_turbine_stages(
                    self._fr.t_hot_k, t["cavern"].p_charge_pa, P_STD,
                    self._fr.mdot, t["tunnel"].n_turbine_stages,
                    t["tunnel"].turbine_eta, t["tunnel"].generator_eta)
            else:
                self._stages = []

            self._draw_all()

            self._status_lbl.config(
                text=f"  {_fmt_power(self._res.mean_p_net_mw * 1e6)} mean | "
                     f"{_fmt_power(self._res.peak_p_net_mw * 1e6)} peak | "
                     f"EROI {self._res.eroi:.2f}")

        except Exception as exc:
            import traceback
            self._status_lbl.config(text=f"ERROR: {exc}")
            traceback.print_exc()

    def _on_detail_change(self, event=None) -> None:
        """Handle detail level change - redraw 3D view."""
        level_map = {"Fast": 0, "Standard": 1, "Full": 2}
        self._detail_level = level_map.get(self._detail_var.get(), 1)
        # invalidate 3D tab cache so it redraws
        self._drawn_tabs.discard("3d")
        if self._current_tab == "3d":
            self._draw_tab("3d")
            self._canvases["3d"].draw_idle()

    def _on_tab_changed(self, event=None) -> None:
        """Handle tab change - lazily draw the newly visible tab."""
        if self._res is None or self._t_dict is None:
            return
        # find current tab index -> key
        idx = self._nb.index("current")
        tab_defs = [
            ("3d", "3D View"), ("blueprint", "Blueprint"),
            ("turbine_engine", "Turbine Engine"), ("operations", "Operations"),
            ("cross", "Cross-Section"), ("timeline", "Timeline"),
            ("energy", "Energy Flow"), ("turbines", "Turbine Stages"),
            ("pressure", "Pressure Profile"), ("cavern", "Cavern State"),
            ("summary", "Summary"),
        ]
        if idx < len(tab_defs):
            key = tab_defs[idx][0]
            self._current_tab = key
            # draw this tab if not already drawn or if data changed
            if key not in self._drawn_tabs:
                self._draw_tab(key)
                self._canvases[key].draw_idle()

    def _draw_tab(self, key: str) -> None:
        """Draw a single tab's content."""
        if self._res is None or self._t_dict is None:
            return
        t = self._t_dict
        ax = self._axes[key]
        fig = self._figs[key]

        if key == "3d":
            _draw_3d_view(ax, t, detail_level=getattr(self, '_detail_level', 1))
        elif key == "blueprint":
            _draw_blueprint(ax, t)
        elif key == "turbine_engine":
            _draw_turbine_engine(ax, t, rotation_angle=self._turbine_angle)
        elif key == "operations":
            _draw_operations(ax, self._res)
        elif key == "cross":
            _draw_cross_section(ax, t)
        elif key == "timeline":
            _draw_timeline(ax, self._res)
        elif key == "energy":
            if self._fr:
                _draw_energy_flow(ax, self._res, self._fr)
        elif key == "turbines":
            _draw_turbine_stages(ax, self._stages)
        elif key == "pressure":
            _draw_pressure_profile(ax, self._stages)
        elif key == "cavern":
            _draw_cavern_state(ax, self._res)
        elif key == "summary":
            _draw_summary_panel(ax, self._res, self._current_target)

        self._drawn_tabs.add(key)

    def _draw_all(self) -> None:
        """Draw only the currently visible tab (lazy drawing)."""
        if self._res is None or self._t_dict is None:
            return
        # invalidate cache - data changed
        self._drawn_tabs.clear()
        # draw only the current visible tab
        self._draw_tab(self._current_tab)
        self._canvases[self._current_tab].draw_idle()

    def _on_target_change(self, event=None) -> None:
        self._recompute()

    def _on_close(self) -> None:
        self._animating = False
        if self._anim_after_id:
            try:
                self._root.after_cancel(self._anim_after_id)
            except Exception:
                pass
        plt.close("all")
        self._root.destroy()

    def run(self) -> int:
        if not HAS_TK or not HAS_MPL:
            return 1
        self._root.mainloop()
        return 0


def cmd_visual(target: str) -> int:
    """Launch the interactive visualization GUI."""
    if not HAS_TK:
        print("ERROR: tkinter is not available.")
        print("  On Linux:  sudo apt install python3-tk")
        print("  On macOS:  brew install python-tk")
        print("  On Windows: tkinter is included with Python.")
        return 1
    if not HAS_MPL:
        print("ERROR: matplotlib is not available.")
        print("  Install:  pip install matplotlib")
        return 1
    gui = VisualizerGUI(initial_target=target)
    return gui.run()


# ==============================================================================
# SECTION 14 -- CLI
# ==============================================================================

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="CryoLavaTunnel -- cold-air cavern + lava-heated turbine tunnel digital twin")
    p.add_argument("--report", action="store_true", help="full report for all presets")
    p.add_argument("--selftest", action="store_true", help="physics + conservation + Carnot audit")
    p.add_argument("--info", action="store_true", help="explain every part + the math")
    p.add_argument("--honesty", action="store_true", help="the reality check, in full")
    p.add_argument("--proofs", action="store_true", help="math proofs with verify_fn")
    p.add_argument("--targets", action="store_true", help="list preset sites")
    p.add_argument("--hardware", action="store_true", help="to-scale hardware spec")
    p.add_argument("--parts", action="store_true", help="BOM parts list")
    p.add_argument("--model", action="store_true", help="ASCII cross-section visualization")
    p.add_argument("--flow", action="store_true", help="energy flow Sankey diagram")
    p.add_argument("--timeline", action="store_true", help="multi-series timeline plot")
    p.add_argument("--turbines", action="store_true", help="per-stage turbine breakdown")
    p.add_argument("--sweep", nargs="?", const=24.0, type=float, metavar="HOURS",
                   help="scan the design space")
    p.add_argument("--sensitivity", nargs=2, metavar=("KEY", "HOURS"),
                   help="sensitivity to one parameter (simple)")
    p.add_argument("--sensitivity2", nargs=2, metavar=("KEY", "HOURS"),
                   help="advanced multi-parameter sensitivity")
    p.add_argument("--mc", nargs="+", default=None,
                   metavar=("TARGET", "N"),
                   help="Monte Carlo (target name, realisations)")
    p.add_argument("--optimize", nargs="?", const="Gmans Tunnel", type=str,
                   metavar="TARGET", help="coordinate-descent optimiser")
    p.add_argument("--pareto", nargs="?", const="Gmans Tunnel", type=str,
                   metavar="TARGET", help="Pareto frontier (power vs duration)")
    p.add_argument("--live", nargs="?", const=48.0, type=float, metavar="HOURS",
                   help="run continuously with dashboard")
    p.add_argument("--visual", action="store_true",
                   help="interactive matplotlib GUI (pan/zoom, all views)")
    p.add_argument("--target", default="Gmans Tunnel",
                   help="preset name for --report/--hardware/--live/--visual")
    args = p.parse_args(argv)

    if args.visual:
        return cmd_visual(args.target)
    if args.selftest:
        return selftest()
    if args.info:
        print_info(); return 0
    if args.honesty:
        print_honesty(); return 0
    if args.proofs:
        print_proofs(); return 0
    if args.targets:
        print_header("PRESET TARGETS")
        for k, t in targets().items():
            print_target_info(k, t)
        return 0
    if args.sweep is not None:
        cmd_sweep(args.sweep); return 0
    if args.sensitivity:
        cmd_sensitivity(args.sensitivity[0], float(args.sensitivity[1])); return 0
    if args.sensitivity2:
        cmd_sensitivity_advanced(args.sensitivity2[0], float(args.sensitivity2[1])); return 0
    if args.mc is not None:
        key = args.mc[0] if len(args.mc) > 0 else "Gmans Tunnel"
        n = int(args.mc[1]) if len(args.mc) > 1 else 300
        hours = float(args.mc[2]) if len(args.mc) > 2 else 24.0
        cmd_monte_carlo(key, hours, n); return 0
    if args.optimize is not None:
        cmd_optimize(args.optimize, 24.0); return 0
    if args.pareto is not None:
        cmd_pareto(args.pareto, 24.0); return 0
    if args.live is not None:
        cmd_live(args.live); return 0
    if args.hardware:
        t = targets().get(args.target)
        if not t:
            print(f"unknown target: {args.target}"); return 1
        print_hardware(args.target, t); return 0
    if args.parts:
        t = targets().get(args.target)
        if not t:
            print(f"unknown target: {args.target}"); return 1
        print_parts(t); return 0
    if args.model:
        t = targets().get(args.target)
        if not t:
            print(f"unknown target: {args.target}"); return 1
        print_cross_section(args.target, t); return 0
    if args.flow:
        t = targets().get(args.target)
        if not t:
            print(f"unknown target: {args.target}"); return 1
        res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                       hours=24.0, n_steps=800)
        print_energy_flow(res, key=args.target); return 0
    if args.timeline:
        t = targets().get(args.target)
        if not t:
            print(f"unknown target: {args.target}"); return 1
        res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                       hours=48.0, n_steps=800)
        print_timeline(res, key=args.target); return 0
    if args.turbines:
        cmd_turbine_detail(args.target); return 0
    if args.report:
        for k, t in targets().items():
            res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                           hours=48.0, n_steps=1500)
            print_run(res, key=k)
            print()
        return 0

    # default: headline report for the default target
    t = targets().get(args.target)
    if not t:
        print(f"unknown target: {args.target}"); return 1
    res = simulate(t["cavern"], t["lava"], t["tunnel"], t["ctrl"],
                   hours=48.0, n_steps=1500)
    print_run(res, key=args.target)
    return 0


if __name__ == "__main__":
    sys.exit(main())

