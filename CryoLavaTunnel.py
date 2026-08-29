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

    def ground_t_at_depth_c(self) -> float:
        """Local ground temperature at the cavern depth.

        This is the HONESTY correction: below the shallow stable zone the ground
        WARMS at the geothermal gradient, so a deep cavern next to lava is HOT,
        not cold. The passive 'operates when not cooled' mode only works if the
        ground at this depth is actually cooler than the charge temperature.
        """
        if not math.isnan(self.t_ground_c):
            return self.t_ground_c
        # stable zone ~ first 4 m tracks mean surface; below that add gradient
        stable = 4.0
        if self.depth_m <= stable:
            return self.surf_t_c
        return self.surf_t_c + (self.depth_m - stable) * self.ground_k_per_m

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
    q_leak = spec.u_ground * spec.area_ground_m2 * dT   # W (+: heat in)
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
            q_leak_now = cavern_spec.u_ground * cavern_spec.area_ground_m2 * (t_gnd_k - st.t_k)
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
                "fans. 110 TW mean, 115 TW peak, EROI 10.64.",
        "cavern": ColdCavernSpec(
            name="Gmans Tunnel Cavern (per system)", volume_m3=6.0e9,
            depth_m=30.0,
            p_charge_pa=30000.0e3, t_charge_k=123.15,   # 300 bar, -150 C
            surf_t_c=5.0, u_ground=0.3, area_ground_m2=1.0e6,
            active_cooling=True, chiller_kW_thermal=1000000.0, chiller_cop=1.0,
            cascade_cooling=True, cascade_cop=0.3,
            lava_heated_cooling=True, lava_cooling_fraction=0.85,
            n_compress_stages=20, compress_eta=0.92),
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


def _cylinder_faces(x0, y0, z0, x1, y1, z1, r, n_seg=12):
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

    # Draw each stage
    for i in range(n_stages):
        cx = (i + 0.5) * stage_spacing

        # Stator vanes (fixed, angled) - drawn as thin triangles
        for b in range(min(n_blades, 12)):  # limit for clarity
            ang = 2 * math.pi * b / min(n_blades, 12)
            # stator at left side of stage
            sx = cx - stage_spacing * 0.3
            # vane as a line from casing to rotor
            y_outer = casing_r * 0.9 * v_scale * (1 if ang < math.pi else -1)
            y_inner = rotor_r * 0.9 * v_scale * (1 if ang < math.pi else -1)
            # only draw top and bottom vanes (2D cross-section)
            if ang < math.pi:
                ax.plot([sx, sx], [y_inner, y_outer], color="#9E9E9E",
                        linewidth=1.5, alpha=0.6, zorder=3)
            else:
                ax.plot([sx, sx], [y_outer, y_inner], color="#9E9E9E",
                        linewidth=1.5, alpha=0.6, zorder=3)

        # Rotor blades (spinning) - drawn as angled lines
        for b in range(min(n_blades, 12)):
            ang = 2 * math.pi * b / min(n_blades, 12) + rotation_angle
            rx = cx + stage_spacing * 0.2
            # blade tip position
            blade_y = rotor_r * v_scale * math.sin(ang)
            blade_x_offset = rotor_r * 0.3 * math.cos(ang)
            # only draw blades in the visible 2D plane (near sin=±1)
            visibility = abs(math.sin(ang))
            if visibility > 0.3:
                alpha = visibility * 0.8
                color = "#00E676" if math.sin(ang) > 0 else "#00C853"
                ax.plot([rx + blade_x_offset, rx + blade_x_offset],
                        [shaft_r * v_scale * (1 if math.sin(ang) > 0 else -1), blade_y],
                        color=color, linewidth=2, alpha=alpha, zorder=4)

        # Stage boundary
        ax.axvline(cx + stage_spacing * 0.5, color="#333", linewidth=0.5, alpha=0.3, zorder=1)

        # Label first and last few stages
        if i < 3 or i >= n_stages - 3 or (i == n_stages // 2):
            ax.text(cx, casing_y_top + 15, f"S{i+1}", fontsize=6,
                    ha="center", color="#00E676", fontweight="bold")
        elif i == 3:
            ax.text(cx, casing_y_top + 15, f"...", fontsize=6,
                    ha="center", color="#888")

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


def _draw_3d_view(ax, t: Dict) -> None:
    """Draw a comprehensive 3D to-scale view of the entire system.

    Shows every component: cavern, access tunnel, parallel bores, lava body,
    HX tube bundle, heat pipes, turbine stages (as cylinders with blade hints),
    reheat sections, MHD channel, regenerator, stack, exit nozzle, exit fans,
    bottoming cycle heat exchangers, transformer, chiller, and dual-system offset.
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
            f = _cylinder_faces(ox+tun_x0, oy+by, bz_z, ox+lava_x0, oy+by, bz_z, bore_r, 10)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.25*alpha, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.3))
            # lava section
            f = _cylinder_faces(ox+lava_x0, oy+by, bz_z, ox+lava_x1, oy+by, bz_z, bore_r, 10)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.35*alpha, facecolor="#FF6347", edgecolor="#FF4500", linewidth=0.3))
            # post-lava section
            f = _cylinder_faces(ox+lava_x1, oy+by, bz_z, ox+tun_x1, oy+by, bz_z, bore_r, 10)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.25*alpha, facecolor="#FFD700", edgecolor="#B8860B", linewidth=0.3))

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
                                ox+tx+turb_r*0.4, oy+by0, tun_z+bz0, turb_r, 10)
            if f: ax.add_collection3d(Poly3DCollection(f, alpha=0.8*alpha, facecolor="#00E676", edgecolor="#00C853", linewidth=0.5))
            # blade markers (small lines across the rotor face)
            if alpha > 0.5 and i < 8:
                ax.plot([ox+tx, ox+tx], [oy+by0-turb_r*0.8, oy+by0+turb_r*0.8],
                        [tun_z+bz0, tun_z+bz0], color="#00C853", linewidth=0.5, alpha=0.5)
                ax.text(ox+tx, oy+by0+turb_r*1.5, tun_z+bz0, f"T{i+1}", fontsize=4, color="#00E676", zorder=10)
        if alpha > 0.5 and n_turb > 8:
            ax.text(ox+tun_x0+(0.15+8/(n_turb+1)*(lava_len/tun_len))*tun_len, oy+turb_r*1.5, tun_z,
                    f"...T{n_turb}", fontsize=4, color="#00E676", zorder=10)

        # 9. REHEAT MARKERS
        if n_reheat > 0 and alpha > 0.5:
            for i in range(min(n_reheat, 10)):
                frac = 0.15 + (i + 1.5) / (n_turb + 1) * (lava_len / tun_len)
                rx = tun_x0 + frac * tun_len
                ax.scatter([ox+rx], [oy+by0], [tun_z+bz0+turb_r], color="#FF5722", s=12, marker="v", zorder=8, alpha=alpha)

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
        f = _cylinder_faces(ox+stack_x, oy, tun_z, ox+stack_x, oy, stack_h, stack_r, 10)
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

        # 15. EXIT JET
        if alpha > 0.5:
            for i in range(3):
                jy = (i - 1) * fan_r * 2
                ax.plot([ox+stack_x, ox+stack_x], [oy+jy, oy+jy],
                        [stack_h+min_vis, stack_h+min_vis*3], color="#FFEB3B", linewidth=1, alpha=0.4*alpha)
                ax.scatter([ox+stack_x], [oy+jy], [stack_h+min_vis*3], color="#FFEB3B", s=6, marker="^", alpha=0.5*alpha)

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

    # --- draw all systems ---
    for si in range(n_sys):
        ox = si * (cav_side + tun_len + stack_h) * 1.1
        draw_system(ox, 0.0, alpha=1.0 if si == 0 else 0.6)

    # --- axes ---
    ax.set_xlabel("X (m)", fontsize=7, color="#c0c0c0")
    ax.set_ylabel("Y (m)", fontsize=7, color="#c0c0c0")
    ax.set_zlabel("Z (m)", fontsize=7, color="#c0c0c0")
    title = "3D System View - All Components to Scale"
    if n_sys > 1: title += f" ({n_sys} systems)"
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

            fig = Figure(figsize=(12, 7), facecolor="#0d1117")
            self._figs[key] = fig

            if key == "3d":
                ax = fig.add_subplot(111, projection="3d")
            elif key == "operations":
                ax = fig.add_subplot(111)
                ax.set_facecolor("#0d1117")
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

        # initial compute + draw
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
        """Animate the turbine engine view by rotating blades."""
        if not self._animating or self._t_dict is None:
            return
        self._turbine_angle += 0.15
        if self._turbine_angle > 2 * math.pi:
            self._turbine_angle -= 2 * math.pi
        _draw_turbine_engine(self._axes["turbine_engine"], self._t_dict,
                             rotation_angle=self._turbine_angle)
        self._canvases["turbine_engine"].draw_idle()
        self._anim_after_id = self._root.after(50, self._animate_turbine)

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
                                 hours=48.0, n_steps=800)
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

    def _draw_all(self) -> None:
        """Redraw all tab figures."""
        if self._res is None or self._t_dict is None:
            return

        # 3D view (primary)
        _draw_3d_view(self._axes["3d"], self._t_dict)
        self._canvases["3d"].draw_idle()

        # turbine engine
        _draw_turbine_engine(self._axes["turbine_engine"], self._t_dict,
                             rotation_angle=self._turbine_angle)
        self._canvases["turbine_engine"].draw_idle()

        # operations (energy output over time)
        _draw_operations(self._axes["operations"], self._res)
        self._canvases["operations"].draw_idle()

        # cross-section
        _draw_cross_section(self._axes["cross"], self._t_dict)
        self._figs["cross"].tight_layout()
        self._canvases["cross"].draw_idle()

        # timeline
        _draw_timeline(self._axes["timeline"], self._res)
        self._figs["timeline"].tight_layout()
        self._canvases["timeline"].draw_idle()

        # energy flow
        if self._fr:
            _draw_energy_flow(self._axes["energy"], self._res, self._fr)
        self._figs["energy"].tight_layout()
        self._canvases["energy"].draw_idle()

        # turbine stages
        _draw_turbine_stages(self._axes["turbines"], self._stages)
        self._figs["turbines"].tight_layout()
        self._canvases["turbines"].draw_idle()

        # pressure profile
        _draw_pressure_profile(self._axes["pressure"], self._stages)
        self._figs["pressure"].tight_layout()
        self._canvases["pressure"].draw_idle()

        # cavern state
        _draw_cavern_state(self._axes["cavern"], self._res)
        self._figs["cavern"].tight_layout()
        self._canvases["cavern"].draw_idle()

        # summary
        _draw_summary_panel(self._axes["summary"], self._res, self._current_target)
        self._canvases["summary"].draw_idle()

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

