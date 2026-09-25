# Voltage Stability Assessment of Transmission and Distribution Systems with Solar PV and Wind Integration

**Author:** Hussain Tak — [imhussaintak@gmail.com](mailto:imhussaintak@gmail.com)  
**Project:** M.Tech project (2025), EPES specialization, Department of Electrical Engineering, National Institute of Technology (NIT) Srinagar  
**Language:** MATLAB (R2021a or later, no toolboxes required)  
**License:** GNU GPL v2.0

---

## 1. Overview

This project studies how solar PV and wind generation affect the **voltage stability** of two standard test systems:

| System | Type | Renewables connected at |
|---|---|---|
| IEEE 14-bus | Meshed transmission network | Bus 14 |
| IEEE 33-bus (Baran & Wu) | Radial distribution feeder | Bus 33 |

Each system is analyzed in four scenarios:

1. **Base case**: no renewable generation
2. **PV only**: solar PV plant connected
3. **Wind only**: wind farm connected
4. **PV + Wind**: both plants connected

For every scenario the project computes:

- the **load flow**: Newton–Raphson for transmission, backward/forward sweep for distribution
- the **continuation power flow (CPF)**, which traces the full loadability (P–V / nose) curve and finds the **maximum loadability** λ_max
- the **voltage profile** and **network losses** at base load
- a **voltage stability index**: the L-index for the transmission system and the Fast Voltage Stability Index (FVSI) for the distribution feeder

The PV plant and the wind farm are sized from physical models. The PV plant uses a single-diode panel model, series/parallel array sizing, a boost converter and an inverter. The wind farm uses an aerodynamic turbine model, generator efficiency, a boost converter and an inverter.

---

## 2. Getting Started

### 2.1 Run the complete study

1. Open MATLAB and set the current folder to this project folder.
2. Open **`run_all.m`** and press **F5 (Run)**.

The full study takes about 10 seconds. When it finishes:

- **The Command Window** shows three result tables: the renewable plant ratings, the transmission results, and the distribution results.
- **Twelve figures** open as tabs in one Figures window. Every figure compares Base, PV, Wind and PV + Wind using the same colors.
- **All results are saved** to the `results/` folder: figures as PNG in `results/figures/`, tables as CSV.

| Fig. | Content |
|---|---|
| 1 | Transmission: loadability (P–V) curve at bus 14 |
| 2 | Transmission: bus voltage profile at base load |
| 3 | Transmission: network losses and slack generation from no load to voltage collapse |
| 4 | Transmission: L-index of the load buses |
| 5 | Transmission: effect of generator reactive-power limits |
| 6 | Distribution: loadability curves at bus 33 (renewable bus) and bus 18 (weakest bus) |
| 7 | Distribution: bus voltage profile at base load |
| 8 | Distribution: network losses and substation power from no load to voltage collapse |
| 9 | Distribution: FVSI of every branch |
| 10 | Active and reactive power injected by each renewable plant (both systems) |
| 11 | Maximum loadability λ_max of every scenario (both systems) |
| 12 | PV panel characteristics and plant sizes |

### 2.2 Loadability curve of any bus

After running `run_all`, type the following in the Command Window:

```matlab
explore_buses
```

The program asks for:

- the system (transmission or distribution)
- the x-axis (load multiplier λ or total load in MW)
- whether to save the figures
- the buses to plot, such as `9`, `4 9 14`, `10-14`, or `all`

It keeps asking for buses until you press Enter. It can also be called directly:

```matlab
explore_buses('transmission', [4 9 14])
explore_buses('distribution', [18 33], XAxis="MW", Save=true)
```

### 2.3 Changing the study

All study settings are in **`functions/study_config.m`**:

- the renewable connection bus
- the PV and wind plant ratings and power factors
- the buses plotted
- the CPF step sizes
- whether generator reactive limits are enforced

### 2.4 Verifying the code

```matlab
runtests('tests')
```

This runs 24 automated tests (about 20 s) that check the solvers against published results and against independent methods (Section 5).

---

## 3. Project Structure

```
run_all.m          Main script: runs the complete study (press F5)
explore_buses.m    Interactive loadability curves for any bus
functions/         All functions; study_config.m holds every setting
data/              IEEE 14-bus and IEEE 33-bus system data
results/           Output: figures/*.png, *.csv and *.mat
docs/              BUG_REPORT.md: corrections made to the earlier scripts
legacy/            Original scripts, kept unchanged for reference
```

The main functions in `functions/` are:

| Function | Purpose |
|---|---|
| `run_study` | Solves the four scenarios of one system |
| `run_cpf` | Continuation power flow (predictor–corrector) with generator Q-limits |
| `NRLF`, `nrlf_qlim` | Newton–Raphson load flow, without and with reactive limits |
| `build_bibc_bcbv`, `bfs_load_flow` | Backward/forward sweep load flow for radial feeders |
| `pv_farm_model`, `wind_farm_model` | PV plant and wind farm models |
| `l_index`, `fvsi_index` | Voltage stability indices |
| `make_figures` | All result figures |

---

## 4. Methodology

### 4.1 Load flow and continuation power flow

| Step | Transmission (IEEE 14-bus) | Distribution (IEEE 33-bus) |
|---|---|---|
| Load flow | Newton–Raphson (polar form) | Backward/forward sweep using the BIBC/BCBV matrices (Teng) |
| Loadability curve | Predictor–corrector CPF in three phases: λ → \|V\| → λ | Same CPF, which traces the full nose curve |
| Stability index | L-index (Kessel–Glavitsch) | FVSI (Musirin) |

The loading of the system is increased with the load multiplier λ:

```
P_injected(λ) = P_renewable + λ · (P_generation − P_load)
Q_injected(λ) = Q_renewable + λ · (Q_generation − Q_load)
```

λ = 1 is the base case, and λ_max is the loading at which voltage collapse occurs. The renewable plants inject **constant** power, which does not increase with λ.

### 4.2 Modelling assumptions

- Balanced, positive-sequence, steady-state (RMS) model
- Constant-power loads
- The slack bus (transmission) or substation (distribution) supplies the power mismatch
- Generator reactive-power limits are enforced in the transmission system through PV ↔ PQ bus switching. The slack bus has no limit.
- The renewable plants operate at fixed active power and fixed power factor, with no voltage control
- System base: 100 MVA. Distribution voltage base: 12.66 kV.

---

## 5. Validation

The following checks run automatically in `tests/`:

| Check | Result | Reference |
|---|---|---|
| IEEE 14-bus load flow | max voltage error 1.3×10⁻³ pu, max angle error 0.017° | Published IEEE 14-bus solution |
| Jacobian matrix | error < 10⁻⁶ | Numerical finite differences |
| IEEE 33-bus base case | V_min = 0.9133 pu at bus 18, losses 202.5 kW | Baran & Wu: 0.9131 pu, 202.7 kW |
| Backward/forward sweep vs Newton–Raphson | difference 8×10⁻¹² pu | Independent method |
| IEEE 14-bus CPF without Q-limits | λ_max = 4.060, V = 0.688 pu | Chatterjee's CPF: 4.0585, 0.688 pu |
| IEEE 14-bus CPF with Q-limits | λ_max = 1.77799, V = 0.6166 pu | Repeated Q-limited load flow: 1.77799, 0.6165 pu |
| IEEE 33-bus CPF | λ_max = 3.6247 | Repeated backward/forward sweep: 3.6242 |

---

## 6. Results

### 6.1 Renewable plants

| System | Plant | Size | P (MW) | Q (Mvar) | Power factor |
|---|---|---|---|---|---|
| Transmission | PV plant | 7562 panels (19 series × 398 strings) | 1.30 | 0.95 | 0.807 |
| Transmission | Wind farm | 5 turbines × 1.40 MW | 7.01 | 2.10 | 0.958 |
| Distribution | PV plant | 8721 panels (19 series × 459 strings) | 1.50 | 1.10 | 0.807 |
| Distribution | Wind farm | 3 turbines × 1.40 MW | 4.21 | 2.61 | 0.850 |

### 6.2 Transmission system (IEEE 14-bus, base load 259 MW)

| Scenario | Renewable P + jQ | λ_max | Load margin (MW) | V₁₄ at base load (pu) | Losses (MW) | L-index max | λ_max without Q-limits |
|---|---|---|---|---|---|---|---|
| Base | — | 1.778 | 201.5 | 1.036 | 13.39 | 0.077 | 4.060 |
| PV | 1.30 + j0.95 | 1.791 | 204.8 | 1.039 | 13.20 | 0.073 | 4.073 |
| Wind | 7.01 + j2.10 | 1.829 | 214.8 | 1.046 | 12.46 | 0.061 | 4.119 |
| PV + Wind | 8.32 + j3.06 | 1.841 | 217.8 | 1.050 | 12.30 | 0.060 | 4.130 |

### 6.3 Distribution system (IEEE 33-bus, base load 3.715 MW)

| Scenario | Renewable P + jQ | λ_max | V₁₈ at base load (pu) | V₃₃ at base load (pu) | Losses (kW) | FVSI max |
|---|---|---|---|---|---|---|
| Base | — | 3.625 | 0.913 | 0.917 | 202.5 | 0.067 |
| PV | 1.50 + j1.10 | 4.227 | 0.946 | 1.019 | 86.3 | 0.032 |
| Wind | 4.21 + j2.61 | 4.796 | 0.990 | **1.154** | **417.3** | 0.024 |
| PV + Wind | 5.71 + j3.71 | 5.066 | 1.012 | **1.224** | **798.2** | 0.023 |

### 6.4 Key findings

1. **Renewables increase loadability in both systems.**
   - In the transmission system the gain is modest: +0.7 % (PV) to +3.5 % (PV + Wind). The plants supply only 0.5–3.2 % of the load there.
   - In the distribution feeder the gain is large: +17 % (PV) to +40 % (PV + Wind). The plants there are comparable in size to the whole feeder load.

2. **Generator reactive-power limits determine the transmission loadability limit.** Without them, λ_max is 4.06. With them it falls to 1.78, a 56 % reduction.
   - All four generators reach their Q_max between λ = 1.08 and 1.22.
   - At base load only about 40 Mvar of generator reactive reserve remains (Fig. 5).

3. **The wind and PV + Wind plants are oversized for bus 33 of the distribution feeder.**
   - Their output (113 % and 154 % of the feeder load) causes reverse power flow.
   - Bus 33 rises to 1.15–1.22 pu at base load, above the usual 1.05 pu limit.
   - Losses are 2–4 times the base case.
   - So a higher λ_max alone does not mean a better operating point: plant size and location must also respect voltage limits.

4. **The PV plant alone improves every measure on the feeder:** losses −57 %, minimum voltage +0.033 pu, and λ_max +17 %.

5. **The effect of renewables on feeder losses depends on the loading level.**
   - At light load, the large plants increase losses because of reverse power flow.
   - Above λ ≈ 1.7–2.5 they reduce losses, because they then supply local load (Fig. 8).
   - Negative FVSI values on branches 26–33 show reactive power flowing back towards the substation.

6. **Weakest buses:** bus 14 in the transmission system and bus 18 in the distribution feeder.

---

## 7. Limitations and Future Work

- The renewable plants are modeled as constant PQ injections. Voltage-controlled operation or a volt-var characteristic would be more realistic.
- The PV and wind models represent a single steady-state operating point; irradiance and wind-speed profiles are not included.
- Plant size and location are fixed. Sweeping both would give a hosting-capacity map of each system.
- Slack-bus reactive limits and generator active-power limits are not modeled.

---

## 8. Author

**Hussain Tak**  
M.Tech, EPES specialization (2025)  
Department of Electrical Engineering  
National Institute of Technology (NIT) Srinagar

Email: [imhussaintak@gmail.com](mailto:imhussaintak@gmail.com)

This work was carried out as the author's M.Tech project at NIT Srinagar in 2025.

For questions, suggestions or collaboration, please get in touch by email.

---

## 9. Acknowledgements and License

The Newton–Raphson CPF core (`Calculate_Ybus`, `Jacobian_NRLF`, `NRLF`, `Calculate_PcalcQcalc`) is adapted from S. Chatterjee's open-source [continuation-power-flow](https://github.com/sayonsom/continuation-power-flow) (DOI: 10.5281/zenodo.1209882), released under the GNU GPL v2.0. This project is therefore distributed under the same license; see `LICENSE`.

## 10. References

1. M. Crow, *Computational Methods for Electric Power Systems*, CRC Press.
2. V. Ajjarapu and C. Christy, "The continuation power flow: a tool for steady state voltage stability analysis," *IEEE Trans. Power Systems*, 1992.
3. J.-H. Teng, "A direct approach for distribution system load flow solutions," *IEEE Trans. Power Delivery*, 2003.
4. M. E. Baran and F. F. Wu, "Network reconfiguration in distribution systems for loss reduction and load balancing," *IEEE Trans. Power Delivery*, 1989.
5. P. Kessel and H. Glavitsch, "Estimating the voltage stability of a power system," *IEEE Trans. Power Delivery*, 1986.
6. I. Musirin and T. K. Abdul Rahman, "Novel fast voltage stability index (FVSI) for voltage stability analysis in power transmission system," *SCOReD*, 2002.
