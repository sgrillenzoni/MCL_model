# MCL_model

Simulink/Simscape models for a hybrid mock circulatory loop (MCL) with an anatomical
aortic model (Swiss Vascular), one heart chamber and three pressure-controlled
Windkessel chambers (upper body, abdominal, lower body).

The repository contains:
- **vascular blocks** that receive the flow measured at the inlet of a chamber and
  return the pressure target that the chamber must generate;
- a **left-ventricle waveform generator** (Ferreira et al. 2005) that drives the
  heart chamber;
- a **digital twin** of the whole loop, used to produce reference curves to compare
  with the testbench measurements.

```
Q_in [mL/s]  →  vascular model (Simscape)  →  P_target [mmHg]  →  chamber pressure control
```

## Contents

```
MCL_model/
├── MCL_model.prj                  MATLAB project (open this first)
├── models/
│   ├── MCL_vascular_lib.slx               library: UpperBody, Abdominal, LowerBody, Ferreira_LV
│   ├── mockloop_numerical.slx             digital twin of the full loop
│   ├── windkessel_upperbody_test.slx      test model: synthetic inflow → UpperBody → scope
│   ├── windkessel_abdominal_test.slx      test model: synthetic inflow → Abdominal → scope
│   ├── windkessel_lowerbody_test.slx      test model: synthetic inflow → LowerBody → scope
│   └── elastance_model.slx                development model of the Ferreira LV generator
├── parameters/
│   ├── init_model.m               master script: loads all parameter files, sets MODE
│   ├── parameters_stergiopulos.m  R, L, C of the arterial tree + terminal loads + P_venous
│   └── parameters_ferreira.m      parameters of the Ferreira LV model
└── data/
    └── stergiopulos_RLC.xlsx      reference spreadsheet used to derive the parameters
```

## Blocks

| Block | Status | Content |
|---|---|---|
| `UpperBody` | done, tested | innominate, left common carotid, left subclavian branches (Stergiopulos segments 3–13, 15–17, 19–25) |
| `Abdominal` | done, tested | renal arteries + splanchnic branches (celiac, hepatic, gastric, splenic, superior and inferior mesenteric) + intercostals |
| `LowerBody` | done, tested | iliac arteries and legs |
| `Ferreira_LV` | done, tested | left-ventricle waveform generator (see below) |

### Vascular blocks (UpperBody, Abdominal, LowerBody)

**Interface**
- **Input** `Q_in` (Simulink): inflow to the chamber [mL/s], from the flow sensor.
- **Physical port** `Q_numerical_*`: direct electrical connection to the aorta (digital twin only).
- **Output** `P_target`: pressure at the inlet of the vascular bed [mmHg], i.e. the chamber setpoint.

**Model**
- Arterial tree from Stergiopulos et al. (1992), Table 1 (segments) and Table 2 (terminal impedances).
- Electrical analogy: voltage = pressure [mmHg], current = flow [mL/s]
  (in Simulink-PS / PS-Simulink converters: 1 V = 1 mmHg, 1 A = 1 mL/s).
- Each segment: R and L in series, C to ground. Each terminal segment is closed by a
  2-element Windkessel (RT to the venous pressure `P_venous` = 5 mmHg, CT to ground).
- R from Poiseuille, L = ρl/(πr²), with blood properties (μ = 0.04 P, ρ = 1.05 g/cm³);
  tapered segments are integrated assuming a linearly varying radius.
  Details in the header of `parameters_stergiopulos.m`.
- Input filter 5 ms on the Simulink-PS converter (Q is imposed on the branch
  inertances, L·dQ/dt): avoids numerical spikes.
- Output filter: 2nd-order Butterworth, 10 Hz: removes the L-C resonances of the
  lumped network before the signal is sent to the chamber.

**Operating mode** (mask parameter `mode`, default value = workspace variable `MODE`
set in `init_model.m`)

| MODE | Use | Flow enters the network via | Internal Solver Configuration |
|---|---|---|---|
| 1 | numerical (digital twin, ideal) | physical port, direct connection | off (main model provides it) |
| 2 | testbench / standalone (Q from the flow sensor, test models) | `Q_in` → controlled current source | on |

Implementation:
- the physical port reaches the network through a Simscape Switch, closed when
  `double(mode==1)` > 0.5;
- the current source is always connected; its command is `Q_in * (mode>=2)`, so in
  MODE 1 it imposes zero current (open circuit). A switch is not used on this branch
  because an ideal current source in series with an open switch gives a diverging voltage;
- the internal Solver Configuration is commented out by the mask initialization code
  when `mode ~= 2` (option "Allow library block to modify its contents" enabled).
  A physical network must have exactly one Solver Configuration: in MODE 1 the
  block network is part of the main network.
- the mask initialization code is a plain script (not the function or class-based
  format of recent releases) for R2018 compatibility;
- `OpenFcn = open_system(gcb,'force')`: double click opens the block content.

Constants that feed a Simulink-PS converter must be `double` (e.g. `double(mode==1)`),
not boolean.

### Ferreira_LV

Left-ventricle waveform generator, Ferreira et al. (2005), reduced model.
- Open loop: it reproduces a physiological waveform and does not react to the
  downstream circuit.
- 3 states (z1 = LVV − V0, arterial pressure AP, aortic flow QA); left atrial pressure
  constant, LAP = 8 mmHg (own parameter, independent of `P_venous`).
- R1 and C3 derived from the Stergiopulos tree (`R_periph`, `C_total`).
- V0 = 10 mL, fixed placeholder (to be replaced by height-sensor calibration).
- Outputs: `LVV` [mL], `z1` = LVV − V0 [mL], `LVP` [mmHg].
- Built with R2018-compatible blocks only (Fcn instead of MATLAB Function).
- Initial conditions set to steady state to avoid start-up transients.
- Note: Ferreira parameters use generic names (R1, R2, C3, L, HR, …); check for name
  collisions with other parameter files.

## Digital twin (`mockloop_numerical.slx`)

Full loop in simulation, used as reference for the testbench results.

```
Ferreira_LV (LVP) → controlled pressure source → aortic valve (diode)
  → ascending aorta → aortic arch ─┬─ current sensor → UpperBody
                                   └─ descending aorta ─┬─ current sensor → Abdominal
                                                        └─ … → current sensor → LowerBody
```

- **Heart**: LVP from `Ferreira_LV` imposed by a controlled voltage source (open loop).
- **Aortic valve**: Simscape Diode, forward voltage 0.01 V (= 0.01 mmHg, must be > 0),
  small on-resistance, small off-conductance. Represents the physical check valve at
  the heart-chamber outlet (the valves inside `Ferreira_LV` only shape the waveform).
- **Aorta**: lumped R-L-C cells with Stergiopulos values.
- **Branches**: vascular blocks in MODE 1, connected through their physical port; a
  current sensor in series measures the flow into each bed, a voltage sensor in
  parallel reads the pressure at each branch node.
- One Solver Configuration at the top level; stiff solver (ode15s or ode23t) because
  of the valve.
- Logged signals: raw flows (To Workspace).

**Reference values (steady state, approx.)**
- aortic pressure ≈ 70/135 mmHg; distal pressures with larger pulse pressure and
  a delay of a few tens of ms (peripheral amplification);
- mean flows ≈ 17 mL/s upper body, ≈ 55 mL/s abdominal, ≈ 16 mL/s lower body;
  total ≈ 90 mL/s, to be checked against the mean flow through the valve.
- periodic steady state after ~2 s (upper/lower body) and ~5 s (abdominal): compute
  means from 6 s on, over an integer number of beats.

**Known limitations**
- One cell per vessel: long segments resonate at frequencies inside the physiological
  band (e.g. subclavian B ≈ 2 Hz, radial ≈ 6.5 Hz) and short ones at 10–20 Hz. This
  gives oscillations and diastolic backflow larger than physiological in the branch
  flows. Kept on purpose (model faithful to Stergiopulos); filtering of the logged
  signals for the comparison with the testbench is still to be implemented.
- The aorta uses blood properties; on the testbench the aorta is physical and the
  working fluid is water-glycerol.
- Tubing between the aortic outlets and the chambers is not modelled.
- The hybrid coupling (measured flow → block → imposed pressure, with filter delays)
  was unstable when tested inside the twin: delayed feedback around the aortic
  inertance. Relevant for the design of the testbench control.

## How to use

1. Open `MCL_model.prj`. The project sets the path and runs `init_model.m`, which
   loads all parameters and sets `MODE`.
2. Digital twin: open `models/mockloop_numerical.slx` with `MODE = 1` and run.
3. Single block test: open `models/windkessel_upperbody_test.slx` (or the abdominal /
   lowerbody test); its InitFcn runs `init_model; MODE = 2;`.
4. To use a block in another model, copy it from `MCL_vascular_lib.slx` and add
   `init_model;` to the model's **InitFcn** callback (Model Properties > Callbacks).
   If the block is used without a connected aorta, set `MODE = 2`.

To add a new parameter file, create it in `parameters/` and add one line to `init_model.m`.

## Compatibility and version control

- Developed in R2024b; models must open in R2018 on the lab PC: only blocks available
  in R2018, export with *Export Model to → Previous Version* (models and library), and
  keep the parameter scripts together with the models.
- Version control with GitHub Desktop (custom `.gitignore`).
