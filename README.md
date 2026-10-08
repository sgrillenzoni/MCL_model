# MCL_model

Lumped-parameter vascular models for the pressure-controlled Windkessel chambers
of a mock circulatory loop (MCL).

Each model receives the **flow measured at the inlet of a chamber** and returns the
**pressure target** that the chamber must generate, i.e. the pressure that the
corresponding arterial bed would produce for that inflow.

```
Q_in [mL/s]  →  vascular model (Simscape)  →  P_target [mmHg]  →  chamber pressure control
```

## Contents

```
MCL_model/
├── MCL_model.prj                  MATLAB project (open this first)
├── models/
│   ├── MCL_vascular_lib.slx       library with the vascular blocks
│   └── test_upperbody.slx         test model: synthetic inflow → UpperBody → scope
├── parameters/
│   ├── init_model.m               master script: loads all parameter files
│   └── parameters_stergiopulos.m  R, L, C of the arterial tree + terminal loads
└── data/
    └── stergiopulos_RLC.xlsx      reference spreadsheet used to derive the parameters
```

## Blocks

| Block | Status | Arteries included |
|---|---|---|
| `UpperBody` | done, tested | innominate, left common carotid, left subclavian branches (Stergiopulos segments 3–13, 15–17, 19–25) |
| `Abdominal` | planned | |
| `LowerBody` | planned | |

- **Input:** `Q_in_mLs`, inflow to the chamber [mL/s]
- **Output:** `P_target_mmHg`, pressure setpoint for the chamber [mmHg]

## Model

- Arterial tree from Stergiopulos et al. (1992), Tables 1 (segments) and 2 (terminal impedances).
- Electrical analogy: voltage = pressure [mmHg], current = flow [mL/s].
- Each segment: R and L in series, C to ground. Each terminal segment is closed by a
  2-element Windkessel (RT to the venous pressure `P_venous`, CT to ground).
- R from Poiseuille, L = ρl/(πr²), with blood properties (μ = 0.04 P, ρ = 1.05 g/cm³);
  tapered segments are integrated assuming a linearly varying radius.
  Details in the header of `parameters_stergiopulos.m`.

## How to use

1. Open `MCL_model.prj`. The project sets the path and runs `init_model.m`, which loads all parameters.
2. Open `models/test_upperbody.slx` and run the simulation.
3. To use a block in another model, copy it from `MCL_vascular_lib.slx` and add `init_model;`
   to the model's **InitFcn** callback (Model Properties > Callbacks).

To add a new parameter file, create it in `parameters/` and add one line to `init_model.m`.

## Notes

- The Simulink-PS Converter at the block input filters the flow (time constant 5 ms).
  This is required because the flow is imposed on the branch inertances (L·dQ/dt)
  and it removes numerical spikes and sensor noise.
- Each block contains its own Solver Configuration: the physical network is self-contained.
- MATLAB version: **TODO**

## Reference

Stergiopulos N, Young DF, Rogge TR. Computer simulation of arterial flow with
applications to arterial and aortic stenoses. *J Biomech* 1992;25(12):1477–1488.

## Author

Samantha Grillenzoni, ETH Zurich
