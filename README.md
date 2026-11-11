# MCL_model

Numerical model of a hybrid mock circulatory loop (MCL) with an anatomical aortic
phantom (Swiss Vascular): one heart chamber and three pressure-controlled
Windkessel chambers (upper body, abdominal, lower body). The bench is the
4-reservoir pressure-controlled HMC of Petrou et al. (2019).

The repository contains:

- a **1D arterial network** (PWDB, Charlton et al. 2019, baseline 25-year-old subject)
  solved with an explicit fixed-step scheme (1 ms), suitable for real time;
- the **chamber block `Twin1D`**: measured inflow of each chamber → pressure target
  of that chamber (to be used on the bench);
- the **full numerical twin `TwinFull1D`**: LVP → aortic root pressure, flows and
  pressures at the chamber inlets (reference for the bench);
- the **left-ventricle pressure generator** (Ferreira et al. 2005), calibrated
  against the twin on the Charlton reference values.

```
Ferreira generator --LVP--> [one-way valve] --> aorta --> chambers
                                                  |
bench:  Q measured at chamber inlets --> Twin1D --> P_target --> chamber pressure control
```

## Getting started

1. Open MATLAB in the repository root and run `startup` (adds all folders to the path).
2. Open `models/twin_numerico.slx` and run it. Its InitFcn runs
   `parameters_charlton; parameters_chamber1d; parameters_ferreira;` and then sets
   the calibrated generator values (see below).

Models are saved in **R2024b**. The lab PC runs **R2018**: use
*Save As → Export to previous version* (R2018a/b) and check that the model opens and
compiles there. Only blocks available in R2018 are used.

## Structure

```
params/        parameters_charlton.m    PWDB network (geometry, stiffness, Windkessels, globals)
               parameters_chamber1d.m   builds S_upper, S_abdominal, S_lower, S_full (1D sub-trees)
               parameters_ferreira.m    Ferreira generator, values of the paper
               parameters_backflowcontrol.m  return control (one-beat mean, pump LUT feedforward, level PI)
solver1d/      charlton1d_setup.m       precomputes the 1D network (cells, faces, nodes, Windkessels)
               charlton1d_step.m        one explicit time step (also used inside the blocks)
blocks/        twin_chambers_block.m    body of MATLAB Function "Twin1D"
               twin_full_block.m        body of MATLAB Function "TwinFull1D"
               chamber1d_block.m        body of a single-chamber block
models/        twin_numerico.slx        Ferreira_LV (2e-4 s) -> Rate Transition -> TwinFull1D / Twin1D (1e-3 s)
               ferreira_charlton.slx    calibrated generator alone
               elastance_model.slx      development model of the Ferreira_LV block
models/tests/  test_twin1d.slx, test_twinfull.slx, test_chamber_block.slx
               backflow_test.slx, backflowblocks_test.slx   return control with a simulated chamber
calibration/   calib_ferreira_grid.m, calib_ferreira_emax.m, calib_ferreira_3p.m (+ results/*.mat)
validation/    run_charlton_1d.m        offline run of the full network with the PWDB inflow
               compare_pwdb_1d.m        comparison with the PWDB reference simulation
               test_chamber_1d.m        offline check of the three chamber sub-trees
               run_ferreira_verification.m, ferreira_metrics.m, twin_metrics.m, make_charlton_inflow.m
               figures/*.fig            validation figures
data/          pwdb_geo0001.csv         PWDB geometry, baseline subject (116 segments)
external/pwdb/ AorticFlowWave_orig.m    aortic inflow wave generator from the PWDB code (Charlton et al.)
```

## Signals

| Block | Input | Output | Sample time |
|---|---|---|---|
| `Twin1D` | `Q_upper`, `Q_abdominal`, `Q_lower` [mL/s], flow measured into each chamber (positive entering) | `P_upper`, `P_abdominal`, `P_lower` [mmHg], pressure targets | 1e-3 s |
| `TwinFull1D` | `LVP` [mmHg], upstream of the aortic valve | `y` 8×1: `[P_root; Q_upper; Q_abdominal; Q_lower; P_upper; P_abdominal; P_lower; Q_ao]` [mmHg, mL/s] | 1e-3 s |
| `Ferreira_LV` | – | `LVV` [mL], `LVV-V0` [mL], `LVP` [mmHg] | 2e-4 s |

The structs `S_*` are passed to the MATLAB Function blocks as **non-tunable parameters**.
The aortic valve of the bench is modelled as an ideal diode with `R_V = 0.003 mmHg·s/mL`
(`S_full.Rv`, in `parameters_chamber1d.m`), solved inside `charlton1d_step`.
`K_AORTA` in `parameters_chamber1d.m` scales the stiffness of the aortic segments
(1 = Charlton); it only affects `TwinFull1D`.

## Calibrated generator

The Ferreira generator is calibrated against `TwinFull1D` on the Charlton baseline
values (SV 67.1 mL, peak flow 421 mL/s, LVET 280 ms; weights [1 0.3 1] on the squared
relative errors). Exponents of the double-Hill function and Emin are those of the paper.
Values set in the InitFcn of `twin_numerico` and `ferreira_charlton`, after
`parameters_ferreira` (which keeps the paper values):

```matlab
KT   = 1.44;  Tmax = KT*Tmax;   % Tmax = 1.44*(0.2 + 0.15*tc) = 0.466 s
Emax = 1.40;                    % [mmHg/mL] (paper: 2.0)
LAP  = 8.50;                    % [mmHg]
```

| | SV [mL] | Peak flow [mL/s] | LVET [ms] | Root P [mmHg] |
|---|---|---|---|---|
| Paper parameters | 70.2 | 739 | 212 | 108.0 / 77.7 |
| **Calibrated** | **64.9** | **458** | **290** | **99.1 / 75.5** |
| Charlton 2019 | 67.1 | 421 | 280 | 101.7 / 76.6 |

Note: in the calibration grids `Tmax` was scaled twice, so the "KT" axis of the saved
grids corresponds to an effective factor KT². The final value (1.44) was checked with
single scaling.

## Validation status

- 1D solver vs PWDB reference simulation: `validation/run_charlton_1d.m`, `compare_pwdb_1d.m`.
- `Twin1D` (measured Q → P_target) vs `TwinFull1D`: differences < 1 mmHg on all three chambers.

`compare_pwdb_1d.m` needs the PWDB reference files (`sim_1_*.his`, a few MB each,
not tracked) in `data/pwdb_sim1/`; they can be downloaded from the PWDB Zenodo record
(Charlton et al. 2019).

## Open points

- Save the calibrated generator as a separate library block.
- Solver of the lab models (fixed step 2e-4 s): if discrete, replace the integrators of
  Ferreira_LV with Discrete-Time Integrators. Check computation time on the MicroLabBox.
- Integration with the existing chamber pressure control and the return control (next commits).

## References

- Charlton PH et al. (2019). Modeling arterial pulse waves in healthy aging: a database
  for in silico evaluation of hemodynamics and pulse wave indexes. *Am J Physiol Heart
  Circ Physiol* 317:H1062–H1085.
- Ferreira A, Chen S, Simaan MA, Boston JR, Antaki JF (2005). A nonlinear state-space model
  of a combined cardiovascular system and a rotary pump. *Proc. 44th IEEE CDC-ECC*, 897–902.
- Petrou A, Granegger M, Meboldt M, Schmid Daners M (2019). A versatile hybrid mock
  circulation for hydraulic investigations of active and passive cardiovascular implants.
  *ASAIO J* 65(5):495–502. doi:10.1097/MAT.0000000000000851

## History

The earlier Simscape lumped-parameter version (Stergiopulos / Westerhof networks,
`MCL_vascular_lib.slx`) is in the history of this repository: tag `v0.0-simscape-stergiopulos`.
