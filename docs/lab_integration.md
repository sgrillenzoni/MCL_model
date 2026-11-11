# Integration with the bench (to fill in at the lab)

For each arrow: signal name in the lab model, unit, sample rate.

| From | To | Signal (twin side) | Lab signal name | Unit | Rate |
|---|---|---|---|---|---|
| Flow sensors at chamber inlets | `Twin1D` | `Q_upper`, `Q_abdominal`, `Q_lower` [mL/s] | | | 1 kHz |
| `Twin1D` | Pressure control of the 3 chambers | `P_upper`, `P_abdominal`, `P_lower` [mmHg] | | | 1 kHz |
| Calibrated Ferreira generator | Heart chamber pressure setpoint (open loop) | `LVP` [mmHg] | | | 5 kHz |
| Level sensors + mean flows | Return control | levels [m], flows [mL/s] | | | |
| Return control | 3 return pumps (Jabsco 23680-4103) | voltages [V] | | | |

Questions for the lab:
- Sensor conversion (V -> mmHg, L/min, mm): where, coefficients, offset at start-up.
- Position of the pressure sensor "before" each chamber (where P_target must hold).
- Solver of the lab models: fixed/variable, continuous/discrete, step.
- Existing start-up / priming sequence and safety logic.
