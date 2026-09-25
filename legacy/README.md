# Legacy scripts (read-only)

These are the original case folders, unmodified, kept for traceability. They are **not** on the MATLAB path; `setup_paths` excludes this directory.

| Folder | Original name | Content |
|---|---|---|
| `ieee14_pv_original` | `continuation-power-flow-master with pv` | IEEE-14 CPF + PV |
| `ieee14_wind_original` | `continuation-power-flow-master with wind` | IEEE-14 CPF + wind |
| `ieee14_hybrid_original` | `Trabsmission System` | IEEE-14 CPF + PV + wind, L-index (`ind.m`) |
| `ieee33_pv_original` | `mod clf for bibc pv` | IEEE-33 BIBC/BFS + PV |
| `ieee33_wind_original` | `mod clf for bibc wind` | IEEE-33 BIBC/BFS + wind |
| `ieee33_hybrid_original` | `Distribution System` | IEEE-33 BFS + PV + wind, FVSI (`indccds.m`) |

See `../docs/BUG_REPORT.md` for what was wrong in each one. Delete this folder once you no longer need the comparison.
