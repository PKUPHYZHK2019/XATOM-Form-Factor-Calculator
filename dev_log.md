================================================================================
Benchmark
================================================================================

## 2026-06-16 — [Benchmark] `calculate_formfac_baseline`

### Goal
Establish baseline performance before modifying `calculate_formfac()`.

### Benchmark setup
- Photon energy: `690 eV`
- Number of configs tested: `3`
- Successful runs: `3`
- Failed runs: `0`
- Python version: `3.12.12`
- Platform: `Linux-4.18.0-553.123.1.el8_10.x86_64-x86_64-with-glibc2.28`

### Runtime summary
- Total runtime: `70.987 s`
- Mean runtime per successful config: `23.662 s`
- Min runtime: `19.779 s`
- Max runtime: `26.820 s`

### Per-config results

| Config index | Success | Runtime / s | f0 | f1 | f2 |
|---:|:---:|---:|---:|---:|---:|
| 0 | True | 19.779 | 27 | -11.938194 | -0.00011698452 |
| 1 | True | 24.387 | 28 | -12.977465 | -0.27146792 |
| 2 | True | 26.820 | 29 | -14.439758 | -0.50818083 |

### Decision
Baseline recorded. Future code changes should be compared against this benchmark.

### Output files
- CSV: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_155636.csv`
- JSON: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_155636.json`
- Markdown: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_155636.md`



2026-06-16 — [Performance] Move environment activation outside run_xatom()

Motivation

run_xatom() was repeatedly activating the Anaconda/XATOM environment for every single call. Since calculate_formfac() may call run_xatom() many times, this introduced unnecessary overhead.

Code change

Removed environment activation from inside run_xatom(). The required Anaconda and XATOM environment is now activated once manually in the terminal before running the Python script.


## 2026-06-16 — [Benchmark] `calculate_formfac_baseline`

### Goal
Establish baseline performance before modifying `calculate_formfac()`.

### Benchmark setup
- Photon energy: `690 eV`
- Number of configs tested: `3`
- Successful runs: `3`
- Failed runs: `0`
- Python version: `3.12.12`
- Platform: `Linux-4.18.0-553.123.1.el8_10.x86_64-x86_64-with-glibc2.28`

### Runtime summary
- Total runtime: `48.561 s`
- Mean runtime per successful config: `16.187 s`
- Min runtime: `13.510 s`
- Max runtime: `18.100 s`

### Per-config results

| Config index | Success | Runtime / s | f0 | f1 | f2 |
|---:|:---:|---:|---:|---:|---:|
| 0 | True | 13.510 | 27 | -11.938194 | -0.00011698452 |
| 1 | True | 16.951 | 28 | -12.977465 | -0.27146792 |
| 2 | True | 18.100 | 29 | -14.439758 | -0.50818083 |

### Decision
Baseline recorded. Future code changes should be compared against this benchmark.

### Output files
- CSV: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_160209.csv`
- JSON: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_160209.json`
- Markdown: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_160209.md`



SQL cache for photon-energy-independent XATOM calculations

Implemented an SQLite-based cache layer for XATOM calculations whose inputs are independent of the evaluated photon energy. The cached quantities include bound-bound transition search results, excited-state lifetimes, and non-resonant photoionization-derived (f’’(E)) tables.

Previously, these XATOM calculations were repeated for each photon energy, even when the same electronic configuration had already been evaluated. The new logic uses the electronic configuration as the database key and checks the SQL cache before launching XATOM. If the configuration already exists in the database, the stored result is loaded directly and the expensive XATOM call is skipped.

SQLite is used in WAL mode with a lock table to support safe parallel execution over many photon energies. Each process opens its own database connection, and the lock mechanism prevents multiple processes from writing the same configuration result at the same time.

Result: for configurations already saved in the SQL cache, the recalculation cost becomes negligible compared with a full XATOM run. This significantly accelerates repeated scans and parallel photon-energy evaluations, especially when the same configurations appear multiple times.



## 2026-06-16 — [Benchmark] `calculate_formfac_baseline`

### Goal
Establish baseline performance before modifying `calculate_formfac()`.

### Benchmark setup
- Photon energy: `690 eV`
- Number of configs tested: `3`
- Successful runs: `3`
- Failed runs: `0`
- Python version: `3.12.12`
- Platform: `Linux-4.18.0-553.123.1.el8_10.x86_64-x86_64-with-glibc2.28`

### Runtime summary
- Total runtime: `0.032 s`
- Mean runtime per successful config: `0.011 s`
- Min runtime: `0.005 s`
- Max runtime: `0.022 s`

### Per-config results

| Config index | Success | Runtime / s | f0 | f1 | f2 |
|---:|:---:|---:|---:|---:|---:|
| 0 | True | 0.022 | 27 | -11.938194 | -0.00011698452 |
| 1 | True | 0.006 | 28 | -12.977465 | -0.27146792 |
| 2 | True | 0.005 | 29 | -14.439758 | -0.50818083 |

### Decision
Baseline recorded. Future code changes should be compared against this benchmark.

### Output files
- CSV: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_175529.csv`
- JSON: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_175529.json`
- Markdown: `/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles/benchmarks/calculate_formfac_baseline_20260616_175529.md`