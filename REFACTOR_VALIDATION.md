# Phase 2 modular refactor — validation required

This source tree is a **candidate** derived from the user's MATLAB Online zip that previously passed 13/13 numerical experiments **before** this refactor.

Changes:
- `phase2_main.m` is a function that returns explicit baseline state to group functions.
- 15 experiment entry points are functions; `run_phase2` uses `feval` instead of workspace `eval`.
- Group and test calls to old processing scripts now use explicit `_process` functions.
- `simulate_adc_cube_process` additionally returns the same generated complex noise.
- `roi_only_cfar_process` additionally returns detected range/doppler indices and sparse range lookup, without changing the CFAR loop.
- Original test equations, RNG seeds, CFAR parameters, margin rules and numerical thresholds were retained.
- Historical, unreferenced processing scripts remain **intentionally** until MATLAB regression passes on this candidate.

In MATLAB Online, upload / replace with this tree, then run:

```matlab
cd('/MATLAB Drive/matlab_fmcw_pointcloud')
[unused, blockers] = check_legacy_dependencies();
all_pass = run_phase2_regression_suite();
assert(unused && all_pass, "Dependency or MATLAB regression failed");
```

The suite runs base, 2, 3A through 5E (16 cases) and prints errors. Group 5C may reuse an existing checkpoint. Reconfirm output CSV numbers versus the pre-refactor outputs, especially Group 3B/3D Monte Carlo seeded results, Group 4B, and Group 5E sparse CFAR accuracy.

Do not delete legacy source files or merge into `main` until the **candidate** passes regression on MATLAB Online. This offline preparation has NOT run MATLAB. PNG export must be verified separately where relevant.
