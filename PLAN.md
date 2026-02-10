  Summary

  Branch: feat/schottky (based on main)

  New file

  - src/schottky.f90 (~230 lines) — Cleaned-up Schottky physics module with:
    - get_normal_dir — contact normal detection from geometry
    - schottky_velocity — Richardson thermionic velocity
    - schottky_n0b — equilibrium injection density with IFBL
    - schottky_tunneling — WKB sub-barrier tunneling (Tsu-Esaki integral)
    - wkb_transmission — transmission coefficient with smooth transition

  Modified files
  ┌───────────────────────┬────────────────────────────────────────────────────────────────────────────┐
  │         File          │                                  Changes                                   │
  ├───────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ src/contact.f90       │ CT_SCHOTTKY=3, Schottky fields, set_phims_schottky                         │
  ├───────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ src/region.f90        │ Schottky parameter parsing from .ini files                                 │
  ├───────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ src/device_params.f90 │ Schottky field propagation + set_phims_schottky call                       │
  ├───────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ src/continuity.f90    │ Robin BC: optional efield param, per-iteration Jacobian, reset+add pattern │
  ├───────────────────────┼────────────────────────────────────────────────────────────────────────────┤
  │ src/device.f90        │ Conditional efield passing, calc_efield in sys_dd                          │
  └───────────────────────┴────────────────────────────────────────────────────────────────────────────┘
  Optimizations vs feat/robin

  - Removed ~30 DEBUG prints, 100+ lines of dead commented code, unused log1p_exp
  - Fixed phi_b - norm(0.2, "eV") debug hack in tunneling
  - Adapted to main's API: get_inv_dist, newton/newton_opt, add()+reset() instead of set()
  - jaco_dens only non-constant when Schottky contacts present

  ★ Insight ─────────────────────────────────────
  The reset()+add() pattern we used instead of feat/robin's set() is actually more explicit about what's happening: each Newton
  iteration, the non-constant Jacobian entries are zeroed, then rebuilt from scratch. This makes the data flow clearer — there's no
  hidden state from a previous iteration leaking through. The set_matr(const=.false., nonconst=.true.) call at the end materializes
  only the non-constant blocks into the sparse matrix, leaving the constant interior stencil untouched.
  ─────────────────────────────────────────────────

