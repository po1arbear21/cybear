# Cybear

A general-purpose drift–diffusion semiconductor device simulator, forked from
the ITHE "Drift–Diffusion" codebase and extended with advanced Schottky contact
physics, heterojunction support, trap-assisted transport, beam / EBIC analysis,
adjoint-based numerics, and NEGF-coupled boundary conditions.

## Build

Uses [fargo](https://git.rwth-aachen.de/ithe/fargo).

```bash
fargo build release
fargo run test
```

`run.ini` holds the default solver parameters.

## Repository layout

```
src/                Fortran sources (Poisson, continuity, contacts, …)
lib/fortran-basic/  Core numerics submodule
devices/            Device configurations
  examples/         Reference geometries
  galene/           GALENE-format geometry files
test/               Test programs and fixtures
fargo.toml          Build configuration
run.ini             Default solver parameters
tables/             Distribution lookup tables
```

## Physics

- Drift–diffusion transport
- Self-consistent Poisson + continuity
- Contact types: Ohmic, Schottky, RealOhmic, Gate
- Steady-state, transient, small-signal AC, and harmonic-balance analyses
- Schottky contact physics: thermionic emission, tunneling, image-force barrier lowering
- Trap-assisted transport
