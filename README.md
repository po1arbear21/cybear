# Cybear — General-Purpose Drift-Diffusion Semiconductor Device Simulator

> *A comprehensive solver for carrier transport in advanced semiconductor devices, forked from the ITHE "Drift-Diffusion" codebase and extensively refactored for novel device architectures.*

## Overview

Cybear is a versatile drift-diffusion (DD) simulator capable of modeling a wide range of semiconductor devices including:

- **Reconfigurable Nanowire FETs** - Quantum-confined channel devices with dynamic polarity control
- **2D Material MOSFETs** - Transition metal dichalcogenide and graphene-based transistors
- **Perovskite Vertical FETs** - Novel architecture with perforated source electrodes
- **Schottky Barrier Devices** - Metal-semiconductor junctions with advanced tunneling models
- **Traditional MOSFETs** - Silicon and compound semiconductor devices

## Key Features

### Physics Models
- Drift-diffusion transport with field-dependent mobility
- Poisson equation for electrostatics
- Thermionic emission and thermionic-field emission at Schottky contacts
- Image force barrier lowering (IFBL)
- Incomplete ionization models
- Generation-recombination processes
- Quantum confinement effects (in development)

### Numerical Methods
- Finite volume discretization
- Newton-Raphson solver with adaptive damping
- Gummel iteration for decoupled solution
- Full Newton for coupled system
- Ramo-Shockley current calculation

### Device Support
- 1D, 2D, and 3D geometries
- Multiple contact types (Ohmic, Gate, Schottky)
- Heterojunctions and material interfaces
- Time-dependent and steady-state analysis
- Small-signal AC analysis

## Build System

Built with Fargo - a modern Fortran build tool designed for scientific computing.

## Author

Chenyang Yu, PhD student at RWTH Aachen University

## License

See LICENSE file for details.

Last updated: 2025-10-06