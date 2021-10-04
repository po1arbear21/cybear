# equation
## precon
- add 2nd jacobian "jaco_precon"
- turned off by default
- could use different stencil!
- include in res_equation
  - jaco_f_precon

# esystem
- test

# variable

# analysis
## stationary

## HB

## eig

## transient
- transient routines: bwE, bdf2, fwE?

## Gummel
- init by following arguments
  - esystems: nonlinear PE, transport model (one for each carrier type)
  - equation: dens2imref
    - array, one for each carrier type
    - call eval to get imref
