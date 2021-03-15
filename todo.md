# grid
## triang
neighbors
  not include itself (otherwise tensor grid broken)

## edge
get edge length
  procedure: get_len(idx)
  idx: edge index

## tensor grid
get_neighb

## test grid
uncomment test region (line 426??)
  aka uncrapify

# equation
## precon
- add 2nd jacobian "jaco_precon"
- turned off by default
- could use different stencil!
- include in res_equation
  - jaco_f_precon

## test
- test
  - test_jacobian (finite differences; overwrite in res_equation)
  - optional
  - loop dependencies
    - relative step width

# simple_eq
## vsel <-> var + tab
- `call dum_eq%init(var, tab, name)`
  - creates vselector internally
  - generic routine
  - name optional

# variable
## output
- csv
- plt

# AC+transient+HB
