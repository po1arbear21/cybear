# grid

## triang_grid
fix get_neighb (exclude self coupling)

## tensor grid
test get_neighb

## test grid
uncomment test region (line 426??)
  aka uncrapify

## submodules
- goal:
  - grid should have grid_table which consists of all idx
    type :: grid
      type(grid_table), alloctable :: full_tab(:,:)
        !! indices: (idx_type, idx_dir)
    end type
- grid_data, grid_table submodules of grid
- vselector
  - vselector_init_var_tab, vselector_init_nvar_tab
  - gtab optional, default: variable%grid%full_tab
- equation
  - equation_depend_variable
  - equation_provide_variable
  - gtab optional, default: variable%grid%full_tab

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

# esystem
## AC+HB+eig
## transient
- transient routines: bwE, bdf2, fwE?
