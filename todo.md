# equation
## precon
- add 2nd jacobian "jaco_precon"
- turned off by default
- could use different stencil!
- include in res_equation
  - jaco_f_precon

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

## dag
- rewrite from scratch
- rename dag to depgraph; dag_node => node; dag_equ => equ
- use integer indices instead of pointers
- use hash of vselector for searching when building dependency graph
