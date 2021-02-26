# equation

## provide/depend
- make equation_provide a function returning iprov
- make equation_depend a function returning idep

## test
- test
  - test_jacobian (finite differences; overwrite in res_equation)
  - optional
  - loop dependencies
    - relative step width
