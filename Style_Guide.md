# Style Guide

A quick guide of coding style conventions used in this project.

## Beginning of a program
### Importing other modules
- Import only modules you need and include out of these just the functions you are using in the program
- Order the modules and the functions preferably in an alphabetic order (or if more appropiate in the order of importance)
- Align the *use* and *only* commands each

Example:
```f90
use example_matrices_m, only: matrix2
use gmres_m,            only: gmres
use matop_m,            only: single_matop_real
use matrix_m,           only: sparse_real
use mkl_ilu_m,          only: mkl_ilu0, mkl_ilut
use test_case_m,        only: test_case
use util_m,             only: int2str
```
### Implicit types
- Turn the implicit type option off
```f90
implicit none
```
## Midprogram
### Functions and subroutines
- Specify the type of the function if needed
- Append the input of the function/ subroutine in a reasonable order
- Specify in the same line the output for the function
- Describe the function briefly in a comment below (use two exclamation marks for the documentation)
- Contrary to statements (e.g. if-statements) do not use a space between the functions name and the brackets

Example:
```f90
pure function linspace(x0, x1, nx) result(x)
  !! create array of linear spaced values
```

### Declaring variables
- Declare only variable you actually use in the program
- Order the data types in an alphabetic order (or if more appropiate in the order of importance)
- Order the variables within a declaration in an alphabetic order (or if more appropiate in the order of importance), too
- Align the colons

Example:
```f90
integer           :: n(4), i, j
real              :: a(4), b(4)
real, allocatable :: list(:), dif(:), dif_exp(:)
```
#### Declaring variables within functions
- Declare the in- and output variables individually at the beginning
- Append to each declaration the intention of the variables
- Describe each variable of the in- and output briefly in a comment below the declaration (use two exclamation marks for the documentation)
- Append in a seperate block the declaration of local variables with all the rules described above

Example:
```f90
real,               intent(in)  :: p(:)
  !! function parameters (can be empty array if not needed)
type(newton1D_opt), intent(in)  :: opt
  !! iteration options
real,               intent(in)  :: x0
  !! first guess for solution
real,               intent(out) :: x
  !! output solution
integer, optional,  intent(in)  :: ipar(:)
  !! integer function parameters
real,    optional,  intent(out) :: dxdp(:)
  !! optional output derivative of x wrt parameters

! local variables
integer :: it
real    :: err, err0, xmin, xmax, fmin, fmax
real    :: f, dfdx, dfdp(size(p)), dx
logical :: bounded
```
### Assigning values to variables
- Use in the context reasonable blocks for the assignments of variables (e.g. use a block for variables of the same datatype or with the same intent)
- Align the equal signs within a block
- For assigning arrays align all entries:
  * for integers align the last digits
  * for floats align the decimal point

Example:
```f90
a = [  0.0,  10.0,  431.23140,   0.0]
b = [100.0, -10.0, 1233.42576,  1e-5]

e0 = log(x0)
e1 = log(x1)
```
### Using mathematical operations
- Use between every variable and operator a space

Example:
```f90
a = (b * c) / (d * e)**2
```
