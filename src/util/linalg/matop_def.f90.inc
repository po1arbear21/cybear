type, abstract :: matop_{}T
  !! Base type for matrix operations on vectors

  integer :: nrows
    !! Number of rows (size of output vector)
  integer :: ncols
    !! Number of columns (size of input vector)
contains
  procedure :: matop_{}T{}_init
  generic   :: init  => matop_{}T{}_init
  procedure :: exec1 => matop_{}T{}_exec1
  procedure :: exec2 => matop_{}T{}_exec2
  generic   :: exec  => exec1, exec2
end type

type, extends(matop_{}T) :: single_matop_{}T
  !! matrix operation using a single matrix

  class(matrix_{}T), pointer :: mat => null()
    !! pointer to matrix
  logical                    :: inv
    !! inverse flag
contains
  procedure :: single_matop_{}T{}_init
  generic   :: init  => single_matop_{}T{}_init
  procedure :: exec1 => single_matop_{}T{}_exec1
  procedure :: exec2 => single_matop_{}T{}_exec2
end type

type, extends(matop_{}T) :: chain_matop_{}T
  !! multiple matrices chainged together

  type(matrix_ptr_{}T), allocatable :: mat(:)
    !! pointers to matrices
  logical,              allocatable :: inv(:)
    !! inverse flags
contains
  procedure :: chain_matop_{}T{}_init
  generic   :: init    => chain_matop_{}T{}_init
  procedure :: set_mat => chain_matop_{}T{}_set_mat
  procedure :: exec1   => chain_matop_{}T{}_exec1
  procedure :: exec2   => chain_matop_{}T{}_exec2
end type
