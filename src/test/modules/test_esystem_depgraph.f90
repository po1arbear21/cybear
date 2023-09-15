m4_include(../../util/macro.f90.inc)

module test_esystem_depgraph_m

  use equation_m,     only: equation
  use error_m,        only: program_error
  use esystem_m,      only: esystem
  use jacobian_m,     only: jacobian
  use grid_m,         only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use grid_data_m,    only: grid_data0_real
  use grid_table_m,   only: grid_table
  use res_equation_m, only: res_equation
  use stencil_m,      only: dirichlet_stencil
  use test_case_m,    only: test_case
  use util_m,         only: fsleep
  use variable_m,     only: variable_real, variable_ptr
  use vector_m,       only: vector_int
  use vselector_m,    only: vselector

  implicit none

  private
  public test_eval_list1, test_eval_list2

  type, extends(variable_real) :: var
  contains
    procedure :: init => var_init
  end type

  type, extends(equation) :: tequation
    integer :: sleep
  contains
    procedure :: init => tequation_init
    procedure :: eval => tequation_eval
  end type

  type, extends(res_equation) :: tres
    integer :: sleep
  contains
    procedure :: init => tres_init
    procedure :: eval => tres_eval
  end type

  logical, parameter :: PARALLEL_PRINT = .false.
    !! For manual purposes

contains

  subroutine test_eval_list1()
    type(test_case) :: tc

    integer             :: ieval, eval_blocks(9)
    type(var)           :: a, b, c, d, e, f, g, h, i, j, k
    type(variable_ptr)  :: empty(0)
    type(tequation)     :: eq1, eq2, eq3, eq4, eq5
    type(tres)          :: res1, res2, res3, res4
    type(esystem)      :: es

    call tc%init("esystem: depgraph evaluation list 1")

    ! init vars
    call a%init("a")
    call b%init("b")
    call c%init("c")
    call d%init("d")
    call e%init("e")
    call f%init("f")
    call g%init("g")
    call h%init("h")
    call i%init("i")
    call j%init("j")
    call k%init("k")

    ! equation system
    call es%init('test eqs system1', parallel_eval = .true.)

    ! equations
    call eq1%init("EQ 1", [e%get_ptr(), f%get_ptr()], [a%get_ptr(), b%get_ptr(), c%get_ptr()], sleep = 3)
    call es%add_equation(eq1)
    call eq2%init("EQ 2", [g%get_ptr()], [a%get_ptr()], sleep = 4)
    call es%add_equation(eq2)
    call eq3%init("EQ 3", [h%get_ptr()], [b%get_ptr(), e%get_ptr()], sleep = 6)
    call es%add_equation(eq3)
    call eq4%init("EQ 4", [i%get_ptr(), j%get_ptr()], [g%get_ptr(), c%get_ptr()], sleep = 8)
    call es%add_equation(eq4)
    call eq5%init("EQ 5", [k%get_ptr()], [d%get_ptr()], sleep = 10)
    call es%add_equation(eq5)

    ! residuals
    call res1%init("RQ 1", a%get_ptr(), [g%get_ptr()], sleep = 3)
    call es%add_equation(res1)
    call res2%init("RQ 2", b%get_ptr(), [e%get_ptr(), i%get_ptr()], sleep = 30)
    call es%add_equation(res2)
    call res3%init("RQ 3", c%get_ptr(), [k%get_ptr(), j%get_ptr()], sleep = 10)
    call es%add_equation(res3)
    call res4%init("RQ 4", d%get_ptr(), empty, sleep = 5)
    call es%add_equation(res4)

    ! finish equation system
    call es%init_final()
    !call es%g%output("src/test/depgraph1", eval_list = .true.)

    ! Expected evaluation list from residuals up
    call tc%assert_eq([2,6,1,4,7,5,8,9], es%g%ieval%to_array(), "ieval list does not match")

    ! Eval blocks of equations
    do ieval = 1, es%g%equs%n
      eval_blocks(ieval) = es%g%equs%d(ieval)%eval_block
    end do
    call tc%assert_eq([1,1,0,2,1,2,3,3,1], eval_blocks, "eval blocks do not match")
    call tc%assert_eq([1,2,5,9,4,6,7,8], es%g%ieval_final, "ieval_final does not match")
    call tc%assert_eq([1,5,7,9], es%g%ieblks, "ieblks list does not match")

    ! parallel execution
    if (PARALLEL_PRINT) then
      call es%eval()
    end if

    ! destruct esystem
    call es%destruct()

    call tc%finish()
  end subroutine

  subroutine test_eval_list2()
    type(test_case) :: tc

    integer             :: ieval, eval_blocks(23), eval_blocks_ex(23), ieval_list(22), feval, ieval0
    type(var)           :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z
    type(variable_ptr)  :: empty(0)
    type(tequation)     :: eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16
    type(tres)          :: res1, res2, res3, res4, res5, res6, res7
    type(esystem)       :: es

    call tc%init("esystem: depgraph evaluation list 2")

    ! init vars
    call a%init("a")
    call b%init("b")
    call c%init("c")
    call d%init("d")
    call e%init("e")
    call f%init("f")
    call g%init("g")
    call h%init("h")
    call i%init("i")
    call j%init("j")
    call k%init("k")
    call l%init("l")
    call m%init("m")
    call n%init("n")
    call o%init("o")
    call p%init("p")
    call q%init("q")
    call r%init("r")
    call s%init("s")
    call t%init("t")
    call u%init("u")
    call v%init("v")
    call w%init("w")
    call x%init("x")
    call y%init("y")
    call z%init("z")

    ! equation system
    call es%init('test eqs system2', parallel_eval = .true.)

    ! equations
    call eq1%init("EQ 1", [h%get_ptr()], [d%get_ptr()], sleep = 10)
    call es%add_equation(eq1)
    call eq2%init("EQ 2", [i%get_ptr()], [a%get_ptr(), b%get_ptr(), c%get_ptr(), h%get_ptr()], sleep = 1)
    call es%add_equation(eq2)
    call eq3%init("EQ 3", [j%get_ptr()], [a%get_ptr(), i%get_ptr()], sleep = 5)
    call es%add_equation(eq3)
    call eq4%init("EQ 4", [k%get_ptr()], [b%get_ptr(), e%get_ptr(), j%get_ptr()], sleep = 20)
    call es%add_equation(eq4)
    call eq5%init("EQ 5", [l%get_ptr()], [g%get_ptr(), j%get_ptr(), k%get_ptr()], sleep = 12)
    call es%add_equation(eq5)
    call eq6%init("EQ 6", [m%get_ptr()], [d%get_ptr()], sleep = 4)
    call es%add_equation(eq6)
    call eq7%init("EQ 7", [n%get_ptr()], [l%get_ptr(), b%get_ptr()], sleep = 3)
    call es%add_equation(eq7)
    call eq8%init("EQ 8", [o%get_ptr()], [a%get_ptr()], sleep = 1)
    call es%add_equation(eq8)
    call eq9%init("EQ 9", [p%get_ptr(), q%get_ptr()], [b%get_ptr(), e%get_ptr()], sleep = 3)
    call es%add_equation(eq9)
    call eq10%init("EQ 10", [r%get_ptr()], [g%get_ptr(), l%get_ptr()], sleep = 30)
    call es%add_equation(eq10)
    call eq11%init("EQ 11", [s%get_ptr()], [i%get_ptr()], sleep = 7)
    call es%add_equation(eq11)
    call eq12%init("EQ 12", [t%get_ptr()], [a%get_ptr(), p%get_ptr(), d%get_ptr()], sleep = 5)
    call es%add_equation(eq12)
    call eq13%init("EQ 13", [u%get_ptr()], [a%get_ptr(), r%get_ptr(), c%get_ptr()], sleep = 5)
    call es%add_equation(eq13)
    call eq14%init("EQ 14", [v%get_ptr(), w%get_ptr()], [t%get_ptr(), e%get_ptr()], sleep = 6)
    call es%add_equation(eq14)
    call eq15%init("EQ 15", [x%get_ptr(), y%get_ptr()], [v%get_ptr(), n%get_ptr()], sleep = 2)
    call es%add_equation(eq15)
    call eq16%init("EQ 16", [z%get_ptr()], [t%get_ptr(), x%get_ptr(), s%get_ptr()], sleep = 5)
    call es%add_equation(eq16)

    ! residuals
    call res1%init("RQ 1", a%get_ptr(), [t%get_ptr()], sleep = 2)
    call es%add_equation(res1)
    call res2%init("RQ 2", b%get_ptr(), [a%get_ptr(), r%get_ptr()], sleep = 1)
    call es%add_equation(res2)
    call res3%init("RQ 3", c%get_ptr(), [t%get_ptr(), j%get_ptr()], sleep = 20)
    call es%add_equation(res3)
    call res4%init("RQ 4", d%get_ptr(), empty, sleep = 1)
    call es%add_equation(res4)
    call res5%init("RQ 5", e%get_ptr(), [l%get_ptr(), h%get_ptr(), v%get_ptr()], sleep = 40)
    call es%add_equation(res5)
    call res6%init("RQ 6", f%get_ptr(), [b%get_ptr(), x%get_ptr(), z%get_ptr()], sleep = 10)
    call es%add_equation(res6)
    call res7%init("RQ 7", g%get_ptr(), [m%get_ptr(), y%get_ptr(), o%get_ptr(), e%get_ptr(), f%get_ptr()], sleep = 10)
    call es%add_equation(res7)

    ! finish equation system
    call es%init_final()
    !call es%g%output("src/test/depgraph2", eval_list = .true.)

    ! Expected evaluation list from residuals up
    ieval_list = [9,12,17,1,2,3,4,5,10,18,19,20,14,21,7,15,11,16,22,6,8,23]
    call tc%assert_eq(ieval_list, es%g%ieval%to_array(), "ieval list does not match")

    ! Eval blocks of equations
    do ieval = 1, es%g%equs%n
      eval_blocks(ieval) = es%g%equs%d(ieval)%eval_block
    end do
    eval_blocks_ex = [1,2,3,4,5,1,6,1,1,6,3,2,0,3,7,8,3,7,4,1,6,9,8]
    call tc%assert_eq(eval_blocks_ex, eval_blocks, "eval blocks do not match")
    ! Unpredictable qsort?
    ! feval_list = [1,6,8,9,20,2,12,3,11,14,17,4,19,5,7,10,21,15,18,16,23,22]
    call tc%assert_eq([1,6,8,12,14,15,18,20,22,23], es%g%ieblks, "ieblks list does not match")

    ! Check no unexpected eval blocks in eval list
    do ieval = 1, es%g%neblks
      ieval0 = es%g%ieval_final(es%g%ieblks(ieval))
      do feval = es%g%ieblks(ieval)+1, es%g%ieblks(ieval+1) - 1
        call tc%assert_eq(es%g%equs%d(ieval0)%eval_block, es%g%equs%d(es%g%ieval_final(feval))%eval_block, "Final evaluation order not as expected")
      end do
    end do

    ! parallel execution
    if (PARALLEL_PRINT) then
      call es%eval()
    end if

    ! destruct esystem
    call es%destruct()

    call tc%finish()
  end subroutine

  subroutine var_init(this, name)
    class(var),   intent(out) :: this
    character(*), intent(in)  :: name

    call this%variable_init(name, '1')
  end subroutine

  subroutine tequation_init(this, name, prov, dep, sleep)
    class(tequation), target, intent(out) :: this
    character(*),             intent(in)  :: name
    type(variable_ptr),       intent(in)  :: prov(:)
      !! provided vars
    type(variable_ptr),       intent(in)  :: dep(:)
      !! dependent vars
    integer, optional,   intent(in)  :: sleep

    integer                 :: i, j
    integer, allocatable    :: iprov(:), idep(:)
    type(jacobian), pointer :: jaco

    ! init base
    call this%equation_init(name)

    ! provides
    allocate(iprov(size(prov)))
    do i = 1, size(prov)
      iprov(i) = this%provide(prov(i)%p)
    end do

    ! add dependency
    allocate(idep(size(dep)))
    do i = 1, size(dep)
      idep(i) = this%depend(dep(i)%p)
    end do

    ! init jacobians
    do i = 1, size(prov)
      do j = 1, size(dep)
        jaco => this%init_jaco(iprov(i), idep(j))
      end do
    end do

    ! finish initialization
    call this%init_final()

    this%sleep = 0
    if (present(sleep)) this%sleep = sleep

    deallocate(iprov)
    deallocate(idep)
  end subroutine

  subroutine tequation_eval(this)
    class(tequation), intent(inout) :: this

    if (this%sleep /= 0 .and. PARALLEL_PRINT) then
      print *, "┌──── ", this%name
      call fsleep(this%sleep)
      print *, "└──── ", this%name
    end if
  end subroutine

  subroutine tres_init(this, name, mvar, dep, sleep)
    class(tres), target, intent(out) :: this
    character(*),        intent(in)  :: name
    type(variable_ptr),  intent(in)  :: mvar
      !! main var
    type(variable_ptr),  intent(in)  :: dep(:)
      !! dependent vars
    integer, optional,   intent(in)  :: sleep

    integer                 :: i, idep
    type(jacobian), pointer :: jaco

    ! init base
    call this%equation_init(name)

    ! main var
    call this%init_f(mvar%p)

    ! dependencies
    idep =  this%depend(mvar%p)
    jaco => this%init_jaco_f(idep)
    do i = 1, size(dep)
      idep =  this%depend(dep(i)%p)
      jaco => this%init_jaco_f(idep)
    end do

    ! finish initialization
    call this%init_final()

    this%sleep = 0
    if (present(sleep)) this%sleep = sleep
  end subroutine

  subroutine tres_eval(this)
    class(tres), intent(inout) :: this

    if (this%sleep /= 0 .and. PARALLEL_PRINT) then
      print *, "┌──── ", this%name
      call fsleep(this%sleep)
      print *, "└──── ", this%name
    end if
  end subroutine
end module
