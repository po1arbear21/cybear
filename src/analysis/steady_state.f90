m4_include(../util/macro.f90.inc)

module steady_state_m

  use array_m,       only: array1_real, array1_cmplx
  use error_m,       only: assert_failed, program_error
  use esystem_m,     only: esystem
  use gmres_m,       only: gmres_options
  use grid_data_m,   only: grid_data_real, grid_data_cmplx
  use input_src_m,   only: input_src
  use json_m,        only: json_array, json_object
  use matrix_m,      only: matrix_real, sparse_real
  use newton_m,      only: newton, newton_opt
  use output_file_m, only: output_file
  use variable_m,    only: variable_real, variable_cmplx

  implicit none

  private
  public steady_state

  type steady_state
    !! steady-state analyzer

    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system
    real,      allocatable :: x(:,:)
      !! steady-state data
  contains
    procedure :: run    => steady_state_run
    procedure :: select => steady_state_select
  end type

  abstract interface
    subroutine gummel()
    end subroutine
  end interface

contains

  recursive subroutine steady_state_run(this, sys, nopt, gopt, input, t_input, gum, ofile, oname)
    !! perform steady-state analysis for one or multiple sets of input parameters
    class(steady_state),           intent(out)   :: this
    type(esystem),       target,   intent(inout) :: sys
      !! equation system
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)
    class(input_src),    optional, intent(in)    :: input
      !! optional input source
    real,                optional, intent(in)    :: t_input(:)
      !! perform quasi-stationary simulations for these time-points (default: [0.0])
    procedure(gummel),   optional                :: gum
      !! optional procedure to run before newton iteration
    type(output_file),   optional, intent(inout) :: ofile
      !! output file
    character(*),        optional, intent(in)    :: oname
      !! output name

    integer                         :: nt, i, j, i0, i1
    integer,            allocatable :: sh(:)
    real                            :: p(0)
    real,               allocatable :: t(:)
    type(array1_real),  allocatable :: rx(:)
    type(array1_cmplx), allocatable :: cx(:)
    type(newton_opt)                :: nopt_
    type(sparse_real),  target      :: df, dfp
    type(json_object),  pointer     :: obj, data

    ! save pointer to equation system for later use in select routine
    this%sys => sys

    ! parse optional newton options
    if (present(nopt)) then
      m4_assert(size(nopt%atol) == sys%n)
      nopt_ = nopt
    else
      call nopt_%init(sys%n)
    end if

    ! input
    if (present(input)) then
      m4_assert(input%n == sys%ninput)
      if (present(t_input)) then
        t = t_input
      else
        t = [0.0]
      end if
      nt = size(t)
    else
      m4_assert(sys%ninput == 0)
      m4_assert(.not. present(t_input))
      nt = 1
    end if

    ! allocate memory to store results
    allocate (this%x(sys%n, nt))

    ! output initialization
    if (present(ofile)) then
      m4_assert(present(oname))
      obj => ofile%new_object("SteadyState")
      call obj%add_string("Name", oname)
      if (nt > 1) then
        call ofile%write(obj, "Time", t)
      end if
      call obj%add_object("Data", p = data)

      allocate (rx(sys%vars%n))
      allocate (cx(sys%vars%n))
      do i = 1, sys%vars%n
        select type (v => sys%vars%d(i)%p)
        class is (variable_real)
          allocate (rx(i)%d(v%data%n * nt))
        class is (variable_cmplx)
          allocate (cx(i)%d(v%data%n * nt))
        end select
      end do
    end if

    ! solve steady-state for each input time
    do i = 1, nt
      if (present(input)) call sys%set_input(input%get(t(i)))
      if (present(gum))   call gum()
      call newton(fun, p, nopt_, sys%get_x(), this%x(:,i), gmres_opt = gopt)

      ! release factorization
      if (nopt_%it_solver) then
        call dfp%destruct()
      else
        call df%destruct()
      end if

      ! collect values
      if (present(ofile)) then
        do j = 1, sys%vars%n
          select type (v => sys%vars%d(j)%p)
          class is (variable_real)
            i0 = v%data%n * (i-1) + 1
            i1 = v%data%n * i
            rx(j)%d(i0:i1) = v%data%get()
          class is (variable_cmplx)
            i0 = v%data%n * (i-1) + 1
            i1 = v%data%n * i
            cx(j)%d(i0:i1) = v%data%get()
          end select
        end do
      end if
    end do

    ! output
    if (present(ofile)) then
      do i = 1, sys%vars%n
        select type (v => sys%vars%d(j)%p)
        class is (variable_real)
          if (nt > 1) then
            allocate (sh(size(v%data%idx_n)+1))
            sh(size(v%data%idx_n)+1) = nt
          else
            allocate (sh(size(v%data%idx_n)))
          end if
          sh(1:size(v%data%idx_n)) = v%data%idx_n

          call ofile%write(data, v%name, rx(i)%d, sh, unit = v%unit)

          deallocate (sh)
        class is (variable_cmplx)
          if (nt > 1) then
            allocate (sh(size(v%data%idx_n)+1))
            sh(size(v%data%idx_n)+1) = nt
          else
            allocate (sh(size(v%data%idx_n)))
          end if
          sh(1:size(v%data%idx_n)) = v%data%idx_n

          call ofile%write(data, v%name, cx(i)%d, sh, unit = v%unit)

          deallocate (sh)
        end select
      end do
    end if

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_ignore(p)
      m4_ignore(dfdp)

      ! set variables
      call this%sys%set_x(x)

      ! evaluate system to compute residuals, jacobian, and possible preconditioner
      if (nopt_%it_solver) then
        call this%sys%eval(f = f, df = df, dfp = dfp)
      else
        call this%sys%eval(f = f, df = df)
      end if

      ! output jacobian pointers
      if (present(dfdx     )) dfdx      => df
      if (present(dfdx_prec)) dfdx_prec => dfp
    end subroutine

  end subroutine

  subroutine steady_state_select(this, i)
    !! save steady-state results (overwrite esystem variables)
    class(steady_state), intent(in) :: this
    integer,             intent(in) :: i
      !! input set index

    call this%sys%set_x(this%x(:,i))
  end subroutine

end module
