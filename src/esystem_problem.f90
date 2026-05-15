module esystem_problem_m
  !! Adapter that lets jacobian_validator_m operate on an esystem.
  !!
  !! Wraps a pointer to an esystem and exposes its residual + Jacobian in the
  !! flat (x, f, J) form the validators expect. Use it to check any cybear
  !! Jacobian (poisson, continuity, schottky_bc, dd) against numerical
  !! references: finite differences, Taylor slope, JVP+Richardson.
  !!
  !! Does NOT extend jacobian_problem_cmplx because the DD residual is not
  !! analytic in general (mobility max(), Fermi-Dirac, sign-conditioned
  !! tunneling). The complex-step validator therefore is not applicable.

  use jacobian_validator_m, only: jacobian_problem
  use esystem_m,            only: esystem
  use dense_m,              only: dense_real
  use vselector_m,          only: vselector

  implicit none

  private
  public :: esystem_problem
  public :: print_esystem_layout

  type, extends(jacobian_problem) :: esystem_problem
    type(esystem), pointer :: sys => null()
  contains
    procedure :: eval_real => esp_eval_real
    procedure :: eval_jac  => esp_eval_jac
  end type

contains

  subroutine esp_eval_real(this, x, f)
    class(esystem_problem), intent(in)  :: this
    real,                   intent(in)  :: x(:)
    real,                   intent(out) :: f(:)

    call this%sys%set_x(x)
    call this%sys%eval(f = f)
  end subroutine

  subroutine esp_eval_jac(this, x, J)
    class(esystem_problem), intent(in)  :: this
    real,                   intent(in)  :: x(:)
    real,                   intent(out) :: J(:,:)

    type(dense_real) :: D

    call this%sys%set_x(x)
    call this%sys%eval()
    call this%sys%get_df(D)
    J = D%d
  end subroutine

  subroutine print_esystem_layout(sys, unit, indent)
    !! Print the row layout of an esystem: each block's row range plus the
    !! main-variable selector name and the grid-table name.
    !!
    !! Use to map row/column indices flagged by print_top_defects back to
    !! a physical equation / variable / mesh subdomain. Rows and columns
    !! share the same layout in esystem (residual_i is the equation governing
    !! the i-th flat variable entry).
    type(esystem), intent(in)           :: sys
    integer,       intent(in), optional :: unit
    integer,       intent(in), optional :: indent

    integer                   :: u, ind, ibl, iimvar, itab
    class(vselector), pointer :: mv
    character(:), allocatable :: pad

    u   = 6
    ind = 8
    if (present(unit))   u   = unit
    if (present(indent)) ind = indent
    allocate (character(ind) :: pad)
    pad = repeat(" ", ind)

    write (u, "(A,A)") pad, "esystem layout (rows = cols, main-var blocks):"
    do ibl = 1, sys%nbl
      iimvar =  sys%block2res(1, ibl)
      itab   =  sys%block2res(2, ibl)
      mv     => sys%get_main_var(iimvar)
      write (u, "(A,A,I3,A,I5,A,I5,A,A,I0,A,A,A,A)")        &
        pad, "  block ", ibl,                               &
        ":  rows [", sys%i0(ibl), " .. ", sys%i1(ibl), "]", &
        "  n=",   sys%i1(ibl) - sys%i0(ibl) + 1,            &
        "  ",     mv%name,                                  &
        " / ",    trim(mv%tab(itab)%p%name)
    end do
  end subroutine

end module
