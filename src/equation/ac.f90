module ac_m

  use esystem_m, only: esystem
  use matrix_m,  only: sparse_cmplx, matrix_add

  implicit none

  type ac
    !! small-signal analyzer

    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system

    integer, allocatable :: isrc(:)
      !! flat indices of cosine sources (e.g. contact voltages)

    complex, allocatable :: s(:)
      !! complex 'frequency': x = x_0 + Re{x_1 * exp(s*t)}

    complex, allocatable :: x(:,:,:)
      !! small-signal data (sys%n, size(isrc), size(s))
    complex, allocatable :: dxds(:,:,:)
      !! derivatives of x wrt s (sys%n, size(isrc), size(s))
  contains
    procedure :: run_analysis => ac_run_analysis
    procedure :: select_real  => ac_select_real
    procedure :: select_imag  => ac_select_imag
    procedure :: select_abs   => ac_select_abs
    procedure :: select_phase => ac_select_phase
  end type

contains

  subroutine ac_run_analysis(this, sys, isrc, s, calc_dxds)
    !! perform small signal analysis
    class(ac),             intent(out)   :: this
    type(esystem), target, intent(inout) :: sys
      !! equation system
    integer,               intent(in)    :: isrc(:)
      !! flat indices indicating where to place cosine sources
      !! make sure corresponding row in jacobian is 1 on main diagonal, 0 elsewhere
      !! sys%dft(isrc,:) must be 0
    complex,               intent(in)    :: s(:)
      !! assume: x(isrc) = x_0(isrc) + Re{exp(s*t)}
    logical, optional,     intent(in)    :: calc_dxds
      !! calculate dxds in addition to x (default: false)

    integer              :: i, nsrc, ns
    logical              :: calc_dxds_
    complex, allocatable :: rhs(:,:), x(:,:), tmp(:,:), dxds(:,:)
    type(sparse_cmplx)   :: df, dft, mat

    ! optional arguments
    calc_dxds_ = .false.
    if (present(calc_dxds)) calc_dxds_ = calc_dxds

    ! dimensions
    nsrc = size(isrc)
    ns   = size(s)

    ! set simple members
    this%sys  => sys
    this%isrc =  isrc
    this%s    =  s

    ! allocate memory
    allocate (this%x(sys%n,nsrc,ns))
    if (calc_dxds_) then
      allocate (this%dxds(sys%n,nsrc,ns))
    end if

    ! allocate temporary memory
    call mat%init(df%nrows)
    allocate (rhs(df%nrows,nsrc), source = (0.0, 0.0))
    allocate (  x(df%nrows,nsrc), source = (0.0, 0.0))
    if (calc_dxds_) then
      allocate ( tmp(df%nrows,nsrc), source = (0.0, 0.0))
      allocate (dxds(df%nrows,nsrc), source = (0.0, 0.0))
    end if

    ! evaluate equation system
    call sys%eval()

    ! get jacobian and jacobian wrt time derivative
    call sys%get_df(df)
    call sys%get_dft(dft)

    ! set right-hand sides
    do i = 1, nsrc
      rhs(isrc(i),i) = 1.0
    end do

    ! get admittance parameters for all frequencies
    do i = 1, ns
      print "(I6,A,I6)", i, "/", ns

      ! construct matrix (mat = df + s * dft)
      call matrix_add(df, dft, mat, fact2 = s(i))

      ! factorize and solve
      call mat%factorize()
      call mat%solve_mat(rhs, x)

      ! extract small-signal data
      this%x(:,:,i) = x(:,:)

      ! derivatives
      if (calc_dxds_) then
        ! tmp = dft * x
        call dft%mul_mat(x, tmp)

        ! dxds = mat \ (- dft * x)
        call mat%solve_mat(-tmp, dxds)

        ! extract derivatives
        this%dxds(:,:,i) = dxds(:,:)
      end if

      ! release memory
      call mat%destruct()
    end do
  end subroutine

  subroutine ac_select_real(this, isrc, is)
    !! save real part of small-signal data (overwrite equation system variables)
    class(ac), intent(in) :: this
    integer,   intent(in) :: isrc
      !! source index
    integer,   intent(in) :: is
      !! frequency index

    call this%sys%set_x(real(this%x(:,isrc,is)))
  end subroutine

  subroutine ac_select_imag(this, isrc, is)
    !! save imaginary part of small-signal data (overwrite equation system variables)
    class(ac), intent(in) :: this
    integer,   intent(in) :: isrc
      !! source index
    integer,   intent(in) :: is
      !! frequency index

    call this%sys%set_x(aimag(this%x(:,isrc,is)))
  end subroutine

  subroutine ac_select_abs(this, isrc, is)
    !! save absolute value of small-signal data (overwrite equation system variables)
    class(ac), intent(in) :: this
    integer,   intent(in) :: isrc
      !! source index
    integer,   intent(in) :: is
      !! frequency index

    call this%sys%set_x(abs(this%x(:,isrc,is)))
  end subroutine

  subroutine ac_select_phase(this, isrc, is)
    !! save phase of small-signal data (overwrite equation system variables)
    class(ac), intent(in) :: this
    integer,   intent(in) :: isrc
      !! source index
    integer,   intent(in) :: is
      !! frequency index

    call this%sys%set_x(atan2(aimag(this%x(:,isrc,is)),real(this%x(:,isrc,is))))
  end subroutine

end module
