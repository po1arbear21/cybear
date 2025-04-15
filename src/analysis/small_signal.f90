m4_include(../util/macro.f90.inc)

module small_signal_m

  use container_m,     only: container, DYNAMIC_EXT, STORAGE_WRITE, DYNAMIC_APP, ERR_ALREADY_EXISTS
  use error_m,         only: assert_failed, program_error
  use esystem_m,       only: esystem
  use gmres_m,         only: gmres_options, gmres
  use grid_data_m,     only: allocate_grid_data0_cmplx, grid_data_cmplx
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use matrix_m,        only: block_cmplx, matrix_add, matrix_convert, SPSOLVER_PARDISO, SOLVER_GMRES
  use matop_m,         only: single_matop_cmplx
  use string_m,        only: string, new_string
  use util_m,          only: int2str
  use variable_m,      only: variable_ptr

  implicit none

  private
  public small_signal

  type small_signal
    !! small-signal simulation -> solver linear equation (s*M + dfdx(x0)) * dx1/dInput = -df/dInput

    ! general data
    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system
    complex, allocatable   :: s(:)
      !! complex 'frequency': x = x_0 + Re{x_1 * exp(s*t)} with s = sigma + j omega
    complex, allocatable   :: x(:,:,:)
      !! small-signal data (sys%n, size(isrc), size(s)), only used if results are stored in RAM memory
    complex, allocatable   :: dxds(:,:,:)
      !! derivatives of x wrt s (sys%n, size(isrc), size(s))
    logical                :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical                :: log
      !! enable/disable logging (default false)

    ! solver parameters
    integer             :: solver
      !! matrix solver (default SPSOLVER_PARDISO)
    type(gmres_options) :: gopt
      !! options for gmres in case an iterative solver is used (default atol 1e-16, default rtol 1e-12)

    ! data for output (default no output)
    type(string), allocatable :: output_vars(:)
      !! list of variables for output
    character(:), allocatable :: varfile
      !! output file for variables
    character(:), allocatable :: cachefile
      !! output file for cache
  contains
    procedure :: init              => small_signal_init
    procedure :: set_solver_params => small_signal_set_solver_params
    procedure :: init_output       => small_signal_init_output
    procedure :: init_cache        => small_signal_init_cache
    procedure :: run               => small_signal_run
    procedure :: get_scalar        => small_signal_get_scalar

    procedure, private :: output   => small_signal_output
  end type

contains

  subroutine small_signal_init(this, sys, use_ram, log)
    !! initialize small-signal simulation
    class(small_signal),    intent(out) :: this
    type(esystem), target,  intent(in)  :: sys
      !! equation system
    logical,      optional, intent(in)  :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical,      optional, intent(in)  :: log
      !! enable/disable logging (default false)

    this%sys => sys
    this%use_ram = .true.
    if (present(use_ram)) this%use_ram = use_ram
    this%log = .false.
    if (present(log)) this%log = log

    this%solver = SPSOLVER_PARDISO
  end subroutine

  subroutine small_signal_set_solver_params(this, solver, gopt)
    !! initialize small-signal simulation
    class(small_signal),           intent(inout) :: this
    integer,             optional, intent(in)    :: solver
      !! matrix solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for gmres in case an iterative solver is used

    if (present(solver)) this%solver = solver
    if (present(gopt)) this%gopt = gopt
  end subroutine

  subroutine small_signal_init_output(this, vars, varfile)
    !! add main variable to output
    class(small_signal), intent(inout) :: this
    type(string),        intent(in)    :: vars(:)
      !! names of the main variables
    character(*),        intent(in)    :: varfile
      !! file to which the variables are written

    integer            :: i, jimvar, k
    logical            :: found

    ! check that input is valid: small-signal parameters are only calculated for main variables
    m4_assert(size(vars) /= 0)
    do i = 1, size(vars)
      found = .false.
      do jimvar = 1, this%sys%g%imvar%n
        associate (n => this%sys%g%nodes%d(this%sys%g%imvar%d(jimvar))%p)
          do k = 1, size(n%v%v)
            if (n%v%v(k)%p%name == vars(i)%s) found = .true.
          end do
        end associate
      end do
      if (.not. found) call program_error("the esystem has no main variable named " // vars(i)%s)
    end do

    ! delete old output settings
    if (allocated(this%output_vars)) deallocate(this%output_vars)
    if (allocated(this%varfile)) deallocate(this%varfile)

    ! new output settings
    this%output_vars = vars
    this%varfile = varfile
  end subroutine

  subroutine small_signal_init_cache(this, cachefile)
    !! init output for solution vector this%sys%get_x()
    class(small_signal), intent(inout) :: this
    character(*),        intent(in)    :: cachefile
      !! file to which the variables are written

    integer :: unit, stat

    ! delete old output settings
    if (allocated(this%cachefile)) deallocate(this%cachefile)

    ! new output settings
    this%cachefile = cachefile
    ! if cachefile already exists, delete it
    open(newunit = unit, file = this%cachefile, iostat = stat, status = "replace")
    close(unit, status="delete")
  end subroutine

  subroutine small_signal_run(this, s, calc_dxds)
    !! perform small-signal analysis
    class(small_signal), intent(inout) :: this
    complex,             intent(in)    :: s(:)
      !! assume: x = x_0 + Re{x_1 * exp(s*t)}
    logical, optional,   intent(in)    :: calc_dxds
      !! calculate dxds in addition to x (default: false)

    complex, allocatable     :: rhs(:,:), x(:,:), tmp(:,:), dxds(:,:)
    integer                  :: nsrc, ns, i, j, k, stat
    logical                  :: calc_dxds_
    type(block_cmplx)        :: dfdx, dft, dfdx_prec, mat, mat_prec
    type(container)          :: c_vars
    type(single_matop_cmplx) :: mulvec, precon
    type(variable_ptr)       :: ptr

    ! optional arguments
    calc_dxds_ = .false.
    if (present(calc_dxds)) calc_dxds_ = calc_dxds

    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%s)) deallocate(this%s)

    ! input
    nsrc = this%sys%ninput
    ns = size(s)
    this%s = s

    ! evaluate equation system
    call this%sys%eval()

    ! get jacobian and jacobian wrt time derivative
    call matrix_convert(this%sys%df, dfdx)
    call matrix_convert(this%sys%dft, dft)
    if (this%solver == SOLVER_GMRES) call matrix_convert(this%sys%dfp, dfdx_prec)

    ! allocate memory
    if (this%use_ram) allocate (this%x(this%sys%n,nsrc,ns))
    if (this%use_ram .and. calc_dxds_) allocate (this%dxds(this%sys%n,nsrc,ns))

    ! allocate temporary memory
    allocate (rhs(dfdx%nrows,nsrc), source = (0.0, 0.0))
    allocate (  x(dfdx%nrows,nsrc), source = (0.0, 0.0))
    if (calc_dxds_) then
      allocate ( tmp(dfdx%nrows,nsrc), source = (0.0, 0.0))
      allocate (dxds(dfdx%nrows,nsrc), source = (0.0, 0.0))
    end if

    ! set right-hand sides
    k = 0
    do i = 1, this%sys%input_equs%n
      do j = this%sys%input_i0(i), this%sys%input_i1(i)
        k = k + 1
        rhs(j,k) = 1.0
      end do
    end do

    ! output the grids on which the output variables are defined
    if (allocated(this%varfile)) then
      call c_vars%open(this%varfile, flag = STORAGE_WRITE)
      do j = 1, size(this%output_vars)
        ptr = this%sys%search_var(this%output_vars(j)%s)
        ! save grid if it is not already saved
        call c_vars%save(ptr%p%g, stat = stat)
        if (stat /= 0 .and. stat /= ERR_ALREADY_EXISTS) call program_error("Error when saving grids")
      end do
      call c_vars%close()
    end if

    ! get small-signal quantities for all frequencies
    do i = 1, ns
      if (this%log) print *, "small-signal step " // int2str(i) // " of " // int2str(ns)

      if (this%solver == SOLVER_GMRES) then
        ! initial solution
        x = 0
        ! construct matrix (mat = s*M + dfdx(x0))
        call matrix_add(dfdx,      dft, mat,      fact2 = s(i))
        call matrix_add(dfdx_prec, dft, mat_prec, fact2 = s(i))
        ! init gmres
        call mat_prec%factorize(solver = this%gopt%solver)
        call mulvec%init(mat)
        call precon%init(mat_prec, inv = .true.)
        ! solve for dx by gmres
        do j = 1, size(x, dim = 2)
          call gmres(rhs(:,j), mulvec, x(:,j), opts = this%gopt, precon = precon)
        end do
      else
        ! construct matrix (mat = s*M + dfdx(x0))
        call matrix_add(dfdx, dft, mat, fact2 = s(i))
        ! factorize and solve
        call mat%factorize(solver = this%solver)
        call mat%solve_mat(rhs, x)
      end if

      ! extract small-signal data
      if (this%use_ram) this%x(:,:,i) = x(:,:)

      ! derivatives
      if (calc_dxds_) then
        ! tmp = dft * x
        call dft%mul_mat(x, tmp)

        ! dxds = mat \ (- dft * x)
        if (this%solver == SOLVER_GMRES) then
          dxds = 0
          do j = 1, size(x, dim = 2)
            call gmres(-tmp(:,j), mulvec, dxds(:,j), opts = this%gopt, precon = precon)
          end do
        else
          call mat%solve_mat(-tmp, dxds)
        end if

        ! extract derivatives
        if (this%use_ram) this%dxds(:,:,i) = dxds(:,:)
      end if

      ! release memory
      call mat%destruct()

      ! output
      call this%output(i, x, dxds)
    end do
  end subroutine

  function small_signal_get_scalar(this, var_name, input_name) result (d)
    !! return small-signal data for given scalar main variable and input
    class(small_signal), intent(inout) :: this
    character(*),        intent(in)    :: var_name
      !! name of the result main variable
    character(*),        intent(in)    :: input_name
      !! name of the input variable
    complex, allocatable               :: d(:)
      !! small-signal data d_var/d_input for each frequency

    integer            :: i_in, i_out, j, k, l, cntr, iimvar, itab, ipoint, ivar
    type(variable_ptr) :: ptr_input

    if (.not. allocated(this%x)) call program_error("no small-signal data stored on RAM!")

    i_in = 0
    i_out = 0

    ! get input variable position
    l = 0
    do j = 1, this%sys%input_equs%n
      associate (vsel => this%sys%g%nodes%d(this%sys%g%equs%d(this%sys%input_equs%d(j))%imain)%p%v)
        do k = this%sys%input_i0(j), this%sys%input_i1(j)
          l = l + 1
          ptr_input = vsel%v(mod(k-this%sys%input_i0(j), vsel%nvar) + 1)
          if (ptr_input%p%name == input_name) then
            if (i_in /= 0) call program_error("Input variable is not scalar")
            i_in = l
          end if
        end do
      end associate
    end do
    if (i_in == 0) call program_error("Input variable not found")
    ! get main variable position
    cntr = 0
    do iimvar = 1, this%sys%g%imvar%n
      associate (v => this%sys%g%nodes%d(this%sys%g%imvar%d(iimvar))%p%v)
        do itab = 1, v%ntab
          do ipoint = 1, v%tab(itab)%p%n
            do ivar = 1, v%nvar
              cntr = cntr + 1
              if (v%v(ivar)%p%name == var_name) then
                if (i_out /= 0) call program_error("Main variable is not scalar")
                i_out = cntr
              end if
            end do
          end do
        end do
      end associate
    end do
    if (i_out == 0) call program_error("Main variable not found")
    d = this%x(i_out, i_in, :)
  end function

  subroutine small_signal_output(this, i, x, dxds)
    !! write x and dxds to grid data sorted by input and output variable and save to this%varfile and this%cachefile
    class(small_signal), intent(inout) :: this
    integer,             intent(in)    :: i
      !! frequency index
    complex,             intent(in)    :: x(:,:)
      !! small-signal data
    complex, optional,   intent(in)    :: dxds(:,:)
      !! derivatives of x wrt s

    integer                             :: iout, j, k, l, cntr, iimvar, itab, ipoint, ivar
    integer, allocatable                :: idx(:)
    type(container)                     :: c_vars, c_cache
    class(grid_data_cmplx), allocatable :: gdata, gdata_dxds
    type(variable_ptr)                  :: ptr, ptr_input

    if (allocated(this%varfile)) then
      call c_vars%open(this%varfile, flag = STORAGE_WRITE)
      call c_vars%write("small-signal/s", [this%s(i)], unit = "Hz", dynamic = DYNAMIC_APP)
      ! iout: loop over output variables
      do iout = 1, size(this%output_vars)
        ptr = this%sys%search_var(this%output_vars(iout)%s)
        ! j, k, l: loop over input variables
        l = 0
        do j = 1, this%sys%input_equs%n
          associate (vsel => this%sys%g%nodes%d(this%sys%g%equs%d(this%sys%input_equs%d(j))%imain)%p%v)
            do k = this%sys%input_i0(j), this%sys%input_i1(j)
              l = l + 1
              if (this%sys%input_i1(j) - this%sys%input_i0(j) + 1 /= vsel%nvar) call program_error("automatic output only implemented for scalar input variables")
              ptr_input = vsel%v(mod(k-this%sys%input_i0(j), vsel%nvar) + 1)
              call allocate_grid_data0_cmplx(gdata, ptr%p%g%idx_dim)
              if (present(dxds)) call allocate_grid_data0_cmplx(gdata_dxds, ptr%p%g%idx_dim)
              call gdata%init(ptr%p%g, ptr%p%idx_type, ptr%p%idx_dir, (0.0, 0.0) * ieee_value(1.0, ieee_quiet_nan))
              if (present(dxds)) call gdata_dxds%init(ptr%p%g, ptr%p%idx_type, ptr%p%idx_dir, (0.0, 0.0) * ieee_value(1.0, ieee_quiet_nan))
              ! iimvar, itab, ipoint, ivar: loop over whole solution vector
              cntr = 0
              do iimvar = 1, this%sys%g%imvar%n
                associate (v => this%sys%g%nodes%d(this%sys%g%imvar%d(iimvar))%p%v)
                  do itab = 1, v%ntab
                    do ipoint = 1, v%tab(itab)%p%n
                      idx = v%tab(itab)%p%get_idx(ipoint)
                      do ivar = 1, v%nvar
                        cntr = cntr + 1
                        if (v%v(ivar)%p%name == this%output_vars(iout)%s) call gdata%set(idx, x(cntr, l))
                        if (v%v(ivar)%p%name == this%output_vars(iout)%s .and. present(dxds)) call gdata_dxds%set(idx, dxds(cntr, l))
                      end do
                    end do
                  end do
                end associate
              end do
              call c_vars%save("small-signal/x/d"//this%output_vars(iout)%s//"_d"//ptr_input%p%name, gdata, unit = ptr%p%unit//"/("//ptr_input%p%unit//")", dynamic = .true.)
              if (present(dxds)) call c_vars%save("small-signal/dxds/d"//this%output_vars(iout)%s//"_d"//ptr_input%p%name, gdata_dxds, unit = ptr%p%unit//"/("//ptr_input%p%unit//")*s", dynamic = .true.)
              deallocate (gdata)
              if (present(dxds)) deallocate (gdata_dxds)
            end do
          end associate
        end do
      end do
      call c_vars%close()
    end if
    if (allocated(this%cachefile)) then
      call c_cache%open(this%cachefile, flag = STORAGE_WRITE)
      call c_cache%write("small-signal/s", [this%s(i)], unit = "Hz", dynamic = DYNAMIC_APP)
      call c_cache%write("small-signal/x", x, dynamic = DYNAMIC_EXT)
      if (present(dxds)) call c_cache%write("small-signal/dxds", dxds, dynamic = DYNAMIC_EXT)
      call c_cache%close()
    end if
  end subroutine

end module
