m4_include(../util/macro.f90.inc)

module eigenvalues_m

  use error_m,   only: assert_failed
  use esystem_m, only: esystem
  use lapack95,  only: ggev
  use math_m,    only: PI
  use matrix_m,  only: sparse_real, dense_real, spbuild_real, block_real, matrix_convert, sparse_cmplx, matrix_add
  m4_divert(m4_ifdef({m4_feast},0,-1))
  use feast_m,  only: feast_option, feast
  m4_divert(0)
  use gmres_m,  only: gmres_options, gmres
  use matop_m,  only: single_matop_cmplx


  implicit none

  private
  public :: eigenvalues

  type eigenvalues
    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system

    integer              :: ns
      !! number of eigenvalues
    complex, allocatable :: s(:)
      !! complex 'frequency': x = x_0 + Re{x_1 * exp( s*t)}
      !! s = sigma + j omega
    complex, allocatable :: x(:,:)
      !! small-signal data (sys%n, ns)
  contains
    m4_divert(m4_ifdef({m4_feast},0,-1))
    procedure :: run_feast   => eigenvalues_run_feast
    m4_divert(0)
    procedure :: run_dense    => eigenvalues_run_dense
    procedure :: schur        => eigenvalues_schur
    procedure :: partition    => eigenvalues_partition
    procedure :: select_real  => eigenvalues_select_real
    procedure :: select_imag  => eigenvalues_select_imag
    procedure :: select_abs   => eigenvalues_select_abs
    procedure :: select_phase => eigenvalues_select_phase
  end type

contains

m4_divert(m4_ifdef({m4_feast},0,-1))
subroutine eigenvalues_run_feast(this, sys, M0, abcd)
    class(eigenvalues),    intent(out)   :: this
    type(esystem), target, intent(inout) :: sys
      !! equation system to analyse
    integer,               intent(in)    :: M0
      !! overestimate of number of eigenvalues
    complex, optional,     intent(in)    :: abcd(4)

    logical                  :: precon
    complex                  :: abcd_(4)
    type(feast_option)       :: feast_opt
    type(gmres_options)      :: gmres_opt
    type(sparse_cmplx)       :: df, dfp, dft, A, Ap, B, Bp, ZeB_A, ZeB_Ap
    type(single_matop_cmplx) :: ZeB_A_mul, ZeB_Ap_solve

    abcd_(1) = cmplx(1.0)
    abcd_(2) = cmplx(-1.0)
    abcd_(3) = cmplx(1.0)
    abcd_(4) = cmplx(1.0)

    if (present(abcd)) abcd_ = abcd

    precon = allocated(sys%dfp)

    call sys%get_df(df)
    call sys%get_dft(dft)

    call matrix_add(df,  dft, A,  abcd_(1), -abcd_(2))
    call matrix_add(df,  dft, B,  abcd_(3), -abcd_(4))

    call feast_opt%init(df%nrows, M0)

    if (precon) then
      call sys%get_dfp(dfp)
      call matrix_add(dfp, dft, Ap, abcd_(1), -abcd_(2))
      call matrix_add(dfp, dft, Bp, abcd_(3), -abcd_(4))

      gmres_opt%print_msg = .true.
      gmres_opt%rtol = 1e-6

      call ZeB_A_mul%init(ZeB_A)
      call ZeB_Ap_solve%init(ZeB_Ap, inv=.true.)
    end if

    call feast(feast_opt, fact_Az, solve_Az, mulvec_A, this%s, evecR = this%x, mulvec_B = mulvec_B)
    this%s  = (abcd_(4) * this%s - abcd_(2)) / (abcd_(1) - abcd_(3) * this%s)
    this%ns = size(this%s)

    if (ZeB_A%factorized) call ZeB_A%destruct()
    if (ZeB_Ap%factorized) call ZeB_Ap%destruct()

  contains

    subroutine fact_Az(Ze, ctrans)
      complex,           intent(in) :: Ze
      logical, optional, intent(in) :: ctrans
        !! factorize (Ze*B-A)^H. default: .false.

      if (ZeB_A%factorized) call ZeB_A%destruct()
      if (ZeB_Ap%factorized) call ZeB_Ap%destruct()

      ! mat = Ze * B - A
      call matrix_add(B,  A,  ZeB_A,  Ze, cmplx(-1.0))
      if (precon) then
        call matrix_add(Bp, Ap, ZeB_Ap, Ze, cmplx(-1.0))
        call ZeB_Ap%factorize()
      else
        call ZeB_A%factorize()
      end if
    end subroutine

    subroutine solve_Az(xb, ctrans)
      !! Solve (Ze*B-A)*x=b
      complex,           intent(inout) :: xb(:,:)
        !! input:  b
        !! output: x
      logical, optional, intent(in)    :: ctrans
        !! solve (Ze*B-A)^H*x=b. default: .false.

      integer              :: i
      complex, allocatable :: b(:,:)

      b = xb
      if (precon) then
        do i = 1, size(xb, dim = 2)
          call gmres(b(:,i), ZeB_A_mul, xb(:,i), precon = ZeB_Ap_solve)
        end do
      else
        call ZeB_A%solve_mat(b, xb)
      end if
    end subroutine

    subroutine mulvec_A(x, y, ctrans)
      !! Matrix-vector multiplication: y <- A*x
      complex,           intent(in)  :: x(:,:)
        !! input vector x
      complex,           intent(out) :: y(:,:)
        !! output vector y
      logical, optional, intent(in)  :: ctrans
        !! y <- A^H*x. default: .false.

      call A%mul_mat(x, y)
    end subroutine

    subroutine mulvec_B(x, y, ctrans)
      !! Matrix-vector multiplication: y <- B*x
      complex,           intent(in)  :: x(:,:)
        !! input vector x
      complex,           intent(out) :: y(:,:)
        !! output vector y
      logical, optional, intent(in)  :: ctrans
        !! y <- B^H*x. default: .false.

      call B%mul_mat(x, y)
    end subroutine

  end subroutine
  m4_divert(0)

  subroutine eigenvalues_run_dense(this, sys, schur)
    !! run eigenvalue analysis
    class(eigenvalues),    intent(out)   :: this
    type(esystem), target, intent(inout) :: sys
      !! equation system to analyse
    logical, optional,     intent(in)    :: schur

    type(sparse_real)    :: df, dft
    type(dense_real)     :: dfs, dfts, df1112
    integer, allocatable :: ipart(:)
    real,    allocatable :: ar(:), ai(:), bb(:), x1(:,:), x2(:,:)
    integer              :: i, j, k, n

    logical :: schur_

    this%sys => sys

    call sys%eval()
    call sys%get_df(df)
    call sys%get_dft(dft)

    schur_ = .false.
    if (present(schur)) schur_ = schur
    if (schur_) then
      ! eliminate zero rows from dft
      call this%schur(df, dft, dfs, dfts, ipart, df1112)
    else
      call matrix_convert(df, dfs)
      call matrix_convert(dft, dfts)
    end if

    ! system size
    n = dfs%nrows

    ! solve generalized eigenvalue problem
    allocate (ar(n), source = 0.0)
    allocate (ai(n), source = 0.0)
    allocate (bb(n), source = 0.0)
    allocate (x2(n,n), source = 0.0)
    dfts%d = - dfts%d
    call ggev(dfs%d, dfts%d, ar, ai, bb, x2)

    ! get number of s values
    this%ns = count(bb /= 0.0)

    ! construct complex eigenvalues
    allocate (this%s(this%ns), source = (0.0, 0.0))
    j = 0
    do i = 1, n
      if (bb(i) == 0.0) cycle
      j = j + 1
      this%s(j) = cmplx(ar(i) / bb(i), ai(i) / bb(i))
    end do

    ! eigenvectors
    allocate (this%x(sys%n,this%ns))

    if (schur_) then
      allocate (x1(sys%n-n,n))
      call df1112%mul_mat(-x2, x1)
      do i = 1, sys%n
        if (ipart(i) < 0) then
          k = 0
          do j = 1, n
            if (bb(j) == 0) cycle
            k = k + 1
            this%x(i,k) = x1(-ipart(i),j)
          end do
        else
          k = 0
          do j = 1, n
            if (bb(j) == 0) cycle
            k = k + 1
            this%x(i,k) = x2(ipart(i),j)
          end do
        end if
      end do
    else
      k = 0
      do j = 1, n
        if (bb(j) == 0) cycle
        k = k + 1
        this%x(:,k) = x2(:,j)
      end do
    end if
  end subroutine

  subroutine eigenvalues_schur(this, df, dft, dfs, dfts, ipart, df1112)
    !! get Schur-Complement of matrices by removing zero rows of dft
    class(eigenvalues),   intent(in)  :: this
    type(sparse_real),    intent(in)  :: df
      !!! jacobian
    type(sparse_real),    intent(in)  :: dft
      !!! jacobian wrt time derivatives
    type(dense_real),     intent(out) :: dfs
      !!! output schur complement of df  (remove zero rows of dft)
    type(dense_real),     intent(out) :: dfts
      !!! output schur complement of dft (remove zero rows of dft)
    integer, allocatable, intent(out) :: ipart(:)
      !!! output partition table
    type(dense_real),     intent(out) :: df1112
      !! df11 \ df12

    type(sparse_real)    :: df11, df12, df21, df22, dft11, dft12, dft21, dft22
    integer              :: i, j, n, m
    integer, allocatable :: rtab(:)
    logical              :: zero
    type(dense_real)     :: rhs

    ! init partition table (ipart)
    allocate (ipart(df%nrows))
    n = 0 ! number of rows to transform away
    m = 0 ! number of rows to keep
    do i = 1, dft%nrows
      ! check if row of dft is empty
      zero = .true.
      do j = dft%ia(i), dft%ia(i+1)-1
        if (dft%a(j) /= 0.0) then
          zero = .false.
          exit
        end if
      end do

      ! save entry in partition table
      if (zero) then
        n = n + 1
        ipart(i) = -n
      else
        m = m + 1
        ipart(i) = +m
      end if
    end do

     ! partition df and dft (dft11, dft12 should be empty)
    call this%partition(df , ipart, df11 , df12 , df21 , df22 )
    call this%partition(dft, ipart, dft11, dft12, dft21, dft22)

    ! get rhs = df12
    call rhs%init(df12%nrows, ncols = df12%ncols)
    call matrix_convert(df12, rhs)

    ! compress rhs (remove zero columns)
    allocate (rtab(m), source = -1)
    j = 0
    do i = 1, m
      if (any(rhs%d(:,i) /= 0)) then
        j = j + 1
        if (i > j) rhs%d(:,j) = rhs%d(:,i)
        rtab(i) = j
      end if
    end do

    ! get df11 \ df12 (dense matrix)
    call df1112%init(n, ncols = m)
    if (j > 0) then
      call df11%factorize()
      call df11%solve_mat(rhs%d(:,1:j), df1112%d(:,1:j))
      call df11%destruct()
    end if

    ! decompress df1112
    do i = m, 1, -1
      if (rtab(i) == -1) then
        df1112%d(:,i) = 0
      else
        df1112%d(:,i) = df1112%d(:,rtab(i))
      end if
    end do

    ! get dfs = df22 - df21 * (df11 \ df12) = df22 - df21 * df1112
    call dfs%init(df22%nrows, ncols = df22%ncols)
    call matrix_convert(df22, dfs)
    call df21%mul_mat(-df1112%d, dfs%d, fact_y = 1.0)

    ! get dfts_dense = dft22 - dft21 * (df11 \ df12) = dft22 - dft21 * df1112
    call dfts%init(dft22%nrows, ncols = dft22%ncols)
    call matrix_convert(dft22, dfts)
    call dft21%mul_mat(-df1112%d, dfts%d, fact_y = 1.0)
  end subroutine

  subroutine eigenvalues_partition(this, s, ipart, s11, s12, s21, s22)
    class(eigenvalues), intent(in)  :: this
    type(sparse_real),  intent(in)  :: s
      !! sparse matrix
    integer,            intent(in)  :: ipart(:)
      !! partition indices (s%nrows; <0: -ith row of S11, S12; >0: ith row of S21, S22)
    type(sparse_real),  intent(out) :: s11
      !! output matrix for ipart(rows)<0, ipart(cols)<0
    type(sparse_real),  intent(out) :: s12
      !! output matrix for ipart(rows)<0, ipart(cols)>0
    type(sparse_real),  intent(out) :: s21
      !! output matrix for ipart(rows)>0, ipart(cols)<0
    type(sparse_real),  intent(out) :: s22
      !! output matrix for ipart(rows)>0, ipart(cols)>0

    integer            :: n, m, i, j, si, sj, row1, row2
    type(spbuild_real) :: sb11, sb12, sb21, sb22

    m4_ignore(this)

    n = count(ipart < 0)
    m = s%nrows - n

    ! init matrices
    call S11%init(n, ncols = n)
    call S12%init(n, ncols = m)
    call S21%init(m, ncols = n)
    call S22%init(m, ncols = m)

    ! partition matrix
    row1 = 0
    row2 = 0
    call sb11%init(S11)
    call sb12%init(S12)
    call sb21%init(S21)
    call sb22%init(S22)
    do i = 1, s%nrows
      si = ipart(i)
      if (si < 0) then
        row1 = row1 + 1
      else
        row2 = row2 + 1
      end if

      do j = s%ia(i), s%ia(i+1)-1
        sj = ipart(s%ja(j))

        if ((si < 0) .and. (sj < 0)) then
          call sb11%set(row1, -sj, s%a(j), search = .false.)
        else if (si < 0) then
          call sb12%set(row1,  sj, s%a(j), search = .false.)
        else if (sj < 0) then
          call sb21%set(row2, -sj, s%a(j), search = .false.)
        else
          call sb22%set(row2,  sj, s%a(j), search = .false.)
        end if
      end do
    end do
    call sb11%save()
    call sb12%save()
    call sb21%save()
    call sb22%save()
  end subroutine

  subroutine eigenvalues_select_real(this, is)
    !! save real part of small-signal data (overwrite equation system variables)
    class(eigenvalues), intent(in) :: this
    integer,            intent(in) :: is
      !! frequency index

    call this%sys%set_x(real(this%x(:,is)))
  end subroutine

  subroutine eigenvalues_select_imag(this, is)
    !! save imaginary part of small-signal data (overwrite equation system variables)
    class(eigenvalues),  intent(in) :: this
    integer,             intent(in) :: is
      !! frequency index

    call this%sys%set_x(aimag(this%x(:,is)))
  end subroutine

  subroutine eigenvalues_select_abs(this, is)
    !! save absolute value of small-signal data (overwrite equation system variables)
    class(eigenvalues),  intent(in) :: this
    integer,             intent(in) :: is
      !! frequency index

    call this%sys%set_x(abs(this%x(:,is)))
  end subroutine

  subroutine eigenvalues_select_phase(this, is)
    !! save phase of small-signal data (overwrite equation system variables)
    class(eigenvalues),  intent(in) :: this
    integer,             intent(in) :: is
      !! frequency index

    call this%sys%set_x(atan2(aimag(this%x(:,is)), real(this%x(:,is))))
  end subroutine

end module
