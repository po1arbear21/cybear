m4_include(../util/macro.f90.inc)

module eigenvalues_m

  use error_m,   only: assert_failed
  use esystem_m, only: esystem
  use lapack95,  only: ggev
  use math_m,    only: PI
  use matrix_m,  only: sparse_real, dense_real, spbuild_real, block_real, matrix_convert

  implicit none

  private
  public :: eigenvalues

  type eigenvalues
    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system

    complex, allocatable :: s(:)
      !! complex 'frequency': x = x_0 + Re{x_1 * exp(s*t)}
      !! s = sigma + j omega
    complex, allocatable :: x(:,:)
      !! small-signal data (sys%n, sys%n)
  contains
    procedure :: run       => eigenvalues_run
    procedure :: schur     => eigenvalues_schur
    procedure :: partition => eigenvalues_partition
  end type

contains

 subroutine eigenvalues_run(this, sys)
    !! run eigenvalue analysis
    class(eigenvalues),    intent(out)   :: this
    type(esystem), target, intent(inout) :: sys
      !! equation system to analyse

    type(sparse_real)    :: df, dft
    type(dense_real)     :: dfs, dfts
    integer, allocatable :: ipart(:)
    real,    allocatable :: ar(:), ai(:), bb(:)
    integer              :: i, j, n, ns

    this%sys => sys

    call sys%eval()
    call sys%get_df(df)
    call sys%get_dft(dft)

    ! eliminate zero rows from dft
    call this%schur(df, dft, dfs, dfts, ipart)

    ! system size
    n = dfs%nrows

    ! solve generalized eigenvalue problem
    allocate (ar(n), source = 0.0)
    allocate (ai(n), source = 0.0)
    allocate (bb(n), source = 0.0)
    dfts%d = - dfts%d
    call ggev(dfs%d, dfts%d, ar, ai, bb)

    ! get number of s values
    ns = count(bb /= 0.0)

    ! construct complex eigenvalues
    allocate (this%s(ns), source = (0.0, 0.0))
    j = 0
    do i = 1, n
      if (bb(i) /= 0.0) then
        j = j + 1
        this%s(j) = cmplx(ar(i) / bb(i), ai(i) / bb(i))
      end if
    end do
  end subroutine

  subroutine eigenvalues_schur(this, df, dft, dfs, dfts, ipart)
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

    type(sparse_real)    :: df21, df22, dft21, dft22
    type(dense_real)     :: df1112
    integer              :: i, j, n, m
    integer, allocatable :: rtab(:)
    logical              :: zero
    type(sparse_real)    :: df11, df12, dft11, dft12
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
    call df11%factorize()
    call df11%solve_mat(rhs%d(:,1:j), df1112%d(:,1:j))
    call df11%destruct()

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

end module
