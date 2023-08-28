m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_spike},0,-1))

module spike_m

  use error_m, only: program_error

  implicit none

  private
  public spike_params, spike_factorize, spike_solve

  type spike_params
    integer, allocatable :: spm(:)
  contains
    procedure :: init => spike_params_init
  end type

  interface spike_factorize
    module procedure :: spike_factorize_D
    module procedure :: spike_factorize_Z
  end interface

  interface spike_solve
    module procedure :: spike_solve_D
    module procedure :: spike_solve_Z
  end interface

  ! interface to external library
  interface
    subroutine spikeinit(spm, n, klu)
      !! initialize spike parameters to sane defaults
      !! uses maximal number of threads
      integer, intent(out) :: spm(64)
        !! output spike parameters
      integer, intent(in)  :: n
        !! matrix size
      integer, intent(in)  :: klu
        !! max of upper and lower bandwidth
    end subroutine

    subroutine spikeinit_nthread(spm, n, klu, nthread)
      !! initialize spike parameters to sane defaults
      integer, intent(out) :: spm(64)
        !! output spike parameters
      integer, intent(in)  :: n
        !! matrix size
      integer, intent(in)  :: klu
        !! max of upper and lower bandwidth
      integer, intent(in)  :: nthread
        !! number of threads to use
    end subroutine

    subroutine spikeinit_default(spm)
      !! initialize spike parameters to sane defaults
      integer, intent(out) :: spm(64)
        !! output spike parameters
    end subroutine

    m4_define({m4_list},{
      m4_X(D,real)
      m4_X(Z,complex)
    })

    m4_define({m4_X},{

    subroutine $1spike_tune(spm)
      !! auto-tune partition size ratios, sets spm(4), spm(5) and spm(7)
      !! slow, but results only depend on hardware
      !! run only once and reuse values
      integer, intent(inout) :: spm(64)
    end subroutine

    subroutine $1spike_gbsv(spm, n, kl, ku, nrhs, A, lda, B, ldb, info)
      !! all-in-one factorization, solving and iterative refinement (A * X = B)
      !! matrix A returned to initial state
      integer, intent(in)    :: spm(64)
        !! spike parameters
      integer, intent(in)    :: n
        !! system size = number of rows and cols of A
      integer, intent(in)    :: kl
        !! lower bandwidth of A
      integer, intent(in)    :: ku
        !! upper bandwidth of A
      integer, intent(in)    :: nrhs
        !! number of right-hand sides
      $2,      intent(in)    :: A(lda,*)
        !! matrix A
      integer, intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,      intent(inout) :: B(ldb,*)
        !! right hand sides, overwritten by solution
      integer, intent(in)    :: ldb
        !! leading dimension of B
      integer, intent(out)   :: info
        !! output information
        !! 0: fact. required no boosting
        !! 1: fact. required boosting, iterative refinement maybe necessary
        !! 2: illegal value in A matrix
    end subroutine

    subroutine $1spike_gbtrf(spm, n, kl, ku, A, lda, work, info)
      !! factorize matrix A
      integer, intent(in)    :: spm(64)
        !! spike parameters
      integer, intent(in)    :: n
        !! system size = number of rows and cols of A
      integer, intent(in)    :: kl
        !! lower bandwidth of A
      integer, intent(in)    :: ku
        !! upper bandwidth of A
      $2,      intent(inout) :: A(lda,*)
        !! matrix A, overwritten by factorization
      integer, intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,      intent(inout) :: work(*)
        !! work array of size klu**2 * spm(10)
      integer, intent(out)   :: info
        !! output information
        !! 0: fact. required no boosting
        !! 1: fact. required boosting, iterative refinement maybe necessary
        !! 2: illegal value in A matrix
    end subroutine

    subroutine $1spike_gbtrs(spm, trans, n, kl, ku, nrhs, A, lda, work, B, ldb)
      !! solution step
      integer,   intent(in)    :: spm(64)
        !! spike parameters
      character, intent(in)    :: trans
        !! 'N': use matrix as is
        !! 'T': use transposed matrix
        !! 'C': use transposed+complex conjugated matrix
      integer,   intent(in)    :: n
        !! system size = number of rows and cols of A
      integer,   intent(in)    :: kl
        !! lower bandwidth of A
      integer,   intent(in)    :: ku
        !! upper bandwidth of A
      integer,   intent(in)    :: nrhs
        !! number of right-hand sides
      $2,        intent(in)    :: A(lda,*)
        !! factorized matrix A
      integer,   intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,        intent(inout) :: work(*)
        !! work array of size klu**2 * spm(10)
      $2,        intent(inout) :: B(ldb,*)
        !! right hand sides, overwritten by solution
      integer,   intent(in)    :: ldb
        !! leading dimension of B
    end subroutine

    subroutine $1spike_gbtrsi(spm, trans, n, kl, ku, nrhs, C, ldc, A, lda, work, B, ldb)
      !! solution step + iterative refinement
      integer,   intent(in)    :: spm(64)
        !! spike parameters
      character, intent(in)    :: trans
        !! 'N': use matrix as is
        !! 'T': use transposed matrix
        !! 'C': use transposed+complex conjugated matrix
      integer,   intent(in)    :: n
        !! system size = number of rows and cols of A
      integer,   intent(in)    :: kl
        !! lower bandwidth of A
      integer,   intent(in)    :: ku
        !! upper bandwidth of A
      integer,   intent(in)    :: nrhs
        !! number of right-hand sides
      $2,        intent(in)    :: C(ldc,*)
        !! original matrix A
      integer,   intent(in)    :: ldc
        !! leading dimension of C (ldc  = kl + ku + 1)
      $2,        intent(in)    :: A(lda,*)
        !! factorized matrix A
      integer,   intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,        intent(inout) :: work(*)
        !! work array of size klu**2 * spm(10)
      $2,        intent(inout) :: B(ldb,*)
        !! right hand sides, overwritten by solution
      integer,   intent(in)    :: ldb
        !! leading dimension of B
    end subroutine

    subroutine $1spike_gbtrfp(spm, n, kl, ku, A, lda, work, ipiv, info)
      !! factorization + pivoting
      integer, intent(in)    :: spm(64)
        !! spike parameters
      integer, intent(in)    :: n
        !! system size = number of rows and cols of A
      integer, intent(in)    :: kl
        !! lower bandwidth of A
      integer, intent(in)    :: ku
        !! upper bandwidth of A
      $2,      intent(inout) :: A(lda,*)
        !! matrix A, overwritten by factorization
      integer, intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,      intent(inout) :: work(*)
        !! work array of size klu**2 * spm(10)
      integer, intent(out)   :: ipiv(*)
        !! pivoting array of size n
      integer, intent(out)   :: info
        !! output information
        !! 0: fact. required no boosting
        !! 1: fact. required boosting, iterative refinement maybe necessary
        !! 2: illegal value in A matrix
    end subroutine

    subroutine $1spike_gbtrsp(spm, n, kl, ku, nrhs, A, lda, work, ipiv, B, ldb)
      !! solution step including pivoting
      integer, intent(in)    :: spm(64)
        !! spike parameters
      integer, intent(in)    :: n
        !! system size = number of rows and cols of A
      integer, intent(in)    :: kl
        !! lower bandwidth of A
      integer, intent(in)    :: ku
        !! upper bandwidth of A
      integer, intent(in)    :: nrhs
        !! number of right-hand sides
      $2,      intent(in)    :: A(lda,*)
        !! factorized matrix A
      integer, intent(in)    :: lda
        !! leading dimension of A (lda  = kl + ku + 1)
      $2,      intent(inout) :: work(*)
        !! work array of size klu**2 * spm(10)
      integer, intent(in)    :: ipiv(*)
        !! pivoting array of size n
      $2,      intent(inout) :: B(ldb,*)
        !! right hand sides, overwritten by solution
      integer, intent(in)    :: ldb
        !! leading dimension of B
    end subroutine

    })
    m4_list

  end interface

contains

  subroutine spike_params_init(this, cmplx, n, klu, nthread)
    !! initialize spike parameters
    class(spike_params), intent(out) :: this
    logical,             intent(in)  :: cmplx
      !! complex or real?
    integer,             intent(in)  :: n
      !! system size (number of rows and cols)
    integer,             intent(in)  :: klu
      !! max of upper and lower bandwidth
    integer, optional,   intent(in)  :: nthread
      !! optional: number of threads to use, default: maximum number is used

    integer, save :: rtune(3) = [-1, -1, -1]
    integer, save :: ctune(3) = [-1, -1, -1]

    allocate (this%spm(64), source = 0)

    ! use library to set sane default values
    if (present(nthread)) then
      call spikeinit_nthread(this%spm, n, klu, nthread)
    else
      call spikeinit(this%spm, n, klu)
    end if

    ! tune
    if (cmplx) then
      if (any(ctune < 0)) then
        call zspike_tune(this%spm)
        ctune(1) = this%spm(4)
        ctune(2) = this%spm(5)
        ctune(3) = this%spm(7)
      else
        this%spm(4) = ctune(1)
        this%spm(5) = ctune(2)
        this%spm(7) = ctune(3)
      end if
    else
      if (any(rtune < 0)) then
        call dspike_tune(this%spm)
        rtune(1) = this%spm(4)
        rtune(2) = this%spm(5)
        rtune(3) = this%spm(7)
      else
        this%spm(4) = rtune(1)
        this%spm(5) = rtune(2)
        this%spm(7) = rtune(3)
      end if
    end if
  end subroutine

  m4_define({m4_X},{

  subroutine spike_factorize_$1(sp, n, kl, ku, A, ipiv, info)
    type(spike_params), intent(in)    :: sp
    integer,            intent(in)    :: n
    integer,            intent(in)    :: kl
    integer,            intent(in)    :: ku
    $2,                 intent(inout) :: A(:,:)
    integer,            intent(out)   :: ipiv(:)
    integer, optional,  intent(out)   :: info

    integer         :: info_
    $2, allocatable :: work(:)

    allocate (work(max(kl,ku)**2 * sp%spm(10)))
    call $1spike_gbtrfp(sp%spm, n, kl, ku, A, max(1,size(A,1)), work, ipiv, info_)

    if (present(info)) then
      info = info_
    else
      if (info == 1) then
        print "(A)", "SPIKE Warning: Iterative refinement maybe necessary"
      end if
      if (info == 2) then
        call program_error("SPIKE: Illegal value in matrix A")
      end if
    end if
  end subroutine

  subroutine spike_solve_$1(sp, n, kl, ku, A, ipiv, b)
    type(spike_params), intent(in)    :: sp
    integer,            intent(in)    :: n
    integer,            intent(in)    :: kl
    integer,            intent(in)    :: ku
    $2,                 intent(in)    :: A(:,:)
    integer,            intent(in)    :: ipiv(:)
    $2,                 intent(inout) :: b(:,:)

    $2, allocatable :: work(:)

    allocate (work(max(kl,ku)**2 * sp%spm(10)))
    call $1spike_gbtrsp(sp%spm, n, kl, ku, size(b,2), A, size(A,1), work, ipiv, b, size(b,1))
  end subroutine

  })
  m4_list

end module

m4_divert(0)
