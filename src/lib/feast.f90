module feast_m
  use iso_fortran_env, only: real64

  implicit none

  private

  interface
    subroutine feastinit(fpm)
      integer, dimension(*) :: fpm
    end subroutine feastinit

    subroutine zfeast_grci(ijob, N, Ze, work, workc, zAq, zBq, fpm, epsout, loop, Emid, r, M0, lambda, Q, mode, res, info)
      import real64
      integer                          :: ijob
      integer                          :: N
      complex(real64)                  :: Ze
      complex(real64), dimension(N,*)  :: work
      complex(real64), dimension(N,*)  :: workc
      complex(real64), dimension(M0,*) :: zAq
      complex(real64), dimension(M0,*) :: zBq
      integer,         dimension(*)    :: fpm
      real(real64)                     :: epsout
      integer                          :: loop
      complex(real64)                  :: Emid
      real(real64)                     :: r
      integer                          :: M0
      complex(real64), dimension(*)    :: lambda
      complex(real64), dimension(N,*)  :: Q
      integer                          :: mode
      real(real64),    dimension(*)    :: res
      integer                          :: info
    end subroutine
  end interface

contains

  subroutine feast_solve()
    integer :: fpm(64)
    integer :: info

    call feastinit(fpm)

    info = 0
    call check_info(info)
  end subroutine

  subroutine check_info(info)
    integer, intent(in) :: info

    character(*), parameter :: ERROR_20X(200:202) = [ &
      "Problem with Emin, Emax or Emid, r",           &
      "Problem with size of subspace M0  ",           &
      "Problem with size of the system N "            &
    ]
    character(*), parameter :: ERROR_10X = "Problem with ith value of the input FEAST parameter (i.e fpm(i))"
    character(*), parameter :: WARNING_X(1:7) = [                                          &
      "No Eigenvalue found in the search interval                                       ", &
      "No Convergence (#iteration loops>fpm(4))                                         ", &
      "Size of the subspace M0 is too small (M0<=M)                                     ", &
      "Only the subspace has been returned using fpm(14)=1                              ", &
      "Only stochastic estimation of #eigenvalues returned fpm(14)=2                    ", &
      "FEAST converges but subspace is not bi-orthonormal                               ", &
      "The search for extreme eigenvalues has failed, search contour must be set by user"  &
    ]
    character(*), parameter :: ERROR_NEG_X(-3:-1) = [                         &
      "Internal error of the reduced eigenvalue solver                     ", &
      "Internal error of the inner system solver in FEAST Driver interfaces", &
      "Internal error conversion single/double                             "  &
    ]
    character(*), parameter :: ERROR_NEG_10X = "Problem with the ith argument of the FEAST interface"

  end subroutine

end module
