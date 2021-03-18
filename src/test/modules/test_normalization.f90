module test_normalization_m
  use test_case_m
  use math_m
  use normalization_m
#ifdef __INTEL_COMPILER
  use ifport
#endif
  implicit none

  private
  public :: test_normalization

contains

  subroutine test_normalization()
    type(test_case) :: tc
    real, parameter :: T_K = 300.0

    call tc%init("normalization")

    call init_normconst(T_K)

    ! check if norming and denorming is exchanged
    block
      real :: ener_eV, ener, ener_exp

      ener_eV  = 5.0                                ! 5eV
      ener_exp = 5.0 / (8.617333262145e-5 * 300.0)  ! = ener_eV / (kB_eV_per_K * T_K) = approx 2e3

      ener = norm(ener_eV, "eV")

      call tc%assert_eq(ener_exp, ener, 2e3*1e-13, "norming 5eV")
    end block

    ! check norm*denorm = identity
    block
      character(:), allocatable         :: unit
      character(*), parameter           :: unit_cm = "cm", unit_eV = "eV", unit_A_cm = "A/cm"
      integer                           :: i, j, n, i_unit
      real                              :: val_scal, val_scal_exp, nval_scal
      real, dimension(:),   allocatable :: val_arr, val_arr_exp, nval_arr
      real, dimension(:,:), allocatable :: val_mat, val_mat_exp, nval_mat

      n = 4
      allocate(val_arr(n),   val_arr_exp(n),   nval_arr(n))
      allocate(val_mat(n,n), val_mat_exp(n,n), nval_mat(n,n))

      allocate (character(0) :: unit)      ! remove gfortran warning

      do i_unit = 1, 3
        select case (i_unit)
          case (1)
            unit = unit_A_cm
          case (2)
            unit = unit_cm
          case (3)
            unit = unit_eV
        end select

        ! scalar
        val_scal_exp = rand()
        nval_scal    = norm(val_scal_exp, unit)
        val_scal     = denorm(nval_scal, unit)
        call tc%assert_eq(val_scal_exp, val_scal, 1e-13, "denorm(norm(scalar))")

        ! array
        do i = 1, n
          val_arr_exp(i) = rand()
        end do
        nval_arr = norm(val_arr_exp, unit)
        val_arr  = denorm(nval_arr, unit)
        call tc%assert_eq(val_arr_exp, val_arr, 1e-13, "denorm(norm(array))")

        ! matrix
        do i = 1, n
          do j = 1, n
            val_mat_exp(i,j) = rand()
          end do
        end do
        nval_mat = norm(val_mat_exp, unit)
        val_mat  = denorm(nval_mat, unit)
        call tc%assert_eq(val_mat_exp, val_mat, 1e-13, "denorm(norm(matrix))")
      end do
    end block

    ! check deg
    block
      real :: deg

      deg = norm(45.0, "deg")
      call tc%assert_eq(PI/4.0, deg, 1e-14, "norming degrees")
    end block

    call destruct_normconst

    call tc%finish
  end subroutine

end module
