submodule(test_matrix_m) test_band_m
  use matrix_m
  implicit none

contains

  module subroutine test_band
    real, dimension(5,5) :: d_band
    type(band_real)      :: tri, band, penta, diag
    type(test_case)      :: tc

    print "(A)", "test_band"
    call tc%init("band")

    ! init examples of band matrices + factorize
    block
      integer           :: i
      real, allocatable :: d0(:,:)

      ! diagonal =
      !   1   0   0   0   0
      !   0  -2   0   0   0
      !   0   0   3   0   0
      !   0   0   0  -4   0
      !   0   0   0   0   5
      d0 = reshape([1.0, -2.0, 3.0, -4.0, 5.0], [1,5])
      call diag%init(5, 0, d0=d0)

      ! tridiag = diag(1:4,1) + diag(5:9) + diag(10:13,-1)
      !   5    1    0    0    0
      !  10    6    2    0    0
      !   0   11    7    3    0
      !   0    0   12    8    4
      !   0    0    0   13    9
      d0      = reshape([(i, i=0,14)], [3,5], order=[2,1])
      d0(3,5) = 0
      call tri%init(5, 1, d0=d0)
      call tri%factorize

      ! penta = diag(7:10,1) + diag(11:15) + diag(16:19,-1) + diag(3:5,2) + diag(21:23,-2)
      !  11    7    3    0    0
      !  16   12    8    4    0
      !  21   17   13    9    5
      !   0   22   18   14   10
      !   0    0   23   19   15
      d0 = reshape([(i, i=1,25)], [5,5], order=[2,1])
      call penta%init(5, 2, d0=d0)
      call penta%factorize

      ! band = diag(7:10,1) + diag(11:15) + diag(16:19,-1) + diag(3:5,2)
      !  11    7    3    0    0
      !  16   12    8    4    0
      !   0   17   13    9    5
      !   0    0   18   14   10
      !   0    0    0   19   15
      d0 = reshape([(i, i=1,20)], [4,5], order=[2,1])
      call band%init(5, 1, nupper=2, d0=d0)
      call band%factorize
      d_band = reshape([11, 16,  0,  0,  0, &
        &                7, 12, 17,  0,  0, &
        &                3,  8, 13, 18,  0, &
        &                0,  4,  9, 14, 19, &
        &                0,  0,  5, 10, 15], [5,5])
    end block

    ! mul_vec
    block
      real :: x(5), y0(5), y(5), y_exp(5), fact = -3

      x  = [-1,2,5, -7, 11]
      y0 = [-10,2,50, -7, 11]

      ! diag
      y_exp = [29, -10, -135, 49, 22]
      y     = y0
      call diag%mul_vec(x, y, fact_y=fact)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: diag")

      ! tri
      y_exp = [27, 6, -114, 69, -25]
      y     = y0
      call tri%mul_vec(x, y, fact_y=fact)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: tridiag")

      ! penta
      y_exp = [48, 14, -80, 167, 114]
      y     = y0
      call penta%mul_vec(x, y, fact_y=fact)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: penta")

      ! band
      y_exp = [48, 14, -59, 123, -1]
      y     = y0
      call band%mul_vec(x, y, fact_y=fact)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: band")
    end block

    ! solve_vec
    block
      real    :: b(5), x(5), x_exp(5)
      integer :: i

      b     = [(i, i=3,7)]

      ! diag
      x_exp = [3.0000, -2.0000, 5.0/3.0, -1.5000, 1.4000]
      call diag%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: diag")

      ! tri
      x_exp = [6.551020408163266e-01,-2.755102040816327e-01,-4.489795918367347e-01,3.724489795918368e+00,-4.602040816326531e+00]
      call tri%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: tri")

      ! penta
      x_exp = [0.0,0.0,1.0,-1.0,0.2]
      call penta%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: penta")

      ! band
      x_exp = [1.694574681848627e-01,1.232417950435365e-01,9.109176155391831e-02,-2.297387809778969e-01,7.576691225720027e-01]
      call band%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: band")
    end block

    ! to dense
    block
      type(dense_real) :: dense
      real             :: mat_exp(7,7)

      call dense%init(7)
      call band%to_dense(dense, 2, 2)

      mat_exp          = 0
      mat_exp(2:6,2:6) = d_band

      call tc%assert_eq(mat_exp, dense%d, 1e-12, "to dense")
    end block

    ! to sparse
    block
      type(sparse_real)    :: s
      type(spbuild_real)   :: sb
      real, allocatable    :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call s%init(7)
      call sb%init(s)

      call band%to_sparse(sb, 2, 2)
      call sb%save

      ! s(2:6,2:6) =
      !  11    7    3    0    0
      !  16   12    8    4    0
      !   0   17   13    9    5
      !   0    0   18   14   10
      !   0    0    0   19   15
      a_exp  = [11, 7, 3, 16, 12, 8, 4, 17, 13, 9, 5, 18, 14, 10, 19, 15]
      ia_exp = [1, 1, 4, 8, 12, 15, 17, 17]
      ja_exp = [2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 5, 6]

      call tc%assert_eq(a_exp, s%a, 1e-12, "to sparse: a")
      call tc%assert_eq(ia_exp, int(s%ia), "to sparse: ia")
      call tc%assert_eq(ja_exp, s%ja, "to sparse: ja")
    end block

    ! set_elem
    block
      type(band_real)   :: m
      real, allocatable :: mat_exp(:,:)
      type(dense_real)  :: d

      ! test 1: nlower=nupper=1
      ! m = /1 2 0\
      !     |0 0 1|
      !     \0 3 4/
      call m%init(3, 1)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(3, 2, 3.0)
      call m%set_elem(3, 3, 4.0)

      call d%init(3)
      call m%to_dense(d)

      mat_exp = reshape([1,0,0,2,0,3,0,1,4], [3,3])

      call tc%assert_eq(mat_exp, d%d, 1e-12, "set_elem: nlower=nupper=1")

      ! test 2: nlower=0, nupper=2
      ! m = /1 2 0 0\
      !     |0 0 1 3|
      !     |0 0 1 0|
      !     \0 0 0 4/
      call m%init(4, 0, nupper=2)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(2, 4, 3.0)
      call m%set_elem(3, 3, 1.0)
      call m%set_elem(4, 4, 4.0)

      call d%init(4)
      call m%to_dense(d)

      mat_exp = reshape([1,0,0,0, 2,0,0,0, 0,1,1,0, 0,3,0,4], [4,4])

      call tc%assert_eq(mat_exp, d%d, 1e-12, "set_elem: nlower=0, nupper=2")
    end block

    ! get_elem
    block
      integer           :: i, j
      type(band_real)   :: m
      real, allocatable :: mat_exp(:,:), mat(:,:)
      type(dense_real)  :: d

      ! test 1: nlower=nupper=1
      ! m = /1 2 0\
      !     |0 0 1|
      !     \0 3 4/
      call m%init(3, 1)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(3, 2, 3.0)
      call m%set_elem(3, 3, 4.0)

      call d%init(3)
      call m%to_dense(d)

      mat_exp = reshape([1,0,0,2,0,3,0,1,4], [3,3])

      allocate(mat(3,3), source=0.0)
      do i = 1, 3
        do j = max(1,i-1),  min(3,i+1)
          mat(i,j) = m%get_elem(i, j)
        end do
      end do

      call tc%assert_eq(mat_exp, mat, 1e-12, "get_elem: nlower=nupper=1")

      ! test 2: nlower=0, nupper=2
      ! m = /1 2 0 0\
      !     |0 0 1 3|
      !     |0 0 1 0|
      !     \0 0 0 4/
      call m%init(4, 0, nupper=2)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(2, 4, 3.0)
      call m%set_elem(3, 3, 1.0)
      call m%set_elem(4, 4, 4.0)

      call d%init(4)
      call m%to_dense(d)

      mat_exp = reshape([1,0,0,0, 2,0,0,0, 0,1,1,0, 0,3,0,4], [4,4])

      deallocate(mat)
      allocate(mat(4,4), source=0.0)
      do i = 1, 4
        do j = i, min(4,i+2)
          mat(i,j) = m%get_elem(i, j)
        end do
      end do

      call tc%assert_eq(mat_exp, mat, 1e-12, "get_elem: nlower=0, nupper=2")
    end block

    ! add_elem
    block
      type(band_real)   :: m
      real, allocatable :: mat_exp(:,:)
      type(dense_real)  :: d

      ! test 1: nlower=nupper=1
      !     /1 2 0\   / 2 3 0 \   / 3 5 0 \
      ! m = |0 0 1| + | 4 0 5 | = | 4 0 6 |
      !     \0 3 4/   \ 0 6 1 /   \ 0 9 5 /
      call m%init(3, 1)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(3, 2, 3.0)
      call m%set_elem(3, 3, 4.0)

      call m%add_elem(1, 1, 2.0)
      call m%add_elem(1, 2, 3.0)
      call m%add_elem(2, 1, 4.0)
      call m%add_elem(2, 3, 5.0)
      call m%add_elem(3, 2, 6.0)
      call m%add_elem(3, 3, 1.0)

      call d%init(3)
      call m%to_dense(d)

      mat_exp = reshape([3,4,0, 5,0,9, 0,6,5], [3,3])

      call tc%assert_eq(mat_exp, d%d, 1e-12, "add_elem: nlower=nupper=1")

      ! test 2: nlower=0, nupper=2
      !     /1 2 0 0\   / 0 1 2 0 \   / 1 3 2 0 \
      ! m = |0 0 1 3| + | 0 3 0 4 | = | 0 3 1 7 |
      !     |0 0 1 0|   | 0 0 5 0 |   | 0 0 6 0 |
      !     \0 0 0 4/   \ 0 0 0 5 /   \ 0 0 0 9 /
      call m%init(4, 0, nupper=2)

      call m%set_elem(1, 1, 1.0)
      call m%set_elem(1, 2, 2.0)
      call m%set_elem(2, 3, 1.0)
      call m%set_elem(2, 4, 3.0)
      call m%set_elem(3, 3, 1.0)
      call m%set_elem(4, 4, 4.0)

      call m%add_elem(1, 2, 1.0)
      call m%add_elem(1, 3, 2.0)
      call m%add_elem(2, 2, 3.0)
      call m%add_elem(2, 4, 4.0)
      call m%add_elem(3, 3, 5.0)
      call m%add_elem(4, 4, 5.0)

      call d%init(4)
      call m%to_dense(d)

      mat_exp = reshape([1,0,0,0, 3,3,0,0, 2,1,6,0, 0,7,0,9], [4,4])

      call tc%assert_eq(mat_exp, d%d, 1e-12, "add_elem: nlower=0, nupper=2")
    end block

    ! reset_row
    block
      integer          :: i, j, k
      real             :: val(5), val_exp(5)
      type(band_real)  :: b0, b
      type(dense_real) :: d, d0

      call  d%init(5)
      call d0%init(5)

      do i = 1, 4             ! loop 4 basic band matrices
        select case(i)
          case (1)
            b0 = diag
          case (2)
            b0 = tri
          case (3)
            b0 = penta
          case (4)
            b0 = band
        end select

        call d0%reset
        call b0%to_dense(d0)                        ! save their dense representation for value comparisons

        do j = 1, 5                                 ! reset row j of b0
          b = b0                                    ! always start w/ fresh b0 s.t. we just reset row j
          call b%reset_row(j)

          call d%reset
          call b%to_dense(d)                        ! d will be dense representation of b0 where row j was resetted

          do k = 1, 5                               ! compare all rows k of b and b0 (take care to overwrite row j of b0)
            val = d%d(k,:)

            if (j == k) then
              val_exp = 0
            else
              val_exp = d0%d(k,:)
            end if

            call tc%assert_eq(val_exp, val, 1e-12, "reset_row")
          end do
        end do
      end do
    end block

    ! set diag
    block
      integer :: i
      type(band_real) :: b, b_exp

      ! test 1: build diagonal diag
      call b%init(5, 0)
      call b%set_diag(real([1, -2, 3, -4, 5]))
      call tc%assert_eq(diag%d, b%d, 1e-12, "set_diag_arr: build diag")

      ! test 2: build tridiagonal tri
      call b%init(5, 1)
      call b%set_diag([(real(i), i= 1,  4)], k=-1)
      call b%set_diag([(real(i), i= 5,  9)]      )
      call b%set_diag([(real(i), i=10, 13)], k= 1)
      call tc%assert_eq(tri%d, b%d, 1e-12, "set_diag_arr: build tri")

      ! test 3: build
      ! [1 2 0]
      ! [4 1 3]
      ! [5 4 1]
      call b%init(    3, 2, nupper=1)
      call b_exp%init(3, 2, nupper=1, d0=reshape(real([0,2,3, 1,1,1, 4,4,0, 5,0,0]), [4,3], order=[2,1]))

      call b%set_diag([2.0, 3.0], k=-1)
      call b%set_diag(1.0)
      call b%set_diag(4.0, k=1)
      call b%set_diag(5.0, k=2)
      call tc%assert_eq(b_exp%d, b%d, 1e-12, "set_diag_{arr,val}: test 3")
    end block

    call tc%finish
  end subroutine

end submodule
