submodule(test_matrix_m) test_sparse_m
  use matrix_m
  implicit none

contains
  module subroutine test_sparse()
    type(test_case) :: tc

    print "(1A)", "test_sparse"
    call tc%init("sparse")

    ! sparse builder
    block
      integer            :: i, j
      type(sparse_real)  :: s
      type(spbuild_real) :: sb

      call s%init(5, 5)
      call sb%init(s)

      call s%init(5)
      call sb%init(s)

      call sb%add(1, 1, 2.5)
      call sb%add(1, 2, 1.7)
      call sb%add(2, 3, 3.2)
      call sb%add(3, 4, 1.1)
      call sb%add(4, 4, 0.1)
      call sb%add(5, 5, 9.0)
      call sb%add(5, 1, 6.7)
      call sb%add(5, 2, 3.8)

      call sb%save()

      do i = 1, s%nrows
        do j = s%ia(i), s%ia(i+1)-1
          print *, i, s%ja(j), s%a(j)
        end do
      end do
    end block

    ! sparse matrix matrix multiplication
    block
      integer            :: i, j
      type(sparse_real)  :: s1, s2, s3
      type(spbuild_real) :: sb

      call s1%init(4, ncols = 3)
      call sb%init(s1)
      call sb%add(1, 1, 7.0)
      call sb%add(1, 3, 1.0)
      call sb%add(2, 2, 3.0)
      call sb%add(2, 3, 4.0)
      call sb%add(3, 1, 5.0)
      call sb%add(4, 3, 2.0)
      call sb%save()

      call s2%init(3, ncols = 5)
      call sb%init(s2)
      call sb%add(1, 3, 4.0)
      call sb%add(2, 1, 8.0)
      call sb%add(2, 4, 4.0)
      call sb%add(3, 1, 2.0)
      call sb%add(3, 3, 3.0)
      call sb%save()

      call s1%mul_sparse(s2, s3)

      print *
      print *
      print *
      do i = 1, s3%nrows
        do j = s3%ia(i), s3%ia(i+1)-1
          print *, i, s3%ja(j), s3%a(j)
        end do
      end do
    end block

    call tc%finish()
  end subroutine
end submodule