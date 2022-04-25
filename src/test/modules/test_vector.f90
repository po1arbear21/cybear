module test_vector_m

  use test_case_m, only: test_case
  use vector_m,    only: vector_int

  implicit none

  private
  public test_vector

contains

  subroutine test_vector()
    type(test_case)  :: tc
    type(vector_int) :: vec

    call tc%init("vector")

    ! init
    call vec%init(0, c = 4)
    call tc%assert_eq(0, vec%n, "init: n")
    call tc%assert(allocated(vec%d), "init: allocated(d)")
    call tc%assert_eq(4, size(vec%d), "init: capacity")

    ! init with initial values
    call vec%init(3, c = 7, x = [1, 2, 3])
    call tc%assert_eq(3, vec%n, "init_x: n")
    call tc%assert(allocated(vec%d), "init_x: allocated(d)")
    call tc%assert_eq(7, size(vec%d), "init_x: capacity")
    call tc%assert_eq([1, 2, 3], vec%d(1:vec%n), "init_x: values")

    ! reset
    call vec%init(2, c = 3)
    call vec%reset()
    call tc%assert_eq(0, vec%n, "reset: n")

    ! reserve
    call vec%reserve(50)
    call tc%assert_eq(50, size(vec%d), "reserve: capacity")

    ! resize
    call vec%resize(10)
    call tc%assert_eq(10, vec%n, "resize: n")

    ! front/back
    call vec%init(3, c = 4, x = [1, 2, 3])
    call tc%assert_eq(1, vec%front(), "front")
    call tc%assert_eq(3, vec%back(), "back")

    ! push_elem
    call vec%init(0, c = 1)
    call vec%push(1)
    call vec%push(2)
    call vec%push(3)
    call vec%push(4)
    call vec%push(5)
    call vec%push(6)
    call vec%push(7)
    call tc%assert_eq(7, vec%n, "push_elem: n")
    call tc%assert_eq(1, vec%front(), "push_elem: front")
    call tc%assert_eq(7, vec%back(), "push_elem: back")

    ! push_elems
    call vec%push([8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    call tc%assert_eq(20, vec%n, "push_elems: n")
    call tc%assert_eq([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], vec%d(1:vec%n), "push_elems: d")

    ! pop
    call vec%pop
    call tc%assert_eq(19, vec%n, "pop")
    call tc%assert_eq(19, vec%back(), "pop: back")

    ! shrink
    call vec%shrink()
    call tc%assert_eq(19, size(vec%d), "shrink: capacity")
    call tc%assert_eq([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], vec%d(1:vec%n), "shrink: d")

    ! to_array
    call tc%assert_eq([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], vec%to_array(), "to_array")

    ! from_array
    call vec%resize(3)
    call vec%shrink()
    call vec%from_array([1, 2, 3, 4, 5])
    call tc%assert_eq(5, vec%n, "from_array: n")
    call tc%assert_eq([1, 2, 3, 4, 5], vec%to_array(), "from_array: d")

    call tc%finish
  end subroutine

end module
