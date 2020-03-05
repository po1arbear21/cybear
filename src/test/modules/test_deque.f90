module test_deque_m
  use test_case_m
  use deque_m
  implicit none

contains

  subroutine test_deque()
    type(test_case) :: tc
    type(deque_int) :: deq

    print "(1A)", "test_deque"
    call tc%init("deque")

    ! test1: init
    call deq%init(0, c = 4)
    call tc%assert_eq(0, deq%n, "init: n")
    call tc%assert(allocated(deq%d), "init: allocated(d)")
    call tc%assert_eq(4, size(deq%d), "init: capacity")

    ! test2: init with initial values
    call deq%init(3, c = 7, x = [1, 2, 3])
    call tc%assert_eq(3, deq%n, "init_x: n")
    call tc%assert(allocated(deq%d), "init_x: allocated(d)")
    call tc%assert_eq(7, size(deq%d), "init_x: capacity")
    call tc%assert_eq([1, 2, 3], deq%d(1:deq%n), "init_x: values")

    ! test3: destruct
    call deq%destruct()
    call tc%assert(.not. allocated(deq%d), "destruct")

    ! test4: reset
    call deq%init(2, c = 3)
    call deq%reset()
    call tc%assert_eq(0, deq%n, "reset: n")

    ! test5: reserve
    call deq%reserve(50)
    call tc%assert_eq(50, size(deq%d), "reserve: capacity")

    ! test6: resize
    call deq%resize(10)
    call tc%assert_eq(10, deq%n, "resize: n")

    ! test7: front/back
    call deq%init(3, c = 4, x = [1, 2, 3])
    call tc%assert_eq(1, deq%front(), "front")
    call tc%assert_eq(3, deq%back(), "back")

    ! test8: push_front
    call deq%init(0, c = 1)
    call deq%push_front(1)
    call deq%push_front(2)
    call deq%push_front(3)
    call deq%push_front(4)
    call deq%push_front(5)
    call tc%assert_eq(5, deq%n, "push_front: n")
    call tc%assert_eq(5, deq%front(), "push_front: front")
    call tc%assert_eq(1, deq%back(), "push_front: back")

    ! test9: push_back
    call deq%push_back(101)
    call deq%push_back(102)
    call deq%push_back(103)
    call deq%push_back(104)
    call deq%push_back(105)
    call deq%push_back(106)
    call deq%push_back(107)
    call deq%push_back(108)
    call tc%assert_eq(13, deq%n, "push_back: n")
    call tc%assert_eq(5, deq%front(), "push_back: front")
    call tc%assert_eq(108, deq%back(), "push_back: back")

    ! test10: pop_front
    call deq%pop_front()
    call deq%pop_front()
    call deq%pop_front()
    call deq%pop_front()
    call deq%pop_front()
    call deq%pop_front()
    call deq%push_back(109)
    call tc%assert_eq(8, deq%n, "pop_front: n")
    call tc%assert_eq(102, deq%front(), "pop_front: front")
    call tc%assert_eq(109, deq%back(), "pop_front: back")

    ! test11: shrink
    call deq%shrink()
    call tc%assert_eq(8, size(deq%d), "shrink: capacity")
    call tc%assert_eq([102, 103, 104, 105, 106, 107, 108, 109], deq%d(1:deq%n), "shrink: d")

    ! test12: to_array
    call tc%assert_eq([102, 103, 104, 105, 106, 107, 108, 109], deq%to_array(), "to_array")

    ! test13: from_array
    call deq%resize(3)
    call deq%shrink()
    call deq%from_array([1, 2, 3, 4, 5])
    call tc%assert_eq(5, deq%n, "from_array: n")
    call tc%assert_eq([1, 2, 3, 4, 5], deq%to_array(), "from_array: d")

    call tc%finish()
  end subroutine

end module
