module test_map_m

  use map_m
  use string_m,    only: string
  use test_case_m, only: test_case
  use util_m,      only: int2str

  implicit none

  private
  public test_map

contains

  subroutine test_map()
    integer                            :: i
    type(map_string_real)              :: map
    type(mapnode_string_real), pointer :: node
    real                               :: r(5),   r_arr(5)
    type(string)                       :: str(5), str_arr(5)
    type(test_case)                    :: tc

    call tc%init("map")

    r = [-1, 0, 1, 2, 3]
    str(1)%s = "a"
    str(2)%s = "b"
    str(3)%s = "c"
    str(4)%s = "d"
    str(5)%s = "e"

    call map%init()

    do i = 1, 5
      call map%insert(str(i), r(i))
    end do

    ! test copy
    block
      type(map_string_real) :: map_copy

      call map_copy%copy(map)
      call tc%assert_eq(r(2), map_copy%get(str(2)), 0.0, "copy")
    end block


    ! test get
    do i = 1, 5
      call tc%assert_eq(r(i), map%get(str(i)), 0.0, "get "//int2str(i))
    end do

    ! test find
    do i = 1, 5
      node => map%find(str(i))
      call tc%assert_eq(r(i), node%value, 0.0, "find "//int2str(i))
    end do

    ! test to array
    call map%to_array(keys = str_arr, values = r_arr)
    call tc%assert_eq(r,   r_arr,   0.0,  "to array: values")
    call tc%assert_eq(str, str_arr,       "to array: keys")

    ! test set
    call map%set(string("f"), 5.0)
    call tc%assert_eq(5.0, map%get(string("f")), 0.0, "set new")
    call map%set(string("f"), 6.0)
    call tc%assert_eq(6.0, map%get(string("f")), 0.0, "set existing")

    ! test destruct
    call map%destruct()
    node => map%find(str(1))
    call tc%assert(.not. associated(node), "destruct")

    call tc%finish()
  end subroutine

end module
