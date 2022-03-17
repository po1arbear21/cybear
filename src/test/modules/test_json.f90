module test_json_m

  use json_m
  use string_m,    only: string, new_string
  use test_case_m, only: test_case

  implicit none

  private
  public test_json

contains

  subroutine test_json()
    type(test_case)           :: tc

    call tc%init("json")

    ! json_array
    block
      integer                    :: i
      logical                    :: l
      real                       :: r
      type(string)               :: s1, s2
      type(json_array)           :: ar
      type(json_array),  pointer :: ar2
      type(json_object), pointer :: obj
      class(json),       pointer :: js

      call ar%init()
      call ar%add()
      call ar%add(.true.)
      call ar%add(.false.)
      call ar%add(2)
      call ar%add(3.0)
      call ar%add("string1")
      call ar%add("string2")
      allocate (ar2)
      call ar2%init()
      call ar2%add(1)
      call ar2%add(2)
      call ar2%add(3)
      call ar%add(ar2)
      allocate (obj)
      call obj%init()
      call ar%add(obj)

      call ar%get_json(1, js)
      select type (js)
      type is (json_null)
        call tc%assert(.true., "json_array_get => null")
      class default
        call tc%assert(.false., "json_array_get => null")
      end select

      call ar%get_json(2, js)
      select type (js)
      type is (json_bool)
        call tc%assert(js%value, "json_array_get => bool")
      class default
        call tc%assert(.false., "json_array_get => bool")
      end select

      call ar%get(3, l)
      call tc%assert(.not. l, "json_array_get_bool")

      call ar%get(4, i)
      call tc%assert_eq(2, i, "json_array_get_int")

      call ar%get(5, r)
      call tc%assert_eq(3.0, r, 0.0, "json_array_get_real")

      call ar%get(6, s1%s)
      call tc%assert_eq(new_string("string1"), s1, "json_array_get_string 1")

      call ar%get(7, s2%s)
      call tc%assert_eq(new_string("string2"), s2, "json_array_get_string 2")

      call ar%get(8, ar2)
      call ar2%get(1, i)
      call tc%assert_eq(1, i, "json_array_get_array 1")
      call ar2%get(2, i)
      call tc%assert_eq(2, i, "json_array_get_array 2")
      call ar2%get(3, i)
      call tc%assert_eq(3, i, "json_array_get_array 3")

      call ar%get(9, obj)
      call tc%assert_eq(0, obj%names%n,      "json_array_get_object 1")
      call tc%assert_eq(0, obj%properties%n, "json_array_get_object 2")

      call ar%set(1, 5)
      call ar%set(2, 6)
      call ar%get(1, i)
      call tc%assert_eq(5, i, "json_array_set_int 1")
      call ar%get(2, i)
      call tc%assert_eq(6, i, "json_array_set_int 2")

      call ar%destruct()
    end block

    ! json_object
    block
      integer                    :: i
      real                       :: r
      type(string)               :: s1
      type(json_array),  pointer :: ar
      type(json_object)          :: obj
      type(json_object), pointer :: obj2
      class(json),       pointer :: js

      call obj%init()
      call obj%add("name0")
      call obj%add("name1", .true.)
      call obj%add("name2", 2)
      call obj%add("name3", 3.0)
      call obj%add("name4", "value")
      allocate (ar)
      call ar%init()
      call ar%add(1)
      call ar%add(2)
      call ar%add(3)
      call obj%add("name5", ar)
      allocate (obj2)
      call obj2%init()
      call obj2%add("a", 1)
      call obj2%add("b", 2)
      call obj2%add("c", 3)
      call obj%add("name6", obj2)

      call obj%get_json("name0", js)
      select type (js)
      type is (json_null)
        call tc%assert(.true., "json_object_get => json_null")
      class default
        call tc%assert(.false., "json_object_get => json_null")
      end select

      call obj%get_json("name1", js)
      select type (js)
      type is (json_bool)
        call tc%assert(js%value, "json_object_get => json_bool")
      class default
        call tc%assert(.false., "json_object_get => json_bool")
      end select

      call obj%get("name2", i)
      call tc%assert_eq(2, i, "json_object_get_int")

      call obj%get("name3", r)
      call tc%assert_eq(3.0, r, 0.0, "json_object_get_real")

      call obj%get("name4", s1%s)
      call tc%assert_eq(new_string("value"), s1, "json_object_get_string")

      call obj%get("name5", ar)
      call ar%get(1, i)
      call tc%assert_eq(1, i, "json_object_get_array 1")
      call ar%get(2, i)
      call tc%assert_eq(2, i, "json_object_get_array 2")
      call ar%get(3, i)
      call tc%assert_eq(3, i, "json_object_get_array 3")

      call obj%get("name6", obj2)
      call obj2%get("a", i)
      call tc%assert_eq(1, i, "json_object_get_object 1")
      call obj2%get("b", i)
      call tc%assert_eq(2, i, "json_object_get_object 1")
      call obj2%get("c", i)
      call tc%assert_eq(3, i, "json_object_get_object 1")

      call obj2%set("d", 4)
      call obj2%set("a", 5)
      call obj2%get("d", i)
      call tc%assert_eq(4, i, "json_object_set_int 1")
      call obj2%get("a", i)
      call tc%assert_eq(5, i, "json_object_set_int 2")

      call obj%destruct()
    end block

    ! json_load, json_save
    block
      integer               :: funit1, funit2, n1, n2, iostat
      character(80)         :: iomsg
      class(json),  pointer :: js
      type(string)          :: s1, s2

      call json_load(js, "src/test/test1.json")
      call json_save(js, "src/test/test3.json")
      call js%destruct()
      deallocate (js)

      ! try to open saved file
      open (newunit = funit1, file = "src/test/test3.json", access = "stream", form = "unformatted", status = "old", action = "read", iostat = iostat, iomsg = iomsg)
      call tc%assert_eq(0, iostat, "json_save 1")
      if (.not. tc%last_passed) goto 100

      open (newunit = funit2, file = "src/test/test2.json", access = "stream", form = "unformatted", status = "old", action = "read", iostat = iostat, iomsg = iomsg)
      call tc%assert_eq(0, iostat, "json_save 2")
      if (.not. tc%last_passed) goto 100

      ! compare file sizes
      inquire (file = "src/test/test3.json", size = n1)
      inquire (file = "src/test/test2.json", size = n2)
      call tc%assert_eq(n2, n1, "json_save 3")
      if (.not. tc%last_passed) goto 100

      allocate (character(n1) :: s1%s)
      read (funit1, iostat = iostat, iomsg = iomsg) s1%s
      close (funit1)
      call tc%assert_eq(0, iostat, "json_save 4")
      if (.not. tc%last_passed) goto 100

      allocate (character(n2) :: s2%s)
      read (funit2, iostat = iostat, iomsg = iomsg) s2%s
      close (funit2)
      call tc%assert_eq(0, iostat, "json_save 5")
      if (.not. tc%last_passed) goto 100

      call tc%assert_eq(s2, s1, "json_save 6")
    end block

    100 call tc%finish()
  end subroutine

end module
