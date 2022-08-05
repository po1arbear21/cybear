module test_json_m

  use iso_fortran_env, only: int32, int64
  use json_m
  use string_m,        only: string, new_string
  use test_case_m,     only: test_case

  implicit none

  private
  public test_json

contains

  subroutine test_json()
    type(test_case)           :: tc

    call tc%init("json")

    ! json_array
    block
      integer(kind=int64)        :: i
      logical                    :: l
      real                       :: r
      type(string)               :: s1, s2
      type(json_array)           :: ar
      type(json_array),  pointer :: ar2
      type(json_object), pointer :: obj
      class(json),       pointer :: js

      call ar%init()
      call ar%add_null()
      call ar%add_log(.true.)
      call ar%add_log(.false.)
      call ar%add_int(2)
      call ar%add_real(3.0)
      call ar%add_string("string1")
      call ar%add_string("string2")
      call ar%add_array(p = ar2)
      call ar2%add_int(1)
      call ar2%add_int(2)
      call ar2%add_int(3)
      call ar%add_object(p = obj)

      js => ar%values%d(1)%p
      select type (js)
      type is (json_null)
        call tc%assert(.true., "json_array => null")
      class default
        call tc%assert(.false., "json_array => null")
      end select

      js => ar%values%d(2)%p
      select type (js)
      type is (json_log)
        call tc%assert(js%value, "json_array => logical")
      class default
        call tc%assert(.false., "json_array => logical")
      end select

      l = ar%get_log(3)
      call tc%assert(.not. l, "json_array_get_log")

      i = ar%get_int(4)
      call tc%assert_eq(2, int(i), "json_array_get_int")

      r = ar%get_real(5)
      call tc%assert_eq(3.0, r, 0.0, "json_array_get_real")

      s1%s = ar%get_string(6)
      call tc%assert_eq(new_string("string1"), s1, "json_array_get_string 1")

      s2%s = ar%get_string(7)
      call tc%assert_eq(new_string("string2"), s2, "json_array_get_string 2")

      ar2 => ar%get_array(8)
      i = ar2%get_int(1)
      call tc%assert_eq(1, int(i), "json_array_get_array 1")
      i = ar2%get_int(2)
      call tc%assert_eq(2, int(i), "json_array_get_array 2")
      i = ar2%get_int(3)
      call tc%assert_eq(3, int(i), "json_array_get_array 3")

      obj => ar%get_object(9)
      call tc%assert_eq(0, obj%names%n,      "json_array_get_object 1")
      call tc%assert_eq(0, obj%properties%n, "json_array_get_object 2")

      call ar%destruct()
    end block

    ! json_object
    block
      integer(kind=int64)        :: i
      real                       :: r
      type(string)               :: s1
      type(json_array),  pointer :: ar
      type(json_object)          :: obj
      type(json_object), pointer :: obj2
      class(json),       pointer :: js

      call obj%init()
      call obj%add_null("name0")
      call obj%add_log("name1", .true.)
      call obj%add_int("name2", 2)
      call obj%add_real("name3", 3.0)
      call obj%add_string("name4", "value")
      call obj%add_array("name5", p = ar)
      call ar%add_int(1)
      call ar%add_int(2)
      call ar%add_int(3)
      call obj%add_object("name6", p = obj2)
      call obj2%add_int("a", 1)
      call obj2%add_int("b", 2)
      call obj2%add_int("c", 3)

      js => obj%get_json("name0")
      call tc%assert(associated(js), "json_object_get_json 1")
      select type (js)
      type is (json_null)
        call tc%assert(.true., "json_object_get_json 1 => json_null")
      class default
        call tc%assert(.false., "json_object_get_json 1 => json_null")
      end select

      js => obj%get_json("name1")
      call tc%assert(associated(js), "json_object_get_json 2")
      select type (js)
      type is (json_log)
        call tc%assert(js%value, "json_object_get_json 2 => json_log")
      class default
        call tc%assert(.false., "json_object_get_json 2 => json_log")
      end select

      i = obj%get_int("name2")
      call tc%assert_eq(2, int(i), "json_object_get_int")

      r = obj%get_real("name3")
      call tc%assert_eq(3.0, r, 0.0, "json_object_get_real")

      s1%s = obj%get_string("name4")
      call tc%assert_eq(new_string("value"), s1, "json_object_get_string")

      ar => obj%get_array("name5")
      i = ar%get_int(1)
      call tc%assert_eq(1, int(i), "json_object_get_array 1")
      i = ar%get_int(2)
      call tc%assert_eq(2, int(i), "json_object_get_array 2")
      i = ar%get_int(3)
      call tc%assert_eq(3, int(i), "json_object_get_array 3")

      obj2 => obj%get_object("name6")
      i = obj2%get_int("a")
      call tc%assert_eq(1, int(i), "json_object_get_object 1")
      i = obj2%get_int("b")
      call tc%assert_eq(2, int(i), "json_object_get_object 1")
      i = obj2%get_int("c")
      call tc%assert_eq(3, int(i), "json_object_get_object 1")

      call obj2%set_int("d", 4)
      call obj2%set_int("a", 5)
      i = obj2%get_int("d")
      call tc%assert_eq(4, int(i), "json_object_set_int 1")
      i = obj2%get_int("a")
      call tc%assert_eq(5, int(i), "json_object_set_int 2")

      call obj%destruct()
    end block

    ! json load and save
    block
      integer         :: funit1, funit2, n1, n2, iostat
      character(80)   :: iomsg
      type(json_file) :: jsfile
      type(string)    :: s1, s2

      call jsfile%load("src/test/test1.json")
      call jsfile%save("src/test/test3.json")
      call jsfile%destruct()

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
