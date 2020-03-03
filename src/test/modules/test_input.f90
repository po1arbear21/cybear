module test_input_m
  use input_m
  use test_case_m

  implicit none

contains

  subroutine test_input
    type(test_case)  :: tc
    type(input_file) :: f

    print "(1A)", "test_input"
    call tc%init("input")

    ! init input file
    call f%init("src/test/example.inp")

    ! reals (+ check normalization implicitly)
    block
      integer              :: i
      real                 :: r
      real,    allocatable :: r_arr(:)
      integer, allocatable :: sid(:)

      ! check get_name_real + init normalization
      call f%get("", "temperature", r, normalize=.false.)
      call tc%assert_eq(300.0, r, 1e-10, "get_name_real")
      call init_normconst(r)

      ! check get_name_real_arr
      call f%get("grid", "y", r_arr)
      call tc%assert_eq(norm([(real(i)/10.0, i=0,9)], 'um'), r_arr, 1e-10, "get_name_real_arr")

      ! get_real, changing sections
      call f%get_sections("region", sid)
      call f%get(sid(1), "y0", r)
      call tc%assert_eq(0.0, r, 1e-16, "get_real")
      call f%get(sid(2), "y0", r)
      call tc%assert_eq(norm(6e2, 'um'), r, 1e-16, "get_real")

      ! get_real_arr
      call f%get_sections("grid", sid)
      call f%get(sid(1), "y", r_arr)
      call tc%assert_eq(norm([(real(i)/10.0, i=0,9)], 'um'), r_arr, 1e-10, "get_real_arr")
    end block

    ! integers
    block
      integer              :: i
      integer, allocatable :: i_arr(:)
      integer, allocatable :: sid(:)

      ! check get_name_int
      call f%get("grid", "Nx", i)
      call tc%assert_eq(100, i, "get_name_int")

      ! check get_name_int_arr
      call f%get("grid", "idx", i_arr)
      call tc%assert_eq([2, 4, 3, 2, 2, -1], i_arr, "get_name_int_arr")

      ! check get_int
      call f%get_sections("grid", sid)
      call f%get(sid(1), "Nx", i)
      call tc%assert_eq(100, i, "get_int")

      ! check get_int_arr
      call f%get_sections("grid", sid)
      call f%get(sid(1), "idx", i_arr)
      call tc%assert_eq([2, 4, 3, 2, 2, -1], i_arr, "get_int_arr")
    end block

    ! strings
    block
      character(:), allocatable :: c_arr
      type(string), allocatable :: s_arr(:)
      integer,      allocatable :: sid(:)

      ! get_name_string
      call f%get("contact", "name", c_arr)
      call tc%assert_eq(string("GAT"), string(c_arr), "get_name_string")

      ! get_name_string_arr
      call f%get("contact", "names", s_arr)
      call tc%assert_eq([string("123"), string("asd asd"), string("asd 123 x")], s_arr, "get_name_string_arr")

      ! get_string
      call f%get_sections("contact", sid)
      call f%get(sid(1), "name", c_arr)
      call tc%assert_eq(string("GAT"), string(c_arr), "get_string")

      ! get_string_arr
      call f%get_sections("contact", sid)
      call f%get(sid(1), "names", s_arr)
      call tc%assert_eq([string("123"), string("asd asd"), string("asd 123 x")], s_arr, "get_string_arr")
    end block

    ! logicals
    block
      logical              :: tf
      logical, allocatable :: tf_arr(:)
      integer, allocatable :: sid(:)

      ! check get_name_logical
      call f%get("contact", "gate", tf)
      call tc%assert_eq(count([.true.]), count([tf]), "get_name_logical")

      ! check get_name_logical_arr
      call f%get("contact", "tf", tf_arr)
      call tc%assert_eq(count(tf_arr), count([.true., .true., .false., .true.]), "get_name_logical_arr")

      ! check get_logical
      call f%get_sections("contact", sid)
      call f%get(sid(1), "gate", tf)
      call tc%assert_eq(count([.true.]), count([tf]), "get_logical")

      ! check get_logical_arr
      call f%get_sections("contact", sid)
      call f%get(sid(1), "tf", tf_arr)
      call tc%assert_eq(count(tf_arr), count([.true., .true., .false., .true.]), "get_logical_arr")
    end block

    call tc%finish
  end subroutine

end module
