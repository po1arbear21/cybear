module test_case_m
  use error_m
  use string_m
  implicit none

  private
  public :: test_case

  type :: test_case
    !! Test case type

    character(:), allocatable :: name
      !! test name to append to err messages
    integer :: num_tests
      !! count number of tests
    integer :: passed_tests
      !! count how many tests passed
    logical :: last_passed
      !! latest test result
  contains
    ! public procedures
    procedure :: init
    procedure :: finish
    generic   :: assert_eq => assert_eq_int   , assert_eq_int_arr   , assert_eq_int_arr2D   , &
                              assert_eq_real  , assert_eq_real_arr  , assert_eq_real_arr2D  , &
                              assert_eq_string, assert_eq_string_arr, assert_eq_string_arr2D

    ! private procedures
    procedure, private :: assert_eq_int
    procedure, private :: assert_eq_int_arr
    procedure, private :: assert_eq_int_arr2D
    procedure, private :: assert_eq_real
    procedure, private :: assert_eq_real_arr
    procedure, private :: assert_eq_real_arr2D
    procedure, private :: assert_eq_string
    procedure, private :: assert_eq_string_arr
    procedure, private :: assert_eq_string_arr2D
    procedure, private :: check
  end type

contains

  subroutine init(this, name)
    !! Initialize test case.

    class(test_case), intent(out)           :: this
      !! Test case
    character(*),     intent(in)            :: name
      !! Name to append to error messages for this test

    this%name         = name
    this%num_tests    = 0
    this%passed_tests = 0
    this%last_passed  = .true.
  end subroutine

  subroutine finish(this)
    !! Finish test and print number of passed tests.

    class(test_case), intent(in) :: this
      !! Test case

    print "(1A, 1I0, 1A, 1I0, 1A)", "  "//trim(this%name)//": ", this%passed_tests, "/", this%num_tests, " tests passed!"
    print *

    if (this%passed_tests /= this%num_tests) call program_error("Test case "//trim(this%name)//" failed!")
  end subroutine

  subroutine assert_eq_int(this, expected, value, msg)
    !! Check integer value (scalar) for equality to expected value.

    class(test_case), intent(inout) :: this
      !! Test case
    integer,          intent(in)    :: expected
      !! Expected value
    integer,          intent(in)    :: value
      !! Actual value
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    call this%check((value == expected), msg, &
                    msg2 = "scalar values differs", ipar = [expected, value])
  end subroutine

  subroutine assert_eq_int_arr(this, expected, values, msg)
    !! Check integer values (1D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    integer,          intent(in)    :: expected(:)
      !! Expected values
    integer,          intent(in)    :: values(:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, failed

    ! check array sizes
    call this%check((size(expected) == size(values)), msg, &
                    msg2 = "array lengths differ")
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check((values(i) == expected(i)), msg, &
                      msg2 = "array values differ; (row, expected, value)", ipar = [i, expected(i), values(i)])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(1A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_eq_int_arr2D(this, expected, values, msg)
    !! Check integer values (2D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    integer,          intent(in)    :: expected(:,:)
      !! Expected values
    integer,          intent(in)    :: values(:,:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, j, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        call this%check((values(i,j) == expected(i,j)), msg, &
                        msg2 = "array values differ; (row, col, expected, value)", ipar = [i, j, expected(i,j), values(i,j)])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(1A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_real(this, expected, value, abs_tol, msg)
    !! Check real value (scalar) for equality to expected value within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected
      !! Expected value
    real,             intent(in)    :: value
      !! Actual value
    real,             intent(in)    :: abs_tol
      !! Absolute tolerance to compare value to
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    call this%check(abs(expected - value) <= abs_tol, msg, &
                    msg2 = "scalar value differs", rpar = [expected, value])
  end subroutine

  subroutine assert_eq_real_arr(this, expected, values, abs_tol, msg)
    !! Check real values (1D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected(:)
      !! Expected values
    real,             intent(in)    :: values(:)
      !! Actual values
    real,             intent(in)    :: abs_tol
      !! Absolute tolerance to compare values to
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, failed

    ! check array sizes
    call this%check((size(expected) == size(values)), msg, &
                    msg2 = "array lengths differ")
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check(abs(expected(i) - values(i)) <= abs_tol, msg, &
                      msg2 = "array values differ; (row; expected, value)", ipar = [i], rpar = [expected(i), values(i)])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(1A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_eq_real_arr2D(this, expected, values, abs_tol, msg)
    !! Check real values (2D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected(:,:)
      !! Expected values
    real,             intent(in)    :: values(:,:)
      !! Actual values
    real,             intent(in)    :: abs_tol
      !! Absolute tolerance to compare values to
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, j, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        call this%check(abs(expected(i,j) - values(i,j)) <= abs_tol, msg, &
                        msg2 = "array values differ; (row, col; expected, value)", ipar = [i, j], rpar = [expected(i,j), values(i,j)])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(1A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_string(this, expected, value, msg)
    !! Check string (scalar) for equality to expected value.

    class(test_case), intent(inout) :: this
      !! Test case
    type(string),     intent(in)    :: expected
      !! Expected value
    type(string),     intent(in)    :: value
      !! Actual value
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    call this%check((value%s == expected%s), msg, &
                    msg2 = "scalar strings differs")
  end subroutine

  subroutine assert_eq_string_arr(this, expected, values, msg)
    !! Check strings (1D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    type(string),     intent(in)    :: expected(:)
      !! Expected values
    type(string),     intent(in)    :: values(:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, failed

    ! check array sizes
    call this%check((size(expected) == size(values)), msg, &
                    msg2 = "array lengths differ")
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check((values(i)%s == expected(i)%s), msg, &
                      msg2 = "array values differ; (row)", ipar = [i])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(1A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_eq_string_arr2D(this, expected, values, msg)
    !! Check strings (2D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    type(string),     intent(in)    :: expected(:,:)
      !! Expected values
    type(string),     intent(in)    :: values(:,:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    ! local variables
    integer :: i, j, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        call this%check((values(i,j)%s == expected(i,j)%s), msg, &
                        msg2 = "array values differ; (row, col)", ipar = [i, j])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(1A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine check(this, chk, msg, msg2, ipar, rpar)
    !! Perform check.

    class(test_case),       intent(inout) :: this
      !! Test case
    logical,                intent(in)    :: chk
      !! True means test will pass
    character(*),           intent(in)    :: msg
      !! Message to print if test fails
    character(*), optional, intent(in)    :: msg2
      !! 2nd part of message (new line)
    integer,      optional, intent(in)    :: ipar(:)
      !! Additional integer parameters to print (new line)
    real,         optional, intent(in)    :: rpar(:)
      !! Additional real parameters to print (new line)

    ! local variables
    integer :: i

    ! count total number of tests
    this%num_tests = this%num_tests + 1

    if (chk) then
      this%last_passed  = .true.
      this%passed_tests = this%passed_tests + 1
    else
      this%last_passed  = .false.
      print "(1A)", "Test failed in "//trim(this%name)
      print "(1A)", "  "//trim(msg)
      if (present(msg2)) print "(1A)", "    "//trim(msg2)
      if (present(ipar)) then
        write(*, fmt="(1A)", advance="no") "      "
        do i = 1, size(ipar)
          if (i < size(ipar)) then
            write(*, fmt="(1I0, 1A)", advance="no") ipar(i), ", "
          else
            write(*, fmt="(1I0)") ipar(i)
          end if
        end do
        print *
      end if
      if (present(rpar)) then
        write(*, fmt="(1A)", advance="no") "      "
        do i = 1, size(rpar)
          if (i < size(rpar)) then
            write(*, fmt="(1E18.8, 1A)", advance="no") rpar(i), ", "
          else
            write(*, fmt="(1E18.8)") rpar(i)
          end if
        end do
        print *
      end if
    end if
  end subroutine

end module