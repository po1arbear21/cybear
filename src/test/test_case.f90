module test_case_m

  use color_m,  only: COL_DEFAULT, COL_GREEN, COL_MAGENTA, COL_RED, COL_WHITE, COL_YELLOW
  use error_m,  only: program_error
  use string_m, only: string

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
    generic   :: assert    => assert_1, assert_arr, assert_arr2D, assert_arr3D
    generic   :: assert_eq => assert_eq_int      , assert_eq_real      , assert_eq_cmplx      , assert_eq_string      , &
                              assert_eq_int_arr  , assert_eq_real_arr  , assert_eq_cmplx_arr  , assert_eq_string_arr  , &
                              assert_eq_int_arr2D, assert_eq_real_arr2D, assert_eq_cmplx_arr2D, assert_eq_string_arr2D, &
                              assert_eq_int_arr3D, assert_eq_real_arr3D, assert_eq_cmplx_arr3D, assert_eq_string_arr3D

    ! private procedures
    procedure, private :: assert_1
    procedure, private :: assert_arr
    procedure, private :: assert_arr2D
    procedure, private :: assert_arr3D
    procedure, private :: assert_eq_int
    procedure, private :: assert_eq_int_arr
    procedure, private :: assert_eq_int_arr2D
    procedure, private :: assert_eq_int_arr3D
    procedure, private :: assert_eq_real
    procedure, private :: assert_eq_real_arr
    procedure, private :: assert_eq_real_arr2D
    procedure, private :: assert_eq_real_arr3D
    procedure, private :: assert_eq_cmplx
    procedure, private :: assert_eq_cmplx_arr
    procedure, private :: assert_eq_cmplx_arr2D
    procedure, private :: assert_eq_cmplx_arr3D
    procedure, private :: assert_eq_string
    procedure, private :: assert_eq_string_arr
    procedure, private :: assert_eq_string_arr2D
    procedure, private :: assert_eq_string_arr3D
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

    print *
    print "(A)", COL_WHITE//name//": "//COL_DEFAULT//"Begin testing..."
  end subroutine

  subroutine finish(this)
    !! Finish test and print number of passed tests.

    class(test_case), intent(in) :: this
      !! Test case

    character(7) :: color

    if (this%passed_tests == this%num_tests) then
      color = COL_GREEN
    else if (this%passed_tests == 0) then
      color = COL_RED
    else
      color = COL_YELLOW
    end if

    print "(A,I0,A,I0,A)", COL_WHITE//this%name//": "//color, this%passed_tests, "/", this%num_tests, " tests passed!"//COL_DEFAULT

    if (this%passed_tests /= this%num_tests) then
      print *
      call program_error("At least one test failed!")
    end if
  end subroutine

  subroutine assert_1(this, value, msg)
    !! Check logical value for true.

    class(test_case), intent(inout) :: this
      !! Test case
    logical,          intent(in)    :: value
      !! Logical value
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    call this%check(value, msg, msg2 = "assertion for scalar logical failed")
  end subroutine

  subroutine assert_arr(this, values, msg)
    !! Check logical values (1D array) for true.

    class(test_case), intent(inout) :: this
      !! Test case
    logical,          intent(in)    :: values(:)
      !! Logical values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, failed

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check(values(i), msg, msg2 = "assertion for array values failed; (row)", ipar = [i])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_arr2D(this, values, msg)
    !! Check logical values (1D array) for true.

    class(test_case), intent(inout) :: this
      !! Test case
    logical,          intent(in)    :: values(:,:)
      !! Logical values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, failed

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        call this%check(values(i,j), msg, msg2 = "assertion for array values failed; (row, col)", ipar = [i, j])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_arr3D(this, values, msg)
    !! Check logical values (1D array) for true.

    class(test_case), intent(inout) :: this
      !! Test case
    logical,          intent(in)    :: values(:,:,:)
      !! Logical values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, k, failed

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        do k = 1, size(values,3)
          call this%check(values(i,j,k), msg, msg2 = "assertion for array values failed; (row, col, slice)", ipar = [i, j, k])
          if (.not. this%last_passed) failed = failed + 1
          if (failed > 5) then
            print "(A)", "    Failed for more than 5 elements, discontinue test!"
            return
          end if
        end do
      end do
    end do
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
        print "(A)", "    Failed for more than 5 elements, discontinue test!"
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
          print "(A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_int_arr3D(this, expected, values, msg)
    !! Check integer values (3D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    integer,          intent(in)    :: expected(:,:,:)
      !! Expected values
    integer,          intent(in)    :: values(:,:,:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, k, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return
    call this%check((size(expected,3) == size(values,3)), msg, &
                    msg2 = "number of slices differ", ipar = [size(expected,3), size(values,3)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        do k = 1, size(values,3)
          call this%check((values(i,j,k) == expected(i,j,k)), msg, &
                          msg2 = "array values differ; (row, col, slice, expected, value)", ipar = [i, j, k, expected(i,j,k), values(i,j,k)])
          if (.not. this%last_passed) failed = failed + 1
          if (failed > 5) then
            print "(A)", "    Failed for more than 5 elements, discontinue test!"
            return
          end if
        end do
      end do
    end do
  end subroutine

  subroutine assert_eq_real(this, expected, value, rtol, atol, msg)
    !! check real value (scalar) for equality to expected value within tolerance

    class(test_case),  intent(inout) :: this
      !! test case
    real,              intent(in)    :: expected
      !! expected value
    real,              intent(in)    :: value
      !! actual value
    real,              intent(in)    :: rtol
      !! relative tolerance to compare value to
    real,              intent(in)    :: atol
      !! absolute tolerance to compare value to
    character(*),      intent(in)    :: msg
      !! message to print if error detected

    call this%check(abs(expected - value) <= max(atol, rtol * abs(expected)), msg, &
                    msg2 = "scalar value differs", rpar = [expected, value])
  end subroutine

  subroutine assert_eq_real_arr(this, expected, values, rtol, atol, msg)
    !! Check real values (1D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected(:)
      !! Expected values
    real,             intent(in)    :: values(:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, failed

    ! check array sizes
    call this%check((size(expected) == size(values)), msg, &
                    msg2 = "array lengths differ")
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check(abs(expected(i) - values(i)) <= max(atol, rtol * abs(expected(i))), msg, &
                      msg2 = "array values differ; (row; expected, value)", ipar = [i], rpar = [expected(i), values(i)])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_eq_real_arr2D(this, expected, values, rtol, atol, msg)
    !! Check real values (2D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected(:,:)
      !! Expected values
    real,             intent(in)    :: values(:,:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

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
        call this%check(abs(expected(i,j) - values(i,j)) <= max(atol, rtol * abs(expected(i,j))), msg, &
                        msg2 = "array values differ; (row, col; expected, value)", ipar = [i, j], rpar = [expected(i,j), values(i,j)])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_real_arr3D(this, expected, values, rtol, atol, msg)
    !! Check real values (3D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    real,             intent(in)    :: expected(:,:,:)
      !! Expected values
    real,             intent(in)    :: values(:,:,:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, k, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return
    call this%check((size(expected,3) == size(values,3)), msg, &
                    msg2 = "number of slices differ", ipar = [size(expected,3), size(values,3)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        do k = 1, size(values,3)
          call this%check(abs(expected(i,j,k) - values(i,j,k)) <= max(atol, rtol * abs(expected(i,j,k))), msg, &
                          msg2 = "array values differ; (row, col, slice; expected, value)", ipar = [i, j, k], rpar = [expected(i,j,k), values(i,j,k)])
          if (.not. this%last_passed) failed = failed + 1
          if (failed > 5) then
            print "(A)", "    Failed for more than 5 elements, discontinue test!"
            return
          end if
        end do
      end do
    end do
  end subroutine

  subroutine assert_eq_cmplx(this, expected, value, rtol, atol, msg)
    !! Check real value (scalar) for equality to expected value within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    complex,          intent(in)    :: expected
      !! Expected value
    complex,          intent(in)    :: value
      !! Actual value
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    call this%check(abs(expected - value) <= max(atol, rtol * abs(expected)), msg, &
                    msg2 = "scalar value differs", cpar = [expected, value])
  end subroutine

  subroutine assert_eq_cmplx_arr(this, expected, values, rtol, atol, msg)
    !! Check real values (1D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    complex,          intent(in)    :: expected(:)
      !! Expected values
    complex,          intent(in)    :: values(:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, failed

    ! check array sizes
    call this%check((size(expected) == size(values)), msg, &
                    msg2 = "array lengths differ")
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values)
      call this%check(abs(expected(i) - values(i)) <= max(atol, rtol * abs(expected(i))), msg, &
                      msg2 = "array values differ; (row; expected, value)", ipar = [i], cpar = [expected(i), values(i)])
      if (.not. this%last_passed) failed = failed + 1
      if (failed > 5) then
        print "(A)", "    Failed for more than 5 elements, discontinue test!"
        return
      end if
    end do
  end subroutine

  subroutine assert_eq_cmplx_arr2D(this, expected, values, rtol, atol, msg)
    !! Check complex values (2D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    complex,          intent(in)    :: expected(:,:)
      !! Expected values
    complex,          intent(in)    :: values(:,:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

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
        call this%check(abs(expected(i,j) - values(i,j)) <= max(atol, rtol * abs(expected(i,j))), msg, &
                        msg2 = "array values differ; (row, col; expected, value)", ipar = [i, j], cpar = [expected(i,j), values(i,j)])
        if (.not. this%last_passed) failed = failed + 1
        if (failed > 5) then
          print "(A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_cmplx_arr3D(this, expected, values, rtol, atol, msg)
    !! Check complex values (3D array) for equality to expected values within tolerance.

    class(test_case), intent(inout) :: this
      !! Test case
    complex,          intent(in)    :: expected(:,:,:)
      !! Expected values
    complex,          intent(in)    :: values(:,:,:)
      !! Actual values
    real,             intent(in)    :: rtol
      !! relative tolerance
    real,             intent(in)    :: atol
      !! absolute tolerance
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, k, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return
    call this%check((size(expected,3) == size(values,3)), msg, &
                    msg2 = "number of slices differ", ipar = [size(expected,3), size(values,3)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        do k = 1, size(values,3)
          call this%check(abs(expected(i,j,k) - values(i,j,k)) <= max(atol, rtol * abs(expected(i,j,k))), msg, &
                          msg2 = "array values differ; (row, col, slice; expected, value)", ipar = [i, j, k], cpar = [expected(i,j,k), values(i,j,k)])
          if (.not. this%last_passed) failed = failed + 1
          if (failed > 5) then
            print "(A)", "    Failed for more than 5 elements, discontinue test!"
            return
          end if
        end do
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
        print "(A)", "    Failed for more than 5 elements, discontinue test!"
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
          print "(A)", "    Failed for more than 5 elements, discontinue test!"
          return
        end if
      end do
    end do
  end subroutine

  subroutine assert_eq_string_arr3D(this, expected, values, msg)
    !! Check strings (3D array) for equality to expected values.

    class(test_case), intent(inout) :: this
      !! Test case
    type(string),     intent(in)    :: expected(:,:,:)
      !! Expected values
    type(string),     intent(in)    :: values(:,:,:)
      !! Actual values
    character(*),     intent(in)    :: msg
      !! Message to print if error detected

    integer :: i, j, k, failed

    ! check array sizes
    call this%check((size(expected,1) == size(values,1)), msg, &
                    msg2 = "number of rows differ", ipar = [size(expected,1), size(values,1)])
    if (.not. this%last_passed) return
    call this%check((size(expected,2) == size(values,2)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,2), size(values,2)])
    if (.not. this%last_passed) return
    call this%check((size(expected,3) == size(values,3)), msg, &
                    msg2 = "number of cols differ", ipar = [size(expected,3), size(values,3)])
    if (.not. this%last_passed) return

    ! check values one by one
    failed = 0
    do i = 1, size(values,1)
      do j = 1, size(values,2)
        do k = 1, size(values,3)
          call this%check((values(i,j,k)%s == expected(i,j,k)%s), msg, &
                          msg2 = "array values differ; (row, col, slice)", ipar = [i, j, k])
          if (.not. this%last_passed) failed = failed + 1
          if (failed > 5) then
            print "(A)", "    Failed for more than 5 elements, discontinue test!"
            return
          end if
        end do
      end do
    end do
  end subroutine

  subroutine check(this, chk, msg, msg2, ipar, rpar, cpar)
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
    complex,      optional, intent(in)    :: cpar(:)
      !! Additional complex parameters to print (new line)

    integer :: i

    ! count total number of tests
    this%num_tests = this%num_tests + 1

    if (chk) then
      this%last_passed  = .true.
      this%passed_tests = this%passed_tests + 1
    else
      this%last_passed  = .false.
      print "(A)", COL_WHITE//trim(this%name)//": "//COL_MAGENTA//"Failed test"//COL_DEFAULT
      print "(A)", "  "//trim(msg)
      if (present(msg2)) print "(A)", "    "//trim(msg2)
      if (present(ipar)) then
        write(*, fmt="(A)", advance="no") "      "
        do i = 1, size(ipar)
          if (i < size(ipar)) then
            write(*, fmt="(I0,A)", advance="no") ipar(i), ", "
          else
            write(*, fmt="(I0)") ipar(i)
          end if
        end do
        print *
      end if
      if (present(rpar)) then
        write(*, fmt="(1A)", advance="no") "      "
        do i = 1, size(rpar)
          if (i < size(rpar)) then
            write(*, fmt="(ES25.16E3,A)", advance="no") rpar(i), ", "
          else
            write(*, fmt="(ES25.16E3)") rpar(i)
          end if
        end do
        print *
      end if
      if (present(cpar)) then
        write(*, fmt="(1A)", advance="no") "      "
        do i = 1, size(cpar)
          if (i < size(cpar)) then
            write(*, fmt="(2ES25.16E3,A)", advance="no") cpar(i), ", "
          else
            write(*, fmt="(2ES25.16E3)") cpar(i)
          end if
        end do
        print *
      end if
    end if
  end subroutine

end module
