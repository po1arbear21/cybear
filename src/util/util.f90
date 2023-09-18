m4_include(macro.f90.inc)

module util_m

  use error_m,      only: assert_failed, program_error
  use iso_c_binding
  use string_m,     only: string
  use vector_m,     only: vector_int, vector_real

  implicit none

  private
  public fsleep
  public get_hostname
  public strlen, cstrlen, c2fstring, f2cstring
  public hash
  public int2str, log2str, real2str
  public is_digit, is_letter, is_whitespace
  public split_string, split_folder_file
  public load_array
  public select_int
  public get_memory_usage

  interface hash
    module procedure :: hash_int32, hash_int32_array, hash_int64, hash_int64_array
  end interface

  ! interfaces to C routines
  interface
    function c_hash_int32(i) result(h) bind(c)
      import :: c_int32_t
      integer(c_int32_t), value, intent(in) :: i
      integer(c_int32_t)                    :: h
    end function
    function c_hash_int32_array(a, n) result(h) bind(c)
      import :: c_int32_t, c_size_t
      integer(c_int32_t),       intent(in) :: a(*)
      integer(c_size_t), value, intent(in) :: n
      integer(c_int32_t)                   :: h
    end function
    function c_hash_int64(i) result(h) bind(c)
      import :: c_int64_t
      integer(c_int64_t), value, intent(in) :: i
      integer(c_int64_t)                    :: h
    end function
    function c_hash_int64_array(a, n) result(h) bind(c)
      import :: c_int64_t, c_size_t
      integer(c_int64_t),       intent(in) :: a(*)
      integer(c_size_t), value, intent(in) :: n
      integer(c_int64_t)                   :: h
    end function
  end interface

  interface
     function strlen(s) result(len) bind(c)
       import :: c_size_t, c_ptr
       type(c_ptr), value, intent(in) :: s
       integer(c_size_t)              :: len
     end function strlen
  end interface

contains

  subroutine fsleep(seconds)
    !! sleep for a number of seconds
    m4_ifdef({m4_intel},use ifport)
    integer, intent(in) :: seconds

    m4_ifdef({m4_intel},{
    m4_ifelse(m4_intsize,32,{
    call sleep(seconds)
    },{
    integer(kind=4) :: s
    s = int(seconds, kind = 4)
    call sleep(s)
    })
    })

    m4_ifdef({m4_gnu},{
    call sleep(seconds)
    })
  end subroutine

  function get_hostname() result(hostname)
    !! get host name
    m4_ifdef({m4_intel},use ifport)
    character(:), allocatable :: hostname

    character(256) :: hn

    m4_ifdef({m4_intel},{
    integer(4) :: istat

    istat = hostnam(hn)
    })

    m4_ifdef({m4_gnu},{
    call hostnm(hn)
    })

    hostname = trim(hn)
  end function

  pure function cstrlen(cstr) result(len)
    !! get length of c string
    character(1), intent(in) :: cstr(*)
    integer                  :: len

    ! go through cstr until null character is found
    len = 1
    do while (cstr(len) /= c_null_char)
      len = len + 1
    end do

    ! discard trailing null character
    len = len - 1
  end function

  pure function c2fstring(cstr) result(fstr)
    !! convert c string to fortran string
    character(1), intent(in) :: cstr(*)
    character(cstrlen(cstr)) :: fstr

    fstr = transfer(cstr(1:len(fstr)), fstr)
  end function

  pure function f2cstring(fstr) result(cstr)
    !! convert fortran string to c string
    character(*), intent(in) :: fstr
    character(1)             :: cstr(len_trim(fstr) + 1)

    ! local variables
    integer :: i

    do i = 1, size(cstr) - 1
      cstr(i) = fstr(i:i)
    end do
    cstr(size(cstr)) = c_null_char
  end function

  function int2str(i, fmt) result(str)
    !! convert integer to string
    integer,                intent(in)  :: i
    character(*), optional, intent(in)  :: fmt
    character(:),           allocatable :: str

    character(256) :: tmp

    if (present(fmt)) then
      write(tmp, fmt) i
    else
      write(tmp, "(I24)") i
    end if
    str = trim(adjustl(tmp))
  end function

  function log2str(l) result(str)
    !! convert real to string
    logical,     intent(in)  :: l
    character(1)             :: str

    write(str, "(L)") l
  end function

  function real2str(r, fmt) result(str)
    !! convert real to string
    real,                   intent(in)  :: r
    character(*), optional, intent(in)  :: fmt
    character(:),           allocatable :: str

    character(256) :: tmp

    if (present(fmt)) then
      write(tmp, fmt) r
    else
      write(tmp, "(ES25.16E3)") r
    end if

    str = trim(adjustl(tmp))
  end function

  pure function is_digit(ch) result(ret)
    character(1), intent(in) :: ch
    logical                  :: ret

    ret = ((ch >= "0") .and. (ch <= "9"))
  end function

  pure function is_letter(ch) result(ret)
    character(1), intent(in) :: ch
    logical                  :: ret

    ! extended ASCII or UTF8 characters are interpreted as letters
    ret = (((ch >= "a") .and. (ch <= "z")) .or. ((ch >= "A") .and. (ch <= "Z")) .or. (iachar(ch) > 127))
  end function

  pure function is_whitespace(ch) result(ret)
    character(1), intent(in) :: ch
    logical                  :: ret

    integer :: ascii

    ascii = iachar(ch)
    ret   = ((ascii == 32) .or. (ascii == 9) .or. (ascii == 13) .or. (ascii == 10))
  end function

  function split_string(str, delim) result(res)
    !! split string using delimiter
    character(*), intent(in)  :: str
      !! string to split
    character(1), intent(in)  :: delim(:)
      !! list of delimiter characters (are removed from result)
    type(string), allocatable :: res(:)
      !! return substrings

    integer          :: i, j, nlen, ndelim, nres
    logical          :: token, split
    type(vector_int) :: i0, i1

    nlen   = len_trim(str)
    ndelim = size(delim)

    ! get indices of split strings
    call i0%init(0, c = nlen)
    call i1%init(0, c = nlen)
    token = .false.
    nres = 0
    do i = 1, nlen
      do j = 1, ndelim
        if (str(i:i) == delim(j)) exit
      end do
      split = (j <= ndelim)

      if (token .and. split) then
        call i1%push(i - 1)
        token = .false.
      elseif (.not. token .and. .not. split) then
        call i0%push(i)
        token = .true.
      end if
    end do
    if (token) call i1%push(nlen)

    ! copy to result
    nres = i0%n
    allocate (res(nres))
    do i = 1, nres
      res(i)%s = str(i0%d(i):i1%d(i))
    end do
  end function

  subroutine split_folder_file(full, folder, file)
    !! split full file name into folder + file
    character(*),              intent(in)  :: full
      !! full filename with folder
    character(:), allocatable, intent(out) :: folder
      !! output folder (including trailing "/")
    character(:), allocatable, intent(out) :: file
      !! output file name without folder

    integer :: i

    i = scan(full, "/", back = .true.)

    if (i == 0) then
      folder = ""
      file = full
    else
      folder = full(1:i)
      file   = full(i+1:len_trim(full))
    end if
  end subroutine

  function select_int(flags, ints) result(t)
    !! select an integer according to the flag that is set (no flags set => 0, multiple set => - 1)
    logical, intent(in) :: flags(:)
    integer, intent(in) :: ints(:)
    integer             :: t

    ! local variables
    integer :: i, j

    t = 0
    j = 0
    do i = 1, size(flags)
      if (flags(i)) then
        t = ints(i)
        j = j + 1
      end if
    end do

    ! ambiguous
    if (j > 1) t = -1
  end function

  function get_memory_usage() result(rss)
    !! get memory usage of program (RSS)
    real :: rss
      !! return used memory in GiB

    character(80) :: line
    integer       :: funit, ios, rss_kiB

    rss = -1

    open (newunit = funit, file = "/proc/self/status", status = "old", action = "read", iostat = ios)
    if (ios /= 0) return

    read(funit, "(A)", iostat = ios) line
    do while (ios == 0)
      if (line(1:6) == "VmRSS:") then
        read (line(7:), *) rss_kiB
        rss = real(rss_kiB) / 2**20
        exit
      end if
      read(funit, "(A)", iostat = ios) line
    end do

    close (funit)
  end function

  function hash_int32(i) result(h)
    !! 32-bit integer hash function (nullprogram.com/blog/2018/07/31/)
    integer(kind=4), intent(in) :: i
    integer(kind=4)             :: h

    h = c_hash_int32(i)
  end function

  function hash_int32_array(i) result(h)
    !! 32-bit integer hash function for arrays
    integer(kind=4), intent(in) :: i(:)
    integer(kind=4)             :: h

    integer(kind=c_size_t) :: n

    n = size(i)

    h = c_hash_int32_array(i, n)
  end function

  function hash_int64(i) result(h)
    !! 64-bit integer hash function (splitmix64)
    integer(kind=8), intent(in) :: i
    integer(kind=8)             :: h

    h = c_hash_int64(i)
  end function

  function hash_int64_array(i) result(h)
    !! 64-bit integer hash function (splitmix64)
    integer(kind=8), intent(in) :: i(:)
    integer(kind=8)             :: h

    integer(kind=c_size_t) :: n

    n = size(i)

    h = c_hash_int64_array(i, n)
  end function

  subroutine load_array(file, x)
    !! load 1D real array from file
    character(*),      intent(in)  :: file
      !! file name
    real, allocatable, intent(out) :: x(:)
      !! output array read from file

    ! local variables
    integer           :: funit, iostat
    real              :: r
    type(vector_real) :: vec

    call vec%init(0, c = 16)

    ! read file into vector
    open (newunit = funit, file = file, status = "old", action = "read", iostat = iostat)
    if (iostat /= 0) then
      close (funit)
      call program_error("Could not open file "//file)
    end if
    do while (.true.)
      read (unit = funit, fmt = *, iostat = iostat) r
      if (iostat > 0) then
        close(funit)
        call program_error("IO-Error")
      else if (iostat == 0) then
        call vec%push(r)
      else
        exit ! end of file
      end if
    end do
    close(funit)

    ! output array
    allocate (x(vec%n), source = vec%d(1:vec%n))
  end subroutine

end module
