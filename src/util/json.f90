m4_include(macro.f90.inc)

module json_m

  use error_m,         only: assert_failed, program_error
  use iso_fortran_env, only: int32, int64
  use map_m,           only: map_string_int, mapnode_string_int
  use string_m,        only: string, new_string
  use util_m,          only: int2str
  use vector_m,        only: vector_int, vector_string

  implicit none

  private
  public :: JS_NULL, JS_LOG, JS_INT, JS_REAL, JS_STRING, JS_ARRAY, JS_OBJECT
  public :: json_file, json, json_ptr, json_null, json_log, json_int, json_real, json_string, json_array, json_object

  integer, parameter :: JS_NULL   = 0
  integer, parameter :: JS_LOG    = 1
  integer, parameter :: JS_INT    = 2
  integer, parameter :: JS_REAL   = 3
  integer, parameter :: JS_STRING = 4
  integer, parameter :: JS_ARRAY  = 5
  integer, parameter :: JS_OBJECT = 6

  type json_file
    !! json file object
    class(json), pointer :: js => null()
      !! main object
  contains
    procedure :: init     => json_file_init
    procedure :: destruct => json_file_destruct
    procedure :: load     => json_file_load
    procedure :: save     => json_file_save
  end type

  type, abstract :: json
    !! base type of json data structure
  contains
    procedure :: init     => json_init
    procedure :: destruct => json_destruct
    procedure :: get_ptr  => json_get_ptr

    procedure :: cast_null   => json_cast_null
    procedure :: cast_log    => json_cast_log
    procedure :: cast_int    => json_cast_int
    procedure :: cast_real   => json_cast_real
    procedure :: cast_string => json_cast_string
    procedure :: cast_array  => json_cast_array
    procedure :: cast_object => json_cast_object
  end type

  type json_ptr
    class(json), pointer :: p => null()
  end type

  m4_define({T},{json_ptr})
  m4_include(vector_def.f90.inc)

  type, extends(json) :: json_null
    !! null json value
  end type

  type, extends(json) :: json_log
    !! logical json value
    logical :: value
  end type

  type, extends(json) :: json_int
    !! integer json number
    integer(kind=int64) :: value
  end type

  type, extends(json) :: json_real
    !! real json number
    real :: value
  end type

  type, extends(json) :: json_string
    !! string json value
    character(:), allocatable :: value
  end type

  type, extends(json) :: json_array
    !! json array
    type(vector_json_ptr) :: values
  contains
    procedure :: get_log    => json_array_get_log
    procedure :: get_int    => json_array_get_int
    procedure :: get_real   => json_array_get_real
    procedure :: get_string => json_array_get_string
    procedure :: get_array  => json_array_get_array
    procedure :: get_object => json_array_get_object

    procedure :: add_null   => json_array_add_null
    procedure :: add_log    => json_array_add_log
    procedure :: json_array_add_int32
    procedure :: json_array_add_int64
    generic   :: add_int    => json_array_add_int32, json_array_add_int64
    procedure :: add_real   => json_array_add_real
    procedure :: add_string => json_array_add_string
    procedure :: add_array  => json_array_add_array
    procedure :: add_object => json_array_add_object
  end type

  type, extends(json) :: json_object
    !! json object

    type(vector_string)   :: names
      !! names of properties
    type(vector_json_ptr) :: properties
      !! properties
    type(map_string_int)  :: property_map
      !! property name -> property index
  contains
    procedure :: get_json   => json_object_get_json
    procedure :: get_log    => json_object_get_log
    procedure :: get_int    => json_object_get_int
    procedure :: get_real   => json_object_get_real
    procedure :: get_string => json_object_get_string
    procedure :: get_array  => json_object_get_array
    procedure :: get_object => json_object_get_object

    procedure :: set_null   => json_object_set_null
    procedure :: set_log    => json_object_set_log
    procedure :: json_object_set_int32
    procedure :: json_object_set_int64
    generic   :: set_int    => json_object_set_int32, json_object_set_int64
    procedure :: set_real   => json_object_set_real
    procedure :: set_string => json_object_set_string
    procedure :: set_array  => json_object_set_array
    procedure :: set_object => json_object_set_object

    procedure :: add_null   => json_object_add_null
    procedure :: add_log    => json_object_add_log
    procedure :: json_object_add_int32
    procedure :: json_object_add_int64
    generic   :: add_int    => json_object_add_int32, json_object_add_int64
    procedure :: add_real   => json_object_add_real
    procedure :: add_string => json_object_add_string
    procedure :: add_array  => json_object_add_array
    procedure :: add_object => json_object_add_object

    procedure, private :: add_property => json_object_add_property
  end type

contains

  m4_define({T},{json_ptr})
  m4_include(vector_imp.f90.inc)
  m4_changequote({{,}})

  subroutine json_file_init(this, main_type)
    !! initialize json file
    class(json_file),  intent(out) :: this
    integer, optional, intent(in)  :: main_type
      !! type of main object (default: JS_OBJECT)

    integer :: main_type_

    ! main object type
    main_type_ = JS_OBJECT
    if (present(main_type)) main_type_ = main_type

    ! allocate and initialize main object
    select case (main_type_)
    case (JS_NULL)
      allocate (json_null :: this%js)
    case (JS_LOG)
      allocate (json_log :: this%js)
    case (JS_INT)
      allocate (json_int :: this%js)
    case (JS_REAL)
      allocate (json_real :: this%js)
    case (JS_STRING)
      allocate (json_string :: this%js)
    case (JS_ARRAY)
      allocate (json_array :: this%js)
    case (JS_OBJECT)
      allocate (json_object :: this%js)
    end select
    call this%js%init()
  end subroutine

  subroutine json_file_destruct(this)
    !! destruct json file
    class(json_file), intent(inout) :: this

    if (associated(this%js)) then
      call this%js%destruct()
      deallocate (this%js)
    end if
  end subroutine

  subroutine json_file_load(this, fname)
    !! load json file
    class(json_file),  intent(out) :: this
    character(*),      intent(in)  :: fname
      !! file name

    character(80)             :: iomsg
    character(:), allocatable :: buf
    integer                   :: funit, iostat, n

    ! try to open file
    open (newunit = funit, file = fname, access = "stream", form = "unformatted", status = "old", action = "read", iostat = iostat, iomsg = iomsg)
    if (iostat /= 0) call program_error("Error while opening file: " // iomsg)

    ! read whole file into buffer
    inquire (file = fname, size = n)
    allocate (character(n) :: buf)
    read (funit, iostat = iostat, iomsg = iomsg) buf
    close (funit)
    if (iostat /= 0) call program_error("Error while reading file: " // iomsg)

    ! parse buffer
    call read_json(this%js, buf)
  end subroutine

  subroutine json_file_save(this, fname)
    !! save json file
    class(json_file), intent(in) :: this
    character(*),     intent(in) :: fname
      !! file name

    character(80)             :: iomsg
    character(:), allocatable :: buf
    integer                   :: funit, iostat

    ! write json into buffer
    call write_json(this%js, buf)

    ! try to open file
    open (newunit = funit, file = fname, access = "stream", form = "unformatted", status = "replace", action = "write", iostat = iostat, iomsg = iomsg)
    if (iostat /= 0) call program_error("Error while opening file: " // iomsg)

    ! write to file
    write (funit, iostat = iostat, iomsg = iomsg) buf
    if (iostat /= 0) call program_error("Error while writing to file: " // iomsg)

    ! newline at the end of file
    write (funit) char(10)

    ! close file
    close(funit)
  end subroutine

  subroutine read_json(js, str)
    !! parse string into json data structure
    class(json), pointer, intent(out) :: js
      !! output json structure
    character(*),         intent(in)  :: str
      !! string to parse

    character(1)         :: c, scope
    integer              :: istr, jstr, nstr, line
    integer, allocatable :: isp(:), line_nr(:)
    type(vector_int)     :: stack

    ! length of string
    nstr = len_trim(str)

    ! allocate scope pointers and stack
    allocate (isp(nstr), line_nr(nstr), source = 0)
    call stack%init(0, c = 16)

    ! set scope pointers
    line  = 1
    istr  = 0
    jstr  = 0
    scope = ' '
    do while (istr < nstr)
      ! get next character
      istr = istr + 1
      c = str(istr:istr)

      ! save line number
      line_nr(istr) = line

      select case (scope)
      case (' ')
        select case (c)
        case (' ', char(10), char(13))
        case ('"', '[', '{', 'n', 't', 'f', '-', '0':'9')
          call push()
        case default
          call program_error("Error in line " // int2str(line) // ": unexpected character '" // c // "'")
        end select
      case ('{')
        select case (c)
        case (' ', char(10), char(13))
        case (',', ':')
          isp(istr) = jstr
        case ('"', '[', '{', 'n', 't', 'f', '-', '0':'9')
          call push()
        case ('}')
          call pop()
        case default
          call program_error("Error in line " // int2str(line) // ": unexpected character '" // c // "'")
        end select
      case ('[')
        select case (c)
        case (' ', char(10), char(13))
        case (',')
          isp(istr) = jstr
        case ('"', '[', '{', 'n', 't', 'f', '-', '0':'9')
          call push()
        case (']')
          call pop()
        case default
          call program_error("Error in line " // int2str(line) // ": unexpected character '" // c // "'")
        end select
      case ('"')
        select case (c)
        case (char(10), char(13))
          call program_error("Error in line " // int2str(line) // ": unterminated string")
        case ('"')
          call pop()
        case ('\')
          call push()
        end select
      case ('\')
        select case (c)
        case (char(10), char(13))
          call program_error("Error in line " // int2str(line) // ": unexpected end of line")
        case ('"', '\', '/', 'b', 'f', 'n', 'r', 't')
          call pop()
        case default
          call program_error("Error in line " // int2str(line) // ": unknown escape sequence '\" // c // "'")
        end select
      case ('n', 't', 'f', '-', '0':'9')
        select case (c)
        case (char(10), char(13), ',', ']', '}')
          istr = istr - 1
          call pop()
        end select
      end select

      if ((c == char(10)) .or. (c == char(13))) line = line + 1
    end do
    if (stack%n /= 0) call program_error("Error in line " // int2str(line) // ": unterminated '" // scope // "'")

    ! parse recursively
    call parse(js, 1, nstr)

  contains

    subroutine push()
      call stack%push(istr)
      jstr  = istr
      scope = c
    end subroutine

    subroutine pop()
      isp(jstr) = istr
      isp(istr) = jstr
      call stack%pop()
      if (stack%n > 0) then
        jstr  = stack%d(stack%n)
        scope = str(jstr:jstr)
      else
        jstr  = 0
        scope = ' '
      end if
    end subroutine

    recursive subroutine parse(p, i0, i1)
      class(json), pointer, intent(out) :: p
      integer,              intent(in)  :: i0
      integer,              intent(in)  :: i1

      integer                   :: i, j, k, l, iostat
      integer(kind=int64)       :: i64
      logical                   :: first, lname
      character(:), allocatable :: name
      real                      :: r
      type(json_ptr)            :: ptr

      line = line_nr(i0)

      select case (str(i0:i0))
      case ('n')
        if (        i1 /= i0 + 3) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        if (str(i0:i1) /= "null") call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        allocate (json_null :: p)

      case ('t')
        if (        i1 /= i0 + 3) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        if (str(i0:i1) /= "true") call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        allocate (json_log :: p)
        select type (p)
        type is (json_log)
          p%value = .true.
        end select

      case ('f')
        if (        i1 /=  i0 + 4) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        if (str(i0:i1) /= "false") call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        allocate (json_log :: p)
        select type (p)
        type is (json_log)
          p%value = .false.
        end select

      case ('-', '0':'9')
        l = scan(str(i0:i1), '.') + scan(str(i0:i1), 'd') + scan(str(i0:i1), 'e')

        if (l == 0) then
          read (str(i0:i1), *, iostat = iostat) i64
          if (iostat /= 0) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
          allocate (json_int :: p)
          select type (p)
          type is (json_int)
            p%value = i64
          end select
        else
          read (str(i0:i1), *, iostat = iostat) r
          if (iostat /= 0) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
          allocate (json_real :: p)
          select type (p)
          type is (json_real)
            p%value = r
          end select
        end if

      case ('"')
        allocate (json_string :: p)
        select type (p)
        type is (json_string)
          call unescape_str(str(i0+1:i1-1), p%value)
        end select

      case ('[')
        allocate (json_array :: p)
        call p%init()
        select type (p)
        type is (json_array)
          j = i0
          first = .true.
          do while (j < i1)
            i = j + 1
            do j = i + 1, i1
              if (isp(j) == i0) exit
            end do
            do while (((str(i:i) == ' ') .or. (str(i:i) == char(10)) .or. (str(i:i) == char(13))) .and. (i <= j))
              i = i + 1
            end do
            if (i == j) then
              if (first .and. (j == i1)) exit
              call program_error("Error in line " // int2str(line) // ": array element expected")
            end if
            k = j - 1
            do while (((str(k:k) == ' ') .or. (str(k:k) == char(10)) .or. (str(k:k) == char(13))) .and. (k > i))
              k = k - 1
            end do
            call parse(ptr%p, i, k)
            call p%values%push(ptr)
            first = .false.
          end do
        end select

      case ('{')
        allocate (json_object :: p)
        call p%init()
        select type (p)
        type is (json_object)
          j = i0
          first = .true.
          lname = .true.
          do while (j < i1)
            i = j + 1
            if (i == i1) exit
            do j = i + 1, i1
              if (isp(j) == i0) exit
            end do
            if (lname .and. (str(j:j) /= ':')) call program_error("Error in line " // int2str(line) // ": missing ':' between object property name and value")
            if (.not. lname .and. ((str(j:j) /= ',') .and. (str(j:j) /= '}'))) call program_error("Error in line " // int2str(line) // ": missing ',' between object properties")
            do while (((str(i:i) == ' ') .or. (str(i:i) == char(10)) .or. (str(i:i) == char(13))) .and. (i <= j))
              i = i + 1
            end do
            if (i == j) then
              if (first .and. lname .and. (j == i1)) exit
              if (lname) then
                call program_error("Error in line " // int2str(line) // ": missing object property name")
              else
                call program_error("Error in line " // int2str(line) // ": missing object property value")
              end if
            end if
            k = j - 1
            do while (((str(k:k) == ' ') .or. (str(k:k) == char(10)) .or. (str(k:k) == char(13))) .and. (k > i))
              k = k - 1
            end do
            call parse(ptr%p, i, k)
            if (lname) then
              select type (s => ptr%p)
              type is (json_string)
                name = s%value
              class default
                call program_error("Error in line " // int2str(line) // ": object property name must be a string")
              end select
              deallocate(ptr%p)
              lname = .false.
            else
              call p%add_property(name, ptr)
              lname = .true.
            end if
          end do
        end select
      end select
    end subroutine

    subroutine unescape_str(s1, s2)
      character(*),              intent(in)  :: s1
      character(:), allocatable, intent(out) :: s2

      character(:), allocatable :: buf
      integer                   :: i, n, ibuf
      logical                   :: escape

      n = len_trim(s1)
      allocate (character(n) :: buf)
      ibuf = 0
      escape = .false.
      do i = 1, n
        c = s1(i:i)
        if (escape) then
          escape = .false.
          ibuf = ibuf + 1
          select case (c)
          case ('"', '\', '/')
            buf(ibuf:ibuf) = c
          case ('b')
            buf(ibuf:ibuf) = char( 8) ! backspace
          case ('f')
            buf(ibuf:ibuf) = char(12) ! form feed
          case ('n')
            buf(ibuf:ibuf) = char(10) ! new line
          case ('r')
            buf(ibuf:ibuf) = char(13) ! carriage return
          case ('t')
            buf(ibuf:ibuf) = char( 9) ! tab
          case default
            call program_error("Error in line " // int2str(line) // ": unknown escape sequence: '\" // c // "'")
          end select
        else
          if (c == '\') then
            escape = .true.
          else
            ibuf = ibuf + 1
            buf(ibuf:ibuf) = c
          end if
        end if
      end do
      s2 = buf(1:ibuf)
    end subroutine

  end subroutine

  recursive subroutine write_json(js, str, indent)
    !! write json data structure into string
    class(json),               intent(in)  :: js
      !! json data structure
    character(:), allocatable, intent(out) :: str
      !! output json data string
    character(*), optional,    intent(in)  :: indent
      !! current indentation

    character(*), parameter :: INDENT_INC = "  "
    character(*), parameter :: REAL_FMT   = "(ES23.16)"
    integer,      parameter :: REAL_LEN   = 24

    character(:), allocatable :: tmp, indent_
    integer                   :: i

    indent_ = ""
    if (present(indent)) indent_ = indent

    select type (js)
    type is (json_null)
      str = "null"
    type is (json_log)
      if (js%value) then
        str = "true"
      else
        str = "false"
      end if
    type is (json_int)
      allocate (character(32) :: tmp)
      write (tmp, "(I0)") js%value
      str = trim(tmp)
    type is (json_real)
      allocate (character(REAL_LEN) :: tmp)
      write (tmp, REAL_FMT) js%value
      str = trim(tmp)
    type is (json_string)
      call escape_string(js%value, tmp)
      str = tmp
    type is (json_array)
      if (js%values%n == 0) then
        str = "[]"
      else
        str = "[" // char(10)
        do i = 1, js%values%n
          call write_json(js%values%d(i)%p, tmp, indent_ // INDENT_INC)
          str = str // indent_ // INDENT_INC // tmp
          if (i < js%values%n) str = str // ","
          str = str // char(10)
        end do
        str = str // indent_ // "]"
      end if
    type is (json_object)
      if (js%properties%n == 0) then
        str = "{}"
      else
        str = "{" // char(10)
        do i = 1, js%properties%n
          call escape_string(js%names%d(i)%s, tmp)
          str = str // indent_ // INDENT_INC // tmp // ": "
          call write_json(js%properties%d(i)%p, tmp, indent_ // INDENT_INC)
          str = str // tmp
          if (i < js%properties%n) str = str // ","
          str = str // char(10)
        end do
        str = str // indent_ // "}"
      end if
    end select

  contains

    subroutine escape_string(s1, s2)
      character(*),              intent(in)  :: s1
      character(:), allocatable, intent(out) :: s2

      character(1) :: c
      integer      :: n, j, k

      n = len_trim(s1)
      allocate (character(2+2*n) :: s2)
      j = 1
      s2(j:j) = '"'
      do k = 1, n
        c = s1(k:k)
        select case (c)
        case (char( 8)) ! backspace
          j = j + 2
          s2(j-1:j) = "\b"
        case (char(12)) ! form feed
          j = j + 2
          s2(j-1:j) = "\f"
        case (char(10)) ! new line
          j = j + 2
          s2(j-1:j) = "\n"
        case (char(13)) ! carriage return
          j = j + 2
          s2(j-1:j) = "\r"
        case (char( 9)) ! tab
          j = j + 2
          s2(j-1:j) = "\t"
        case ('"')
          j = j + 2
          s2(j-1:j) = '\"'
        case ('\')
          j = j + 2
          s2(j-1:j) = '\\'
        case ('/')
          j = j + 2
          s2(j-1:j) = '\/'
        case default
          j = j + 1
          s2(j:j) = c
        end select
      end do
      j = j + 1
      s2(j:j) = '"'
      s2 = s2(1:j)
    end subroutine

  end subroutine

  subroutine json_init(this)
    !! initialize json data structure
    class(json), intent(out) :: this

    select type (this)
    type is (json_array)
      call this%values%init(0, c = 8)

    type is (json_object)
      call this%names%init(0, c = 8)
      call this%properties%init(0, c = 8)
      call this%property_map%init()

    end select
  end subroutine

  subroutine json_destruct(this)
    !! destruct json data structure
    class(json), intent(inout) :: this

    integer :: i

    select type (this)
    type is (json_array)
      if (allocated(this%values%d)) then
        do i = 1, this%values%n
          if (associated(this%values%d(i)%p)) then
            call this%values%d(i)%p%destruct()
            deallocate (this%values%d(i)%p)
          end if
        end do
        call this%values%destruct()
      end if

    type is (json_object)
      call this%names%destruct()
      if (allocated(this%properties%d)) then
        do i = 1, this%properties%n
          if (associated(this%properties%d(i)%p)) then
            call this%properties%d(i)%p%destruct()
            deallocate (this%properties%d(i)%p)
          end if
        end do
        call this%properties%destruct()
      end if
      call this%property_map%destruct()

    end select
  end subroutine

  function json_get_ptr(this) result(ptr)
    !! get pointer to this json data structure
    class(json), target, intent(in) :: this
    type(json_ptr)                  :: ptr

    ptr%p => this
  end function

  function json_cast_null(this) result(p)
    !! cast json data structure to json_null
    class(json), target, intent(in) :: this
    type(json_null), pointer        :: p
      !! return pointer to json_null (or null if different type)

    select type (this)
    type is (json_null)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_log(this) result(p)
    !! cast json data structure to json_log
    class(json), target, intent(in) :: this
    type(json_log), pointer         :: p
      !! return pointer to json_log (or null if different type)

    select type (this)
    type is (json_log)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_int(this) result(p)
    !! cast json data structure to json_int
    class(json), target, intent(in) :: this
    type(json_int), pointer         :: p
      !! return pointer to json_int (or null if different type)

    select type (this)
    type is (json_int)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_real(this) result(p)
    !! cast json data structure to json_real
    class(json), target, intent(in) :: this
    type(json_real), pointer        :: p
      !! return pointer to json_real (or null if different type)

    select type (this)
    type is (json_real)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_string(this) result(p)
    !! cast json data structure to json_string
    class(json), target, intent(in) :: this
    type(json_string), pointer      :: p
      !! return pointer to json_string (or null if different type)

    select type (this)
    type is (json_string)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_array(this) result(p)
    !! cast json data structure to json_array
    class(json), target, intent(in) :: this
    type(json_array), pointer       :: p
      !! return pointer to json_array (or null if different type)

    select type (this)
    type is (json_array)
      p => this
    class default
      p => null()
    end select
  end function

  function json_cast_object(this) result(p)
    !! cast json data structure to json_object
    class(json), target, intent(in) :: this
    type(json_object), pointer      :: p
      !! return pointer to json_object (or null if different type)

    select type (this)
    type is (json_object)
      p => this
    class default
      p => null()
    end select
  end function

  function json_array_get_log(this, i) result(value)
    !! get logical value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    logical                       :: value

    type(json_log), pointer :: js

    js => this%values%d(i)%p%cast_log()
    m4_assert(associated(js))
    value = js%value
  end function

  function json_array_get_int(this, i) result(value)
    !! get integer value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    integer(int64)                :: value

    type(json_int), pointer :: js

    js => this%values%d(i)%p%cast_int()
    m4_assert(associated(js))
    value = js%value
  end function

  function json_array_get_real(this, i) result(value)
    !! get real value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    real                          :: value

    type(json_real), pointer :: js

    js => this%values%d(i)%p%cast_real()
    m4_assert(associated(js))
    value = js%value
  end function

  function json_array_get_string(this, i) result(value)
    !! get string value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    character(:), allocatable     :: value

    type(json_string), pointer :: js

    js => this%values%d(i)%p%cast_string()
    m4_assert(associated(js))
    value = js%value
  end function

  function json_array_get_array(this, i) result(value)
    !! get array value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    type(json_array), pointer     :: value

    value => this%values%d(i)%p%cast_array()
  end function

  function json_array_get_object(this, i) result(value)
    !! get object value
    class(json_array), intent(in) :: this
    integer,           intent(in) :: i
      !! index
    type(json_object), pointer    :: value

    value => this%values%d(i)%p%cast_object()
  end function

  subroutine json_array_add_null(this, p)
    !! add json_null to array
    class(json_array),                  intent(inout) :: this
    type(json_null), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_null object

    class(json),     pointer :: js
    type(json_null), pointer :: p_

    allocate (json_null :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_null()
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_log(this, value, p)
    !! add json_log to array
    class(json_array),                 intent(inout) :: this
    logical,                           intent(in)    :: value
    type(json_log), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_log object

    class(json),    pointer :: js
    type(json_log), pointer :: p_

    allocate (json_log :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_log()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_int32(this, value, p)
    !! add json_int to array
    class(json_array),                 intent(inout) :: this
    integer(kind=int32),               intent(in)    :: value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    class(json),    pointer :: js
    type(json_int), pointer :: p_

    allocate (json_int :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_int()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_int64(this, value, p)
    !! add json_int to array
    class(json_array),                 intent(inout) :: this
    integer(kind=int64),               intent(in)    :: value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    class(json),    pointer :: js
    type(json_int), pointer :: p_

    allocate (json_int :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_int()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_real(this, value, p)
    !! add json_real to array
    class(json_array),                  intent(inout) :: this
    real,                               intent(in)    :: value
    type(json_real), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_real object

    class(json),     pointer :: js
    type(json_real), pointer :: p_

    allocate (json_real :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_real()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_string(this, value, p)
    !! add json_string to array
    class(json_array),                    intent(inout) :: this
    character(*),                         intent(in)    :: value
    type(json_string), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_string object

    class(json),       pointer :: js
    type(json_string), pointer :: p_

    allocate (json_string :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_string()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_array(this, p)
    !! add json_array to array
    class(json_array),                   intent(inout) :: this
    type(json_array), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_array object

    class(json),      pointer :: js
    type(json_array), pointer :: p_

    allocate (json_array :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_array()
    call p_%init()
    if (present(p)) p => p_
  end subroutine

  subroutine json_array_add_object(this, p)
    !! add json_object to array
    class(json_array),                    intent(inout) :: this
    type(json_object), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_object object

    class(json),       pointer :: js
    type(json_object), pointer :: p_

    allocate (json_object :: js)
    call this%values%push(js%get_ptr())
    p_ => js%cast_object()
    call p_%init()
    if (present(p)) p => p_
  end subroutine

  function json_object_get_json(this, name) result(p)
    !! get json data structure by property name
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    class(json), pointer           :: p
      !! return pointer to property (or null if property name not found)

    type(mapnode_string_int), pointer :: nd

    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then
      p => this%properties%d(nd%value)%p
    else
      p => null()
    end if
  end function

  function json_object_get_log(this, name) result(value)
    !! get logical value
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    logical                        :: value

    class(json),    pointer :: js
    type(json_log), pointer :: js2

    js => this%get_json(name)
    m4_assert(associated(js))
    js2 => js%cast_log()
    m4_assert(associated(js2))
    value = js2%value
  end function

  function json_object_get_int(this, name) result(value)
    !! get integer value
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    integer(kind=int64)            :: value

    class(json),    pointer :: js
    type(json_int), pointer :: js2

    js => this%get_json(name)
    m4_assert(associated(js))
    js2 => js%cast_int()
    m4_assert(associated(js2))
    value = js2%value
  end function

  function json_object_get_real(this, name) result(value)
    !! get real value
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    real                           :: value

    class(json),     pointer :: js
    type(json_real), pointer :: js2

    js => this%get_json(name)
    m4_assert(associated(js))
    js2 => js%cast_real()
    m4_assert(associated(js2))
    value = js2%value
  end function

  function json_object_get_string(this, name) result(value)
    !! get string value
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    character(:), allocatable      :: value

    class(json),       pointer :: js
    type(json_string), pointer :: js2

    js => this%get_json(name)
    m4_assert(associated(js))
    js2 => js%cast_string()
    m4_assert(associated(js2))
    value = js2%value
  end function

  function json_object_get_array(this, name) result(value)
    !! get json_array
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    type(json_array), pointer      :: value

    class(json), pointer :: js

    js => this%get_json(name)
    m4_assert(associated(js))
    value => js%cast_array()
  end function

  function json_object_get_object(this, name) result(value)
    !! get json_object
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name
    type(json_object), pointer     :: value

    class(json), pointer :: js

    js => this%get_json(name)
    m4_assert(associated(js))
    value => js%cast_object()
  end function

  subroutine json_object_set_null(this, name, p)
    !! set/add null property
    class(json_object),                 intent(inout) :: this
    character(*),                       intent(in)    :: name
      !! property name
    type(json_null), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_null object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_null),          pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_null
      allocate (json_null :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_null()
      if (present(p)) p => p_
    else ! new property
      call this%add_null(name, p = p)
    end if
  end subroutine

  subroutine json_object_set_log(this, name, value, p)
    !! set/add logical property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    logical,                           intent(in)    :: value
      !! property value
    type(json_log), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_log object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_log),           pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_log
      allocate (json_log :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_log()
      p_%value = value
      if (present(p)) p => p_
    else ! new property
      call this%add_log(name, value, p = p)
    end if
  end subroutine

  subroutine json_object_set_int32(this, name, value, p)
    !! set/add integer property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    integer(kind=int32),               intent(in)    :: value
      !! property value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_int),           pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_int
      allocate (json_int :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_int()
      p_%value = value
      if (present(p)) p => p_
    else ! new property
      call this%add_int(name, value, p = p)
    end if
  end subroutine

  subroutine json_object_set_int64(this, name, value, p)
    !! set/add integer property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    integer(kind=int64),               intent(in)    :: value
      !! property value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_int),           pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_int
      allocate (json_int :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_int()
      p_%value = value
      if (present(p)) p => p_
    else ! new property
      call this%add_int(name, value, p = p)
    end if
  end subroutine

  subroutine json_object_set_real(this, name, value, p)
    !! set/add real property
    class(json_object),                 intent(inout) :: this
    character(*),                       intent(in)    :: name
      !! property name
    real,                               intent(in)    :: value
      !! property value
    type(json_real), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_real object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_real),          pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_real
      allocate (json_real :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_real()
      p_%value = value
      if (present(p)) p => p_
    else ! new property
      call this%add_real(name, value, p = p)
    end if
  end subroutine

  subroutine json_object_set_string(this, name, value, p)
    !! set/add string property
    class(json_object),                   intent(inout) :: this
    character(*),                         intent(in)    :: name
      !! property name
    character(*),                         intent(in)    :: value
      !! property value
    type(json_string), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_string object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_string),        pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_string
      allocate (json_string :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_string()
      p_%value = value
      if (present(p)) p => p_
    else ! new property
      call this%add_string(name, value, p = p)
    end if
  end subroutine

  subroutine json_object_set_array(this, name, p)
    !! set/add array property
    class(json_object),                  intent(inout) :: this
    character(*),                        intent(in)    :: name
      !! property name
    type(json_array), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new empty json_array object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_array),         pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_array
      allocate (json_array :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_array()
      call p_%init()
      if (present(p)) p => p_
    else ! new property
      call this%add_array(name, p = p)
    end if
  end subroutine

  subroutine json_object_set_object(this, name, p)
    !! set/add array property
    class(json_object),                   intent(inout) :: this
    character(*),                         intent(in)    :: name
      !! property name
    type(json_object), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new empty json_object

    integer                           :: i
    type(mapnode_string_int), pointer :: nd
    class(json),              pointer :: js
    type(json_object),        pointer :: p_

    ! lookup property name
    nd => this%property_map%find(new_string(name))
    if (associated(nd)) then ! property with this name already exists
      i = nd%value

      ! delete old data
      if (associated(this%properties%d(i)%p)) then
        call this%properties%d(i)%p%destruct()
        deallocate (this%properties%d(i)%p)
      end if

      ! reallocate to json_object
      allocate (json_object :: js)
      this%properties%d(i)%p => js
      p_ => js%cast_object()
      call p_%init()
      if (present(p)) p => p_
    else ! new property
      call this%add_object(name, p = p)
    end if
  end subroutine

  subroutine json_object_add_null(this, name, p)
    !! add null property
    class(json_object),                 intent(inout) :: this
    character(*),                       intent(in)    :: name
      !! property name
    type(json_null), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_null object

    class(json),     pointer :: js
    type(json_null), pointer :: p_

    allocate (json_null :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_null()
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_log(this, name, value, p)
    !! add logical property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    logical,                           intent(in)    :: value
      !! property value
    type(json_log), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_log object

    class(json),    pointer :: js
    type(json_log), pointer :: p_

    allocate (json_log :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_log()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_int32(this, name, value, p)
    !! add integer property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    integer(kind=int32),               intent(in)    :: value
      !! property value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    class(json),    pointer :: js
    type(json_int), pointer :: p_

    allocate (json_int :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_int()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_int64(this, name, value, p)
    !! add integer property
    class(json_object),                intent(inout) :: this
    character(*),                      intent(in)    :: name
      !! property name
    integer(kind=int64),               intent(in)    :: value
      !! property value
    type(json_int), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_int object

    class(json),    pointer :: js
    type(json_int), pointer :: p_

    allocate (json_int :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_int()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_real(this, name, value, p)
    !! add real property
    class(json_object),                 intent(inout) :: this
    character(*),                       intent(in)    :: name
      !! property name
    real,                               intent(in)    :: value
      !! property value
    type(json_real), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_real object

    class(json),     pointer :: js
    type(json_real), pointer :: p_

    allocate (json_real :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_real()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_string(this, name, value, p)
    !! add string property
    class(json_object),                   intent(inout) :: this
    character(*),                         intent(in)    :: name
      !! property name
    character(*),                         intent(in)    :: value
      !! property value
    type(json_string), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new json_string object

    class(json),       pointer :: js
    type(json_string), pointer :: p_

    allocate (json_string :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_string()
    p_%value = value
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_array(this, name, p)
    !! add array property
    class(json_object),                  intent(inout) :: this
    character(*),                        intent(in)    :: name
      !! property name
    type(json_array), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new empty json_array object

    class(json),      pointer :: js
    type(json_array), pointer :: p_

    allocate (json_array :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_array()
    call p_%init()
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_object(this, name, p)
    !! add object property
    class(json_object),                   intent(inout) :: this
    character(*),                         intent(in)    :: name
      !! property name
    type(json_object), optional, pointer, intent(out)   :: p
      !! optional: output pointer to new empty json_object object

    class(json),       pointer :: js
    type(json_object), pointer :: p_

    allocate (json_object :: js)
    call this%add_property(name, js%get_ptr())
    p_ => js%cast_object()
    call p_%init()
    if (present(p)) p => p_
  end subroutine

  subroutine json_object_add_property(this, name, ptr)
    !! add allocated json value pointer to properties
    class(json_object),   intent(inout) :: this
    character(*),         intent(in)    :: name
      !! property name
    type(json_ptr),       intent(in)    :: ptr
      !! allocated json value pointer

    logical :: status

    call this%names%push(new_string(name))
    call this%properties%push(ptr)
    call this%property_map%insert(new_string(name), this%properties%n, status = status)
    if (.not. status) call program_error("property with name '"//name//"' already exists")
  end subroutine

end module
