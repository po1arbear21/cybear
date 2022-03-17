module json_m

  use error_m,  only: program_error
  use map_m,    only: map_string_int, mapnode_string_int
  use string_m, only: string, new_string
  use util_m,   only: int2str
  use vector_m, only: vector_int, vector_string

  implicit none

  private
  public :: json, json_null, json_bool, json_int, json_real, json_string, json_array, json_object
  public :: json_load, json_save, json_read, json_write

  character(*), parameter :: INDENT_INC = "  "
  character(*), parameter :: REAL_FMT   = "(ES23.16)"
  integer,      parameter :: REAL_LEN   = 24

  type, abstract :: json
  contains
    procedure :: init     => json_init
    procedure :: destruct => json_destruct
  end type

  type json_ptr
    class(json), pointer :: p => null()
  end type

#define T json_ptr
#define TT type(json_ptr)
#include "vector_def.f90.inc"

  type, extends(json) :: json_null
    !! null json value
  end type

  type, extends(json) :: json_bool
    !! logical json value
    logical :: value
  end type

  type, extends(json) :: json_int
    !! integer json number
    integer :: value
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
    procedure :: get_json => json_array_get_json
    generic   :: get      => json_array_get_null, &
      &                      json_array_get_bool, &
      &                      json_array_get_int, &
      &                      json_array_get_real, &
      &                      json_array_get_string, &
      &                      json_array_get_array, &
      &                      json_array_get_object
    generic   :: set      => json_array_set_null, &
      &                      json_array_set_bool, &
      &                      json_array_set_int, &
      &                      json_array_set_real, &
      &                      json_array_set_string, &
      &                      json_array_set_array, &
      &                      json_array_set_object
    generic   :: add      => json_array_add_null, &
      &                      json_array_add_bool, &
      &                      json_array_add_int, &
      &                      json_array_add_real, &
      &                      json_array_add_string, &
      &                      json_array_add_array, &
      &                      json_array_add_object

    procedure, private :: json_array_get_json
    procedure, private :: json_array_get_null
    procedure, private :: json_array_get_bool
    procedure, private :: json_array_get_int
    procedure, private :: json_array_get_real
    procedure, private :: json_array_get_string
    procedure, private :: json_array_get_array
    procedure, private :: json_array_get_object
    procedure, private :: json_array_set_null
    procedure, private :: json_array_set_bool
    procedure, private :: json_array_set_int
    procedure, private :: json_array_set_real
    procedure, private :: json_array_set_string
    procedure, private :: json_array_set_array
    procedure, private :: json_array_set_object
    procedure, private :: json_array_add_null
    procedure, private :: json_array_add_bool
    procedure, private :: json_array_add_int
    procedure, private :: json_array_add_real
    procedure, private :: json_array_add_string
    procedure, private :: json_array_add_array
    procedure, private :: json_array_add_object
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
    procedure :: get_json => json_object_get_json
    generic   :: get      => json_object_get_null, &
      &                      json_object_get_bool, &
      &                      json_object_get_int, &
      &                      json_object_get_real, &
      &                      json_object_get_string, &
      &                      json_object_get_array, &
      &                      json_object_get_object
    generic   :: set      => json_object_set_null, &
      &                      json_object_set_bool, &
      &                      json_object_set_int, &
      &                      json_object_set_real, &
      &                      json_object_set_string, &
      &                      json_object_set_array, &
      &                      json_object_set_object
    generic   :: add      => json_object_add_null, &
      &                      json_object_add_bool, &
      &                      json_object_add_int, &
      &                      json_object_add_real, &
      &                      json_object_add_string, &
      &                      json_object_add_array, &
      &                      json_object_add_object

    procedure, private :: json_object_get_json
    procedure, private :: json_object_get_null
    procedure, private :: json_object_get_bool
    procedure, private :: json_object_get_int
    procedure, private :: json_object_get_real
    procedure, private :: json_object_get_string
    procedure, private :: json_object_get_array
    procedure, private :: json_object_get_object
    procedure, private :: json_object_set_null
    procedure, private :: json_object_set_bool
    procedure, private :: json_object_set_int
    procedure, private :: json_object_set_real
    procedure, private :: json_object_set_string
    procedure, private :: json_object_set_array
    procedure, private :: json_object_set_object
    procedure, private :: json_object_add_null
    procedure, private :: json_object_add_bool
    procedure, private :: json_object_add_int
    procedure, private :: json_object_add_real
    procedure, private :: json_object_add_string
    procedure, private :: json_object_add_array
    procedure, private :: json_object_add_object

    procedure, private :: add_property => json_object_add_property
  end type

contains

#define T json_ptr
#define TT type(json_ptr)
#include "vector_imp.f90.inc"

  subroutine json_load(js, file)
    !! load json data object from file
    class(json), pointer, intent(out) :: js
      !! output pointer to new json data object
    character(*),         intent(in)  :: file
      !! file name

    character(80)             :: iomsg
    character(:), allocatable :: buf
    integer                   :: funit, iostat, n

    ! try to open file
    open (newunit = funit, file = file, access = "stream", form = "unformatted", status = "old", action = "read", iostat = iostat, iomsg = iomsg)
    if (iostat /= 0) call program_error("Error while opening file: " // iomsg)

    ! read whole file into buffer
    inquire (file = file, size = n)
    allocate (character(n) :: buf)
    read (funit, iostat = iostat, iomsg = iomsg) buf
    close (funit)
    if (iostat /= 0) call program_error("Error while reading file: " // iomsg)

    ! parse
    call json_read(js, buf)
  end subroutine

  subroutine json_save(js, file)
    !! save json data object to file
    class(json),  intent(in) :: js
    character(*), intent(in) :: file

    character(80)             :: iomsg
    character(:), allocatable :: buf
    integer                   :: funit, iostat

    ! write json into buffer
    call json_write(js, buf)

    ! try to open file
    open (newunit = funit, file = file, access = "stream", form = "unformatted", status = "replace", action = "write", iostat = iostat, iomsg = iomsg)
    if (iostat /= 0) call program_error("Error while opening file: " // iomsg)

    ! write to file
    write (funit, iostat = iostat, iomsg = iomsg) buf
    if (iostat /= 0) call program_error("Error while writing to file: " // iomsg)

    ! newline at the end of file
    write (funit) char(10)

    ! close file
    close(funit)
  end subroutine

  subroutine json_read(js, str)
    !! parse string to json data object
    class(json), pointer, intent(out) :: js
      !! output pointer to new json data object
    character(*),         intent(in)  :: str
      !! json data string

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

    ! recursive parsing
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
        allocate (json_bool :: p)
        select type (p)
        type is (json_bool)
          p%value = .true.
        end select

      case ('f')
        if (        i1 /=  i0 + 4) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        if (str(i0:i1) /= "false") call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
        allocate (json_bool :: p)
        select type (p)
        type is (json_bool)
          p%value = .false.
        end select

      case ('-', '0':'9')
        l = scan(str(i0:i1), '.') + scan(str(i0:i1), 'd') + scan(str(i0:i1), 'e')

        if (l == 0) then
          read (str(i0:i1), *, iostat = iostat) i
          if (iostat /= 0) call program_error("Error in line " // int2str(line) // ": unable to parse value '" // str(i0:i1) // "'")
          allocate (json_int :: p)
          select type (p)
          type is (json_int)
            p%value = i
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

  recursive subroutine json_write(js, str, indent)
    class(json),               intent(in)  :: js
      !! json data object
    character(:), allocatable, intent(out) :: str
      !! output json data string
    character(*), optional,    intent(in)  :: indent
      !! current indentation

    character(:), allocatable :: tmp, indent_
    integer                   :: i

    indent_ = ""
    if (present(indent)) indent_ = indent

    select type (js)
    type is (json_null)
      str = "null"
    type is (json_bool)
      if (js%value) then
        str = "true"
      else
        str = "false"
      end if
    type is (json_int)
      allocate (character(32) :: tmp)
      write (tmp, "(I0)") js%value
      str = trim(tmp)
      print *, 'int', js%value
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
          call json_write(js%values%d(i)%p, tmp, indent_ // INDENT_INC)
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
          call json_write(js%properties%d(i)%p, tmp, indent_ // INDENT_INC)
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
    !! initialize json
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
    !! destruct json
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

  subroutine json_array_get_json(this, idx, js)
    !! get idx-th element from json array
    class(json_array),    intent(in)  :: this
    integer,              intent(in)  :: idx
      !! index
    class(json), pointer, intent(out) :: js
      !! output pointer to json value

    js => this%values%d(idx)%p
  end subroutine

  subroutine json_array_get_null(this, idx)
    !! get idx-th element from json array, error if not json_null
    class(json_array), intent(in) :: this
    integer,           intent(in) :: idx
      !! index

    select type (p => this%values%d(idx)%p)
    type is (json_null)
    class default
      call program_error(int2str(idx) // "-th element is not of type json_null")
    end select
  end subroutine

  subroutine json_array_get_bool(this, idx, value)
    !! get idx-th element as logical from json array, error if not json_bool
    class(json_array), intent(in)  :: this
    integer,           intent(in)  :: idx
      !! index
    logical,           intent(out) :: value
      !! output logical element from array

    select type (p => this%values%d(idx)%p)
    type is (json_bool)
      value = p%value
    class default
      call program_error(int2str(idx) // "-th element is not of type json_bool")
    end select
  end subroutine

  subroutine json_array_get_int(this, idx, value)
    !! get idx-th element as integer from json array, error if not json_int
    class(json_array), intent(in)  :: this
    integer,           intent(in)  :: idx
      !! index
    integer,           intent(out) :: value
      !! output integer element from array

    select type (p => this%values%d(idx)%p)
    type is (json_int)
      value = p%value
    class default
      call program_error(int2str(idx) // "-th element is not of type json_int")
    end select
  end subroutine

  subroutine json_array_get_real(this, idx, value)
    !! get idx-th element as real from json array, error if not json_real
    class(json_array), intent(in)  :: this
    integer,           intent(in)  :: idx
      !! index
    real,              intent(out) :: value
      !! output real element from array

    select type (p => this%values%d(idx)%p)
    type is (json_real)
      value = p%value
    class default
      call program_error(int2str(idx) // "-th element is not of type json_real")
    end select
  end subroutine

  subroutine json_array_get_string(this, idx, value)
    !! get idx-th element as string from json array, error if not json_string
    class(json_array),         intent(in)  :: this
    integer,                   intent(in)  :: idx
      !! index
    character(:), allocatable, intent(out) :: value
      !! output string element from array

    select type (p => this%values%d(idx)%p)
    type is (json_string)
      value = p%value
    class default
      call program_error(int2str(idx) // "-th element is not of type json_string")
    end select
  end subroutine

  subroutine json_array_get_array(this, idx, ar)
    !! get idx-th element as array from json array, error if not json_array
    class(json_array),         intent(in)  :: this
    integer,                   intent(in)  :: idx
      !! index
    type(json_array), pointer, intent(out) :: ar
      !! output array element from array

    select type (p => this%values%d(idx)%p)
    type is (json_array)
      ar => p
    class default
      call program_error(int2str(idx) // "-th element is not of type json_array")
    end select
  end subroutine

  subroutine json_array_get_object(this, idx, obj)
    !! get idx-th element as object from json array, error if not json_object
    class(json_array),          intent(in)  :: this
    integer,                    intent(in)  :: idx
      !! index
    type(json_object), pointer, intent(out) :: obj
      !! output object element from array

    select type (p => this%values%d(idx)%p)
    type is (json_object)
      obj => p
    class default
      call program_error(int2str(idx) // "-th element is not of type json_object")
    end select
  end subroutine

  subroutine json_array_set_null(this, idx)
    !! set idx-th element of json array to null
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: idx
      !! index

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    allocate (json_null :: this%values%d(idx)%p)
  end subroutine

  subroutine json_array_set_bool(this, idx, value)
    !! set idx-th element of json array to logical
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: idx
      !! index
    logical,           intent(in)    :: value
      !! logical value

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    allocate (json_bool :: this%values%d(idx)%p)
    select type (p => this%values%d(idx)%p)
    type is (json_bool)
      p%value = value
    end select
  end subroutine

  subroutine json_array_set_int(this, idx, value)
    !! set idx-th element of json array to integer
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: idx
      !! index
    integer,           intent(in)    :: value
      !! integer value

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    allocate (json_int :: this%values%d(idx)%p)
    select type (p => this%values%d(idx)%p)
    type is (json_int)
      p%value = value
    end select
  end subroutine

  subroutine json_array_set_real(this, idx, value)
    !! set idx-th element of json array to real
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: idx
      !! index
    real,              intent(in)    :: value
      !! real value

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    allocate (json_real :: this%values%d(idx)%p)
    select type (p => this%values%d(idx)%p)
    type is (json_real)
      p%value = value
    end select
  end subroutine

  subroutine json_array_set_string(this, idx, value)
    !! set idx-th element of json array to string
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: idx
      !! index
    character(*),      intent(in)    :: value
      !! string

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    allocate (json_string :: this%values%d(idx)%p)
    select type (p => this%values%d(idx)%p)
    type is (json_string)
      p%value = value
    end select
  end subroutine

  subroutine json_array_set_array(this, idx, ar)
    !! set idx-th element of json array to array
    class(json_array),         intent(inout) :: this
    integer,                   intent(in)    :: idx
      !! index
    type(json_array), pointer, intent(in)    :: ar
      !! pointer to json array

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    this%values%d(idx)%p => ar
  end subroutine

  subroutine json_array_set_object(this, idx, obj)
    !! set idx-th element of json array to object
    class(json_array),          intent(inout) :: this
    integer,                    intent(in)    :: idx
      !! index
    type(json_object), pointer, intent(in)    :: obj
    !! pointer to json array

    if (associated(this%values%d(idx)%p)) then
      call this%values%d(idx)%p%destruct()
      deallocate (this%values%d(idx)%p)
    end if

    this%values%d(idx)%p => obj
  end subroutine

  subroutine json_array_add_null(this)
    !! add null value to json array
    class(json_array), intent(inout) :: this

    type(json_ptr) :: ptr

    allocate (json_null :: ptr%p)
    call this%values%push(ptr)
  end subroutine

  subroutine json_array_add_bool(this, value)
    !! add boolean to json array
    class(json_array), intent(inout) :: this
    logical,           intent(in)    :: value
      !! boolean value to add

    type(json_ptr) :: ptr

    allocate (json_bool :: ptr%p)
    call this%values%push(ptr)

    select type (p => ptr%p)
    type is (json_bool)
      p%value = value
    end select
  end subroutine

  subroutine json_array_add_int(this, value)
    !! add integer to json array
    class(json_array), intent(inout) :: this
    integer,           intent(in)    :: value
      !! integer value to add

    type(json_ptr) :: ptr

    allocate (json_int :: ptr%p)
    call this%values%push(ptr)

    select type (p => ptr%p)
    type is (json_int)
      p%value = value
    end select
  end subroutine

  subroutine json_array_add_real(this, value)
    !! add real to json array
    class(json_array), intent(inout) :: this
    real,              intent(in)    :: value
      !! real value to add

    type(json_ptr) :: ptr

    allocate (json_real :: ptr%p)
    call this%values%push(ptr)

    select type (p => ptr%p)
    type is (json_real)
      p%value = value
    end select
  end subroutine

  subroutine json_array_add_string(this, value)
    !! add string to json array
    class(json_array), intent(inout) :: this
    character(*),      intent(in)    :: value
      !! string to add

    type(json_ptr) :: ptr

    allocate (json_string :: ptr%p)
    call this%values%push(ptr)

    select type (p => ptr%p)
    type is (json_string)
      p%value = value
    end select
  end subroutine

  subroutine json_array_add_array(this, ar)
    !! add array to json array
    class(json_array),         intent(inout) :: this
    type(json_array), pointer, intent(in)    :: ar
      !! pointer to json array

    type(json_ptr) :: ptr

    ptr%p => ar
    call this%values%push(ptr)
  end subroutine

  subroutine json_array_add_object(this, obj)
    !! add object to json array
    class(json_array),          intent(inout) :: this
    type(json_object), pointer, intent(in)    :: obj
      !! pointer to json object

    type(json_ptr) :: ptr

    ptr%p => obj
    call this%values%push(ptr)
  end subroutine

  subroutine json_object_get_json(this, name, js)
    !! get json value from json object
    class(json_object),   intent(in)  :: this
    character(*),         intent(in)  :: name
      !! property name
    class(json), pointer, intent(out) :: js
      !! output pointer to json value

    integer :: idx

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p
  end subroutine

  subroutine json_object_get_null(this, name)
    !! get null from json object
    class(json_object), intent(in) :: this
    character(*),       intent(in) :: name
      !! property name

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_null)
    class default
      call program_error('property with name "' // name // '" is not of type json_null')
    end select
  end subroutine

  subroutine json_object_get_bool(this, name, value)
    !! get logical from json object
    class(json_object), intent(in)  :: this
    character(*),       intent(in)  :: name
      !! property name
    logical,            intent(out) :: value
      !! output logical value

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_bool)
      value = js%value
    class default
      call program_error('property with name "' // name // '" is not of type json_bool')
    end select
  end subroutine

  subroutine json_object_get_int(this, name, value)
    !! get integer from json object
    class(json_object), intent(in)  :: this
    character(*),       intent(in)  :: name
      !! property name
    integer,            intent(out) :: value
      !! output integer value

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_int)
      value = js%value
    class default
      call program_error('property with name "' // name // '" is not of type json_int')
    end select
  end subroutine

  subroutine json_object_get_real(this, name, value)
    !! get real from json object
    class(json_object), intent(in)  :: this
    character(*),       intent(in)  :: name
      !! property name
    real,               intent(out) :: value
      !! output integer value

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_real)
      value = js%value
    class default
      call program_error('property with name "' // name // '" is not of type json_real')
    end select
  end subroutine

  subroutine json_object_get_string(this, name, value)
    !! get string from json object
    class(json_object),        intent(in)  :: this
    character(*),              intent(in)  :: name
      !! property name
    character(:), allocatable, intent(out) :: value
      !! output integer value

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_string)
      value = js%value
    class default
      call program_error('property with name "' // name // '" is not of type json_string')
    end select
  end subroutine

  subroutine json_object_get_array(this, name, ar)
    !! get array from json object
    class(json_object),        intent(in)  :: this
    character(*),              intent(in)  :: name
      !! property name
    type(json_array), pointer, intent(out) :: ar
      !! output pointer to json array

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_array)
      ar => js
    class default
      call program_error('property with name "' // name // '" is not of type json_array')
    end select
  end subroutine

  subroutine json_object_get_object(this, name, obj)
    !! get object from json object
    class(json_object),         intent(in)  :: this
    character(*),               intent(in)  :: name
      !! property name
    type(json_object), pointer, intent(out) :: obj
      !! output pointer to json array

    integer              :: idx
    class(json), pointer :: js

    idx = this%property_map%get(new_string(name))
    js => this%properties%d(idx)%p

    select type (js)
    type is (json_object)
      obj => js
    class default
      call program_error('property with name "' // name // '" is not of type json_object')
    end select
  end subroutine

  subroutine json_object_set_null(this, name)
    !! set json object property to null
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      allocate (json_null :: this%properties%d(idx)%p)
    else
      call this%add(name)
    end if
  end subroutine

  subroutine json_object_set_bool(this, name, value)
    !! set json object property to logical value
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name
    logical,            intent(in)    :: value
      !! logical value

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      allocate (json_bool :: this%properties%d(idx)%p)
      select type (p => this%properties%d(idx)%p)
      type is (json_bool)
        p%value = value
      end select
    else
      call this%add(name, value)
    end if
  end subroutine

  subroutine json_object_set_int(this, name, value)
    !! set json object property to integer value
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name
    integer,            intent(in)    :: value
      !! integer value

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      allocate (json_int :: this%properties%d(idx)%p)
      select type (p => this%properties%d(idx)%p)
      type is (json_int)
        p%value = value
      end select
    else
      call this%add(name, value)
    end if
  end subroutine

  subroutine json_object_set_real(this, name, value)
    !! set json object property to real value
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name
    real,               intent(in)    :: value
      !! real value

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      allocate (json_real :: this%properties%d(idx)%p)
      select type (p => this%properties%d(idx)%p)
      type is (json_real)
        p%value = value
      end select
    else
      call this%add(name, value)
    end if
  end subroutine

  subroutine json_object_set_string(this, name, value)
    !! set json object property to string
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name
    character(*),       intent(in)    :: value
      !! string

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      allocate (json_string :: this%properties%d(idx)%p)
      select type (p => this%properties%d(idx)%p)
      type is (json_string)
        p%value = value
      end select
    else
      call this%add(name, value)
    end if
  end subroutine

  subroutine json_object_set_array(this, name, ar)
    !! set json object property to array
    class(json_object),        intent(inout) :: this
    character(*),              intent(in)    :: name
      !! property name
    type(json_array), pointer, intent(in)    :: ar
      !! pointer to json array

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      this%properties%d(idx)%p => ar
    else
      call this%add(name, ar)
    end if
  end subroutine

  subroutine json_object_set_object(this, name, obj)
    !! set json object property to object
    class(json_object),         intent(inout) :: this
    character(*),               intent(in)    :: name
      !! property name
    type(json_object), pointer, intent(in)    :: obj
      !! pointer to json object

    integer                           :: idx
    type(mapnode_string_int), pointer :: node

    node => this%property_map%find(new_string(name))
    if (associated(node)) then
      idx = node%value

      if (associated(this%properties%d(idx)%p)) then
        call this%properties%d(idx)%p%destruct()
        deallocate (this%properties%d(idx)%p)
      end if

      this%properties%d(idx)%p => obj
    else
      call this%add(name, obj)
    end if
  end subroutine

  subroutine json_object_add_null(this, name)
    !! add null value to json object
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name (must be unique)

    type(json_ptr) :: ptr

    allocate (json_null :: ptr%p)
    call this%add_property(name, ptr)
  end subroutine

  subroutine json_object_add_bool(this, name, value)
    !! add logical to json object
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name (must be unique)
    logical,            intent(in)    :: value
      !! property value

    type(json_ptr) :: ptr

    allocate (json_bool :: ptr%p)
    call this%add_property(name, ptr)

    select type (p => ptr%p)
    type is (json_bool)
      p%value = value
    end select
  end subroutine

  subroutine json_object_add_int(this, name, value)
    !! add integer to json object
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name (must be unique)
    integer,            intent(in)    :: value
      !! property value

    type(json_ptr) :: ptr

    allocate (json_int :: ptr%p)
    call this%add_property(name, ptr)

    select type (p => ptr%p)
    type is (json_int)
      p%value = value
    end select
  end subroutine

  subroutine json_object_add_real(this, name, value)
    !! add real to json object
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name (must be unique)
    real,               intent(in)    :: value
      !! property value

    type(json_ptr) :: ptr

    allocate (json_real :: ptr%p)
    call this%add_property(name, ptr)

    select type (p => ptr%p)
    type is (json_real)
      p%value = value
    end select
  end subroutine

  subroutine json_object_add_string(this, name, value)
    !! add string to json object
    class(json_object), intent(inout) :: this
    character(*),       intent(in)    :: name
      !! property name (must be unique)
    character(*),       intent(in)    :: value
      !! property value

    type(json_ptr) :: ptr

    allocate (json_string :: ptr%p)
    call this%add_property(name, ptr)

    select type (p => ptr%p)
    type is (json_string)
      p%value = value
    end select
  end subroutine

  subroutine json_object_add_array(this, name, ar)
    !! add array to json object
    class(json_object),        intent(inout) :: this
    character(*),              intent(in)    :: name
      !! property name (must be unique)
    type(json_array), pointer, intent(in)    :: ar
      !! pointer to json array

    type(json_ptr) :: ptr

    ptr%p => ar
    call this%add_property(name, ptr)
  end subroutine

  subroutine json_object_add_object(this, name, obj)
    !! add object to json object
    class(json_object),         intent(inout) :: this
    character(*),               intent(in)    :: name
      !! property name (must be unique)
    type(json_object), pointer, intent(in)   :: obj
      !! pointer to json object

    type(json_ptr) :: ptr

    ptr%p => obj
    call this%add_property(name, ptr)
  end subroutine

  subroutine json_object_add_property(this, name, ptr)
    !! add allocated json value pointer to properties
    class(json_object),   intent(inout) :: this
    character(*),         intent(in)    :: name
      !! property name (must be unique)
    type(json_ptr),       intent(in)    :: ptr
      !! allocated json value pointer

    logical :: status

    call this%names%push(new_string(name))
    call this%properties%push(ptr)
    call this%property_map%insert(new_string(name), this%properties%n, status = status)
    if (.not. status) call program_error("property with name '"//name//"' already exists")
  end subroutine

end module
