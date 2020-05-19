#include "macro.f90.inc"

module input_m
  use error_m
  use normalization_m
  use vector_m
  implicit none

  ! maximum line length
  integer, parameter :: MAX_LINE_LEN = 4096

  ! variable types
  integer, parameter :: INPUT_TYPE_INTEGER = 1
  integer, parameter :: INPUT_TYPE_REAL    = 2
  integer, parameter :: INPUT_TYPE_STRING  = 3
  integer, parameter :: INPUT_TYPE_LOGIC   = 4
  integer, parameter :: INPUT_NUM_TYPES    = 4

  type input_var
    character(:), allocatable :: name
      !! Variable name
    integer                   :: type
      !! Variable type
    integer                   :: size
      !! Size of variable, 1 means scalar, > 1 means array
    integer,      allocatable :: int_data(:)
      !! integer data
    real,         allocatable :: real_data(:)
      !! real data
    type(string), allocatable :: str_data(:)
      !! string data
    logical,      allocatable :: logic_data(:)
      !! logical data
    character(:), allocatable :: unit
      !! Physical unit token (for real variables)
  contains
    procedure :: init => input_var_init
  end type

#define T input_var
#define TT type(input_var)
#include "vector_def.f90.inc"

  type input_section
    character(:), allocatable :: name
      !! Section name
    type(vector_input_var)    :: vars
      !! variables in this section
  contains
    procedure :: init    => input_section_init
    procedure :: add_var => input_section_add_var
  end type

#define T input_section
#define TT type(input_section)
#include "vector_def.f90.inc"

  type input_file
    character(:), allocatable  :: name
      !! File name
    character(:), allocatable  :: default
      !! default file name
    type(vector_input_section) :: sections
      !! Sections in this file
    type(vector_input_section) :: default_sections
      !! Default values for sections
  contains
    private
    procedure :: load         => input_file_load
    procedure :: parse_line   => input_file_parse_line

    procedure :: input_file_get_int
    procedure :: input_file_get_int_arr
    procedure :: input_file_get_name_int
    procedure :: input_file_get_name_int_arr
    procedure :: input_file_get_real
    procedure :: input_file_get_real_arr
    procedure :: input_file_get_name_real
    procedure :: input_file_get_name_real_arr
    procedure :: input_file_get_string
    procedure :: input_file_get_string_arr
    procedure :: input_file_get_name_string
    procedure :: input_file_get_name_string_arr
    procedure :: input_file_get_logical
    procedure :: input_file_get_logical_arr
    procedure :: input_file_get_name_logical
    procedure :: input_file_get_name_logical_arr

    procedure, public :: init         => input_file_init
    procedure, public :: get_section  => input_file_get_section
    procedure, public :: get_sections => input_file_get_sections
    generic,   public :: get          => input_file_get_int,             &
      &                                  input_file_get_int_arr,         &
      &                                  input_file_get_name_int,        &
      &                                  input_file_get_name_int_arr,    &
      &                                  input_file_get_real,            &
      &                                  input_file_get_real_arr,        &
      &                                  input_file_get_name_real,       &
      &                                  input_file_get_name_real_arr,   &
      &                                  input_file_get_string,          &
      &                                  input_file_get_string_arr,      &
      &                                  input_file_get_name_string,     &
      &                                  input_file_get_name_string_arr, &
      &                                  input_file_get_logical,         &
      &                                  input_file_get_logical_arr,     &
      &                                  input_file_get_name_logical,    &
      &                                  input_file_get_name_logical_arr
  end type

contains

#define T input_var
#define TT type(input_var)
#include "vector_imp.f90.inc"

#define T input_section
#define TT type(input_section)
#include "vector_imp.f90.inc"

  subroutine input_var_init(this, data0, valid, err)
    !! initialize variable from data string
    class(input_var), intent(out) :: this
    character(*),     intent(in)  :: data0
      !! data string
    logical,          intent(out) :: valid
      !! valid flag
    character(*),     intent(out) :: err
      !! Error string

    ! local variables
    integer               :: i, j, k, status
    real                  :: r
    character             :: quote
    character(len(data0)) :: data, item

    ! assume no error
    valid = .true.
    err   = ""

    ! work on copy of data string
    data = data0

    ! search for equals sign
    i = scan(data, '=')
    if (i == 0) then
      valid = .false.
      err   = "Variable definition error"
      return
    end if

    ! split var name and data
    this%name = trim(adjustl(data(1:i-1)))
    data      = trim(adjustl(data(i+1:len_trim(data))))
    if (this%name == "") then
      valid = .false.
      err   = "Empty variable name"
      return
    end if
    if (data == "") then
      valid = .false.
      err   = "Empty variable data"
      return
    end if

    ! search for unit
    i = scan(data, ':')
    if (i /= 0) then
      ! implicit real
      this%type = INPUT_TYPE_REAL

      this%unit = trim(adjustl(data(i+1:len_trim(data))))
      data = trim(adjustl(data(1:i-1)))
      if (this%unit == "") then
        valid = .false.
        err   = "Empty variable unit"
        return
      end if
      if (data == "") then
        valid = .false.
        err   = "Empty variable data"
        return
      end if
    else
      this%type = 0
      this%unit = "1"
    end if

    ! get size
    quote = ""
    this%size = 1
    do i = 1, len_trim(data)
      ! check if inside of quote
      if (quote == "") then ! not inside of quote
        ! check for beginning of quote or comment
        if ((data(i:i) == "'") .or. (data(i:i) == '"')) then
          quote = data(i:i)
        elseif (data(i:i) == ",") then
          ! count comma
          this%size = this%size + 1
        end if
      elseif (data(i:i) == quote) then ! end of quote
        quote = ""
      end if
    end do

    ! get type from first item (if not already set)
    if (this%type == 0) then
      i = scan(data, ',')
      if (i == 0) then
        item = data
      else
        item = trim(adjustl(data(1:i-1)))
      end if
      if (item == "") then
        valid = .false.
        err   = "First data item is empty"
        return
      end if

      ! search for quotes
      i = scan(item, '"')
      if (i == 0) then
        i = scan(item, "'")
      end if
      DEDUCE_TYPE: if (i /= 0) then
        ! string
        this%type = INPUT_TYPE_STRING
      else

        ! only try to read as integer if there are no decimal dots or 'd', 'e'
        i = scan(item, '.') + scan(item, 'd') + scan(item, 'e')
        if (i == 0) then
          read (item, *, iostat = status) i
          if (status == 0) then
            this%type = INPUT_TYPE_INTEGER
            exit DEDUCE_TYPE
          end if
        end if

        ! try reading as real
        read (item, *, iostat = status) r
        if (status == 0) then
          this%type = INPUT_TYPE_REAL
          exit DEDUCE_TYPE
        end if

        ! try reading as bool
        if ((trim(item) == "true") .or. (trim(item) == "false")) then
          this%type = INPUT_TYPE_LOGIC
          exit DEDUCE_TYPE
        end if

        ! we reached this far -> couldnt determine data type -> error
        valid = .false.
        err   = "Cannot deduce type of variable"
        return
      end if DEDUCE_TYPE
    end if

    ! read in data
    select case (this%type)
    case (INPUT_TYPE_INTEGER)
      allocate (this%int_data(this%size))

      do i = 1, this%size
        j = scan(data(1:len_trim(data)), ',')
        if (j == 0) then
          item = data
          data = ""
        else
          item = trim(adjustl(data(1:j-1)))
          data = trim(adjustl(data(j+1:len_trim(data))))
        end if

        ! convert to integer
        read (item, *, iostat = status) this%int_data(i)
        if (status /= 0) then
          valid = .false.
          err   = "invalid integer data"
          return
        end if
      end do

    case (INPUT_TYPE_REAL)
      allocate (this%real_data(this%size))

      do i = 1, this%size
        j = scan(data(1:len_trim(data)), ',')
        if (j == 0) then
          item = data
          data = ""
        else
          item = trim(adjustl(data(1:j-1)))
          data = trim(adjustl(data(j+1:len_trim(data))))
        end if

        ! convert to real
        read (item, *, iostat = status) this%real_data(i)
        if (status /= 0) then
          valid = .false.
          err   = "invalid real data"
          return
        end if
      end do

    case (INPUT_TYPE_STRING)
      allocate (this%str_data(this%size))

      quote = ""
      j = 0
      k = 1
      do i = 1, len_trim(data)
        ! check if inside of quote
        if (quote == "") then ! not inside of quote
          ! only quote, comma and space allowed
          if ((data(i:i) == "'") .or. (data(i:i) == '"')) then
            if (j == -1) then
              valid = .false.
              err   = "comma needed between two strings"
              return
            end if

            ! beginning of quote
            quote = data(i:i)
            j = i
          elseif (data(i:i) == ",") then
            ! next variable
            k = k + 1
            j = 0
          elseif (data(i:i) /= " ") then
            valid = .false.
            err   = "invalid string data"
            return
          end if
        elseif (data(i:i) == quote) then ! end of quote
          ! string finished
          quote = ""
          this%str_data(k)%s = data(j+1:i-1)
          j = -1
        end if
      end do

    case (INPUT_TYPE_LOGIC)
      allocate (this%logic_data(this%size))

      do i = 1, this%size
        j = scan(data(1:len_trim(data)), ',')
        if (j == 0) then
          item = data
          data = ""
        else
          item = trim(adjustl(data(1:j-1)))
          data = trim(adjustl(data(j+1:len_trim(data))))
        end if

        ! convert to logical
        if (item == "true") then
          this%logic_data(i) = .true.
        elseif (item == "false") then
          this%logic_data(i) = .false.
        else
          valid = .false.
          err   = "invalid logical data"
          return
        end if
      end do
    end select
  end subroutine

  subroutine input_section_init(this, name)
    !! initialize input section
    class(input_section), intent(out) :: this
    character(*),         intent(in)  :: name
      !! section name

    ! init members
    this%name = name
    call this%vars%init(0, c = 8)
  end subroutine

  subroutine input_section_add_var(this, v)
    !! add/replace variable in section
    class(input_section), intent(inout) :: this
    type(input_var),      intent(in)    :: v
      !! variable to be added/replaced

    ! local variables
    integer :: i

    ! search for existing var
    do i = 1, this%vars%n
      if (this%vars%d(i)%name == v%name) then
        this%vars%d(i) = v
        return
      end if
    end do

    ! add new var
    call this%vars%push(v)
  end subroutine

  subroutine input_file_init(this, name, default)
    !! initialize input file
    class(input_file),      intent(out) :: this
    character(*),           intent(in)  :: name
      !! filename
    character(*), optional, intent(in)  :: default
      !! filename for default sections

    ! save names
    this%name    = name
    this%default = ""
    if (present(default)) this%default = default

    ! load default file
    if (present(default)) then
      call this%load(this%default_sections, default)
    else
      call this%default_sections%init(0, c = 8)
    end if

    ! load file
    call this%load(this%sections, name, default_sections = this%default_sections)
  end subroutine

  subroutine input_file_load(this, sections, name, default_sections)
    !! load file and setup sections
    class(input_file),                    intent(inout) :: this
    type(vector_input_section),           intent(out)   :: sections
      !! setup sections
    character(*),                         intent(in)    :: name
      !! filename to load
    type(vector_input_section), optional, intent(in)    :: default_sections
      !! sections with default values

    ! local variables
    integer                   :: funit, line_number, status
    logical                   :: valid
    character(10)             :: lfmt
    character(:), allocatable :: line, err
    type(vector_string)       :: cont_lines

    ! line format
    write(lfmt, "(A,I0,A)") "(A", MAX_LINE_LEN, ")"

    ! continuation lines vector
    call cont_lines%init(0, c = 8)

    ! init section vector
    call sections%init(0, c = 8)

    ! allocate memory for line
    allocate (character(MAX_LINE_LEN) :: line)

    ! open file
    open(newunit = funit, file = trim(name), status = "old", action = "read")

    ! start at the top
    line_number = 0

    ! read and parse lines
    do
      ! read one line
      read(funit, fmt = lfmt, iostat = status) line
      if (status < 0) then
        ! end of file reached
        close(funit)
        exit
      elseif (status > 0) then
        ! IO error
        close(funit)
        call program_error("IO Error")
      end if
      line_number = line_number + 1

      ! parse line
      call this%parse_line(sections, trim(line), cont_lines, valid, err, default_sections = default_sections)

      ! check if line was valid
      if (.not. valid) then
        close(funit)
        print "(1A, 1I0)", "Input error in line number ", line_number
        call program_error(err)
      end if
    end do

  end subroutine

  subroutine input_file_parse_line(this, sections, line0, cont_lines, valid, err, default_sections)
    !! Parse a single line from input file, update sections
    class(input_file),                    intent(inout) :: this
    type(vector_input_section),           intent(inout) :: sections
      !! Input sections
    character(*),                         intent(in)    :: line0
      !! Line from input file
    type(vector_string),                  intent(inout) :: cont_lines
      !! Collection of lines belonging together
    logical,                              intent(out)   :: valid
      !! Output if parsing was successful (valid syntax)
    character(:), allocatable,            intent(out)   :: err
      !! Error string
    type(vector_input_section), optional, intent(in)    :: default_sections
      !! Input sections

    ! local variables
    character(:), allocatable  :: line, name
    character(1)               :: quote
    type(string)               :: str
    integer                    :: i, j, line_len

    IGNORE(this)

    ! assume valid syntax; change later if invalid
    valid = .true.
    err   = ""

    ! work on copy; remove leading and trailing spaces
    allocate (character(0) :: line) ! remove gfortran warning
    line     = trim(adjustl(line0))
    line_len = len_trim(line)

    ! remove comments
    block
      quote = ""

      ! go through line and search for comment character
      do i = 1, line_len
        ! check if inside of quote
        if (quote == "") then ! not inside of quote
          ! check for beginning of quote or comment
          if ((line(i:i) == "'") .or. (line(i:i) == '"')) then
            quote = line(i:i)
          elseif (line(i:i) == "!") then
            ! trim line up to comment character
            line     = trim(line(1:i-1))
            line_len = len_trim(line)
            exit
          end if
        elseif (line(i:i) == quote) then ! end of quote
          quote = ""
        end if
      end do

      ! make sure line does not end inside quote
      if (quote /= "") then
        valid = .false.
        err   = "quote must be contained in single line"
        return
      end if
    end block

    ! check for line continuation symbol
    if (line_len > 0) then
      if (line(line_len:line_len) == '&') then
        line = trim(adjustl(line(1:line_len-1)))
        if (line /= "") then
          str%s = trim(line)
          call cont_lines%push(str)
        end if
        return
      end if
    end if

    ! check for previous line continuation
    if (cont_lines%n > 0) then ! concat lines in single giant string
      block
        character(:), allocatable :: concat_line
        type(input_var)           :: v

        ! add current line
        str%s = trim(line)
        call cont_lines%push(str)

        ! count length and allocate concat line
        j = 0
        do i = 1, cont_lines%n
          j = j + len_trim(cont_lines%d(i)%s)
        end do
        allocate (character(j) :: concat_line)

        j = 0
        do i = 1, cont_lines%n
          concat_line(j+1:j+len_trim(cont_lines%d(i)%s)) = trim(cont_lines%d(i)%s)
          j = j + len_trim(cont_lines%d(i)%s)
        end do

        call v%init(concat_line, valid, err)
        if (valid) then
          ! if no sections, add one without name
          if (sections%n == 0) then
            call sections%resize(1)
            allocate (character(0) :: name)
            call sections%d(1)%init(name)
          end if

          ! add variable to current section
          call sections%d(sections%n)%add_var(v)
        end if
      end block

      ! remove lines
      call cont_lines%resize(0)
      return
    end if

    ! skip empty lines
    if (line == "") return

    ! check for new section
    block
      character(:), allocatable :: section_name
      type(input_section)       :: sect

      if (line(1:1) == '[') then ! new section
        ! line must end with ]
        if (line(line_len:line_len) /= ']') then
          valid = .false.
          err   = "Section name must end with square bracket"
          return
        end if

        ! get name of section
        allocate (character(0) :: section_name) ! remove gfortran warning
        section_name = trim(adjustl(line(2:line_len-1)))

        ! create new section
        if (present(default_sections)) then
          ! search for section name
          do i = 1, default_sections%n
            if (default_sections%d(i)%name == section_name) then
              ! set section values to default
              sect = default_sections%d(i)
              exit
            end if
          end do
          if (i == default_sections%n + 1) call sect%init(section_name)
        else
          call sect%init(section_name)
        end if
        call sections%push(sect)
        return
      end if
    end block

    ! variable definition (in one line)
    block
      type(input_var) :: v

      call v%init(line, valid, err)
      if (valid) then
        ! if no sections, add one without name
        if (sections%n == 0) then
          call sections%resize(1)
          allocate (character(0) :: name)
          call sections%d(1)%init(name)
        end if

        ! add variable to current section
        call sections%d(sections%n)%add_var(v)
      end if
    end block
  end subroutine

  subroutine input_file_get_section(this, section_name, section_id, status)
    !! get index of section with this name
    class(input_file), intent(in)  :: this
    character(*),      intent(in)  :: section_name
      !! name of section to search
    integer,           intent(out) :: section_id
      !! output section index
    integer, optional, intent(out) :: status
      !! optional: output 0 on success; -1 if not found; >0 if not unique (if not present: error if not found or not unique)

    ! local variables
    integer :: i, st

    st = -1
    do i = 1, this%sections%n
      if (this%sections%d(i)%name == section_name) then
        section_id = i
        st         = st + 1
      end if
    end do

    if (present(status)) status = st

    if (st < 0) then
      if (.not. present(status)) then
        print *, section_name
        call program_error("section name not found")
      end if
    elseif (st > 0) then
      if (.not. present(status)) then
        print *, section_name
        call program_error("multiple sections with this name found")
      end if
    end if
  end subroutine

  subroutine input_file_get_sections(this, section_name, section_ids)
    !! get all section indices with this name
    class(input_file),    intent(in)  :: this
    character(*),         intent(in)  :: section_name
      !! name of section to search
    integer, allocatable, intent(out) :: section_ids(:)
      !! output indices

    ! local variables
    integer          :: i
    type(vector_int) :: sv

    ! collect sections with this name
    call sv%init(0, c = this%sections%n)
    do i = 1, this%sections%n
      if (this%sections%d(i)%name == section_name) then
        call sv%push(i)
      end if
    end do

    ! output section ids
    allocate (section_ids(sv%n), source = sv%d(1:sv%n))
  end subroutine

  subroutine input_file_get_int(this, section_id, name, value, status)
    !! get scalar integer value
    class(input_file), intent(in)  :: this
    integer,           intent(in)  :: section_id
      !! section index
    character(*),      intent(in)  :: name
      !! variable name
    integer,           intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%size /= 1) call program_error("variable is not scalar")
            if (var%type /= INPUT_TYPE_INTEGER) call program_error("variable is not an integer")
            value = var%int_data(1)
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_int_arr(this, section_id, name, values, status)
    !! get array of integers
    class(input_file),    intent(in)  :: this
    integer,              intent(in)  :: section_id
      !! section index
    character(*),         intent(in)  :: name
      !! variable name
    integer, allocatable, intent(out) :: values(:)
      !! output values
    logical, optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%type /= INPUT_TYPE_INTEGER) call program_error("variable is not an integer")
            allocate (values(size(var%int_data)), source = var%int_data)
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_name_int(this, section_name, name, value, status)
    !! get scalar integer value, provide section name instead of index
    class(input_file), intent(in)  :: this
    character(*),      intent(in)  :: section_name
      !! section name
    character(*),      intent(in)  :: name
      !! variable name
    integer,           intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, value, status = status)
  end subroutine

  subroutine input_file_get_name_int_arr(this, section_name, name, values, status)
    !! get array of integers, provide section name instead of index
    class(input_file),    intent(in)  :: this
    character(*),         intent(in)  :: section_name
      !! section name
    character(*),         intent(in)  :: name
      !! variable name
    integer, allocatable, intent(out) :: values(:)
      !! output values
    logical, optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, values, status = status)
  end subroutine

  subroutine input_file_get_real(this, section_id, name, value, normalize, norm_object, status)
    !! get scalar real value
    class(input_file),             intent(in)  :: this
    integer,                       intent(in)  :: section_id
      !! section index
    character(*),                  intent(in)  :: name
      !! variable name
    real,                          intent(out) :: value
      !! output value
    logical,             optional, intent(in)  :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in)  :: norm_object
      !! optional normalization object (default: use global normconst)
    logical,             optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i
    logical :: normalize_

    ! optional values
    normalize_ = .true.
    if (present(normalize)) normalize_ = normalize

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%size /= 1) call program_error("variable is not scalar")
            if (var%type /= INPUT_TYPE_REAL) call program_error("variable is not real")
            if (normalize_) then
              value = norm(var%real_data(1), var%unit, n = norm_object)
            else
              value = var%real_data(1)
            end if
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_real_arr(this, section_id, name, values, normalize, norm_object, status)
    !! get array of reals
    class(input_file),             intent(in)  :: this
    integer,                       intent(in)  :: section_id
      !! section index
    character(*),                  intent(in)  :: name
      !! variable name
    real, allocatable,             intent(out) :: values(:)
      !! output values
    logical,             optional, intent(in)  :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in)  :: norm_object
      !! optional normalization object (default: use global normconst)
    logical,             optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i
    logical :: normalize_

    ! optional values
    normalize_ = .true.
    if (present(normalize)) normalize_ = normalize

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%type /= INPUT_TYPE_REAL) call program_error("variable is not real")
            if (normalize_) then
              allocate (values(size(var%real_data)))
              values = norm(var%real_data, var%unit, n = norm_object)
            else
              allocate (values(size(var%real_data)), source = var%real_data)
            end if
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_name_real(this, section_name, name, value, normalize, norm_object, status)
    !! get scalar real value, provide section name instead of index
    class(input_file),             intent(in)  :: this
    character(*),                  intent(in)  :: section_name
      !! section name
    character(*),                  intent(in)  :: name
      !! variable name
    real,                          intent(out) :: value
      !! output value
    logical,             optional, intent(in)  :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in)  :: norm_object
      !! optional normalization object (default: use global normconst)
    logical,             optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, value, normalize = normalize, norm_object = norm_object, status = status)
  end subroutine

  subroutine input_file_get_name_real_arr(this, section_name, name, values, normalize, norm_object, status)
    !! get array of reals, provide section name instead of index

    class(input_file),             intent(in)  :: this
    character(*),                  intent(in)  :: section_name
      !! section name
    character(*),                  intent(in)  :: name
      !! variable name
    real, allocatable,             intent(out) :: values(:)
      !! output values
    logical,             optional, intent(in)  :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in)  :: norm_object
      !! optional normalization object (default: use global normconst)
    logical,             optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, values, normalize = normalize, norm_object = norm_object, status = status)
  end subroutine

  subroutine input_file_get_string(this, section_id, name, value, status)
    !! get scalar string value

    class(input_file),         intent(in)  :: this
    integer,                   intent(in)  :: section_id
      !! section index
    character(*),              intent(in)  :: name
      !! variable name
    character(:), allocatable, intent(out) :: value
      !! output value
    logical,      optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%size /= 1) call program_error("variable is not scalar")
            if (var%type /= INPUT_TYPE_STRING) call program_error("variable is no string")
            value = var%str_data(1)%s
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_string_arr(this, section_id, name, values, status)
    !! get array of strings

    class(input_file),         intent(in)  :: this
    integer,                   intent(in)  :: section_id
      !! section index
    character(*),              intent(in)  :: name
      !! variable name
    type(string), allocatable, intent(out) :: values(:)
      !! output values
    logical,      optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%type /= INPUT_TYPE_STRING) call program_error("variable is no string")
            allocate (values(size(var%str_data)), source = var%str_data)
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_name_string(this, section_name, name, value, status)
    !! get scalar string value, provide section name instead of index
    class(input_file),         intent(in)  :: this
    character(*),              intent(in)  :: section_name
      !! section name
    character(*),              intent(in)  :: name
      !! variable name
    character(:), allocatable, intent(out) :: value
      !! output value
    logical,      optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, value, status = status)
  end subroutine

  subroutine input_file_get_name_string_arr(this, section_name, name, values, status)
    !! get array of strings, provide section name instead of index
    class(input_file),         intent(in)  :: this
    character(*),              intent(in)  :: section_name
      !! section name
    character(*),              intent(in)  :: name
      !! variable name
    type(string), allocatable, intent(out) :: values(:)
      !! output values
    logical,      optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, values, status = status)
  end subroutine

  subroutine input_file_get_logical(this, section_id, name, value, status)
    !! get scalar logical value
    class(input_file), intent(in)  :: this
    integer,           intent(in)  :: section_id
      !! section index
    character(*),      intent(in)  :: name
      !! variable name
    logical,           intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%size /= 1) call program_error("variable is not scalar")
            if (var%type /= INPUT_TYPE_LOGIC) call program_error("variable is not logical")
            value = var%logic_data(1)
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_logical_arr(this, section_id, name, values, status)
    !! get array of reals
    class(input_file),    intent(in)  :: this
    integer,              intent(in)  :: section_id
      !! section index
    character(*),         intent(in)  :: name
      !! variable name
    logical, allocatable, intent(out) :: values(:)
      !! output values
    logical, optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: i

    associate (sect => this%sections%d(section_id))
      do i = 1, sect%vars%n
        associate (var => sect%vars%d(i))
          if (var%name == name) then
            if (var%type /= INPUT_TYPE_LOGIC) call program_error("variable is not logical")
            allocate (values(size(var%logic_data)), source = var%logic_data)
            if (present(status)) status = .true.
            return
          end if
        end associate
      end do
    end associate

    if (present(status)) then
      status = .false.
    else
      print *, name
      call program_error("variable name not found")
    end if
  end subroutine

  subroutine input_file_get_name_logical(this, section_name, name, value, status)
    !! get scalar logical value, provide section name instead of index
    class(input_file), intent(in)  :: this
    character(*),      intent(in)  :: section_name
      !! section name
    character(*),      intent(in)  :: name
      !! variable name
    logical,           intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, value, status = status)
  end subroutine

  subroutine input_file_get_name_logical_arr(this, section_name, name, values, status)
    !! get array of logicals, provide section name instead of index
    class(input_file),    intent(in)  :: this
    character(*),         intent(in)  :: section_name
      !! section name
    character(*),         intent(in)  :: name
      !! variable name
    logical, allocatable, intent(out) :: values(:)
      !! output values
    logical, optional,    intent(out) :: status
      !! optional output if name was found (if not present: error if not found)

    ! local variables
    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      if (st /= 0) then
        status = .false.
        return
      end if
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, values, status = status)
  end subroutine

end module
