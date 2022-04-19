m4_include(macro.f90.inc)

module input_m

  use error_m,         only: program_error
  use iso_fortran_env, only: int32, int64
  use normalization_m, only: norm, normalization
  use vector_m,        only: vector_int, vector_real, vector_string
  use string_m,        only: string

  implicit none

  ! list of supported data-types
  m4_define({m4_list},{
    m4_X(int32)
    m4_X(int64)
    m4_X(log)
    m4_X(real)
    m4_X(string)
  })

  ! get input type based on typename
  m4_define({m4_get_input_type},{m4_dnl
  m4_ifelse($1,int32,INPUT_TYPE_INTEGER,{m4_dnl
  m4_ifelse($1,int64,INPUT_TYPE_INTEGER,{m4_dnl
  m4_ifelse($1,log,INPUT_TYPE_LOGIC,{m4_dnl
  m4_ifelse($1,real,INPUT_TYPE_REAL,{m4_dnl
  m4_ifelse($1,string,INPUT_TYPE_STRING,)})})})})})

  ! get input_var data variable based on typename
  m4_define({m4_get_data},{m4_dnl
  m4_ifelse($1,int32,int_data,{m4_dnl
  m4_ifelse($1,int64,int_data,{m4_dnl
  m4_ifelse($1,log,logic_data,{m4_dnl
  m4_ifelse($1,real,real_data,{m4_dnl
  m4_ifelse($1,string,str_data,)})})})})})

  private
  public input_file

  ! maximum line length
  integer, parameter :: MAX_LINE_LEN = 4096

  ! variable types
  integer, parameter :: INPUT_TYPE_INTEGER = 1
  integer, parameter :: INPUT_TYPE_REAL    = 2
  integer, parameter :: INPUT_TYPE_STRING  = 3
  integer, parameter :: INPUT_TYPE_LOGIC   = 4
  integer, parameter :: INPUT_NUM_TYPES    = 4

  ! variable type names
  character(7), parameter :: INPUT_TYPE_NAMES(4) = [ &
    "integer", &
    "real   ", &
    "string ", &
    "logical"  &
  ]

  type input_var
    character(:),   allocatable :: name
      !! Variable name
    integer                     :: type
      !! Variable type
    integer                     :: size
      !! Size of variable, 1 means scalar, > 1 means array
    integer(int64), allocatable :: int_data(:)
      !! integer data
    real,           allocatable :: real_data(:)
      !! real data
    type(string),   allocatable :: str_data(:)
      !! string data
    logical,        allocatable :: logic_data(:)
      !! logical data
    character(:),   allocatable :: unit
      !! Physical unit token (for real variables)
  contains
    procedure :: init => input_var_init
  end type

  m4_define({T},{input_var})
  m4_include(vector_def.f90.inc)

  type input_section
    character(:), allocatable :: name
      !! Section name
    type(vector_input_var)    :: vars
      !! variables in this section
  contains
    procedure :: init    => input_section_init
    procedure :: add_var => input_section_add_var
  end type

  m4_define({T},{input_section})
  m4_include(vector_def.f90.inc)

  type input_file
    character(:), allocatable  :: name
      !! File name
    character(:), allocatable  :: default
      !! Default file name
    character(:), allocatable  :: directory
      !! Directory
    type(vector_input_section) :: sections
      !! Sections in this file
    type(vector_input_section) :: default_sections
      !! Default values for sections
  contains
    private
    procedure :: load       => input_file_load
    procedure :: parse_line => input_file_parse_line
    procedure :: load_csv   => input_file_load_csv

    procedure :: get_var    => input_file_get_var

    m4_define({m4_X},{
      procedure :: input_file_get_$1
      procedure :: input_file_get_$1_n
      procedure :: input_file_get_name_$1
      procedure :: input_file_get_name_$1_n
    })
    m4_list

    procedure, public :: init         => input_file_init
    procedure, public :: get_section  => input_file_get_section
    procedure, public :: get_sections => input_file_get_sections

    m4_define({m4_X},{
      generic, public :: get => input_file_get_$1
      generic, public :: get => input_file_get_$1_n
      generic, public :: get => input_file_get_name_$1
      generic, public :: get => input_file_get_name_$1_n
    })
    m4_list
  end type

contains

  m4_define({T},{input_var})
  m4_include(vector_imp.f90.inc)

  m4_define({T},{input_section})
  m4_include(vector_imp.f90.inc)

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

    integer :: i

    ! save names
    this%name    = name
    this%default = ""
    if (present(default)) this%default = default

    ! get directory
    i = scan(name, "/", back = .true.)
    if (i > 1) then
      this%directory = name(1:i)
    else
      this%directory = ""
    end if

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
        print "(A,I0)", "Input error in line number ", line_number
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

    m4_ignore(this)

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

    ! macro (only #load_csv(example.csv) supported right now)
    if (line(1:1) == '#') then
      ! remove # symbol
      if (line_len < 2) then
        valid = .false.
        err = "macro name missing"
        return
      end if
      line = adjustl(line(2:))
      line_len = len_trim(line)

      ! get macro name
      i = scan(line, '(')
      if (i == 0) i = len_trim(line)
      if (i < 2) then
        valid = .false.
        err = "macro name error"
        return
      end if
      allocate (character(0) :: name)      ! remove gfortran warning
      name = line(1:i-1)

      if (name == "load_csv") then
        if ((line(line_len:line_len) /= ')') .or. (i > line_len-2)) then
          valid = .false.
          err = "load_csv macro error"
          return
        end if

        ! get file name
        name = trim(adjustl(line(i+1:line_len-1)))
        call this%load_csv(sections, this%directory//name, valid, err)
      else
        valid = .false.
        err = "unknown macro "//name
        return
      end if

      return
    end if

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

  subroutine input_file_load_csv(this, sections, filename, valid, err)
    !! load csv file and insert variables into current section
    class(input_file),                    intent(inout) :: this
    type(vector_input_section),           intent(inout) :: sections
      !! Input sections
    character(*),                         intent(in)    :: filename
      !! csv filename
    logical,                              intent(out)   :: valid
      !! Output if parsing was successful (valid syntax)
    character(:), allocatable,            intent(out)   :: err
      !! Error string

    integer                        :: funit, line_number, status, nvar, i, j
    character(10)                  :: lfmt
    character(32)                  :: snum
    character(:),      allocatable :: tmp
    character(MAX_LINE_LEN)        :: line
    real                           :: val
    type(vector_int)               :: sep
    type(string),      allocatable :: names(:), units(:)
    type(vector_real), allocatable :: data(:)
    type(input_var)                :: v

    m4_ignore(this)

    ! assume valid syntax; change later if invalid
    valid = .true.
    err   = ""

    ! line format
    write(lfmt, "(A,I0,A)") "(A", MAX_LINE_LEN, ")"

    ! open file
    open(newunit = funit, file = trim(filename), status = "old", action = "read")

    ! read first line (header)
    read(funit, fmt = lfmt, iostat = status) line
    if (status < 0) then
      close(funit)
      valid = .false.
      err = "load_csv("//filename//"): No header line"
      return
    elseif (status > 0) then
      ! IO error
      close(funit)
      valid = .false.
      err = "load_csv("//filename//"): IO Error"
      return
    end if

    ! separate line at each ',' character
    call sep%init(0, c = 8)
    call sep%push(0)
    do i = 1, len_trim(line)
      if (line(i:i) == ',') call sep%push(i)
    end do
    call sep%push(len_trim(line)+1)

    ! number of variables
    nvar = sep%n-1

    ! allocate
    allocate (names(nvar))
    allocate (units(nvar))
    allocate (data(nvar))

    ! get names and units, init data vectors
    do i = 1, nvar
      if (sep%d(i)+1 > sep%d(i+1)-1) then
        close(funit)
        valid = .false.
        err = "load_csv("//filename//"): Error in header line: Name empty"
        return
      end if
      tmp = adjustl(trim(line(sep%d(i)+1:sep%d(i+1)-1)))
      if (len_trim(tmp) < 1) then
        close(funit)
        valid = .false.
        err = "load_csv("//filename//"): Error in header line: Name empty"
        return
      end if

      ! split at ':'
      j = scan(tmp, ':')
      if (j == 0) then
        names(i)%s = tmp
        units(i)%s = "1"
      else
        if (j == 1) then
          close(funit)
          valid = .false.
          err = "load_csv("//filename//"): Error in header line: Name empty"
          return
        end if
        if (j == len_trim(tmp)) then
          close(funit)
          valid = .false.
          err = "load_csv("//filename//"): Error in header line: Unit empty"
          return
        end if
        names(i)%s = adjustl(trim(tmp(1:j-1)))
        units(i)%s = adjustl(trim(tmp(j+1:len_trim(tmp))))
      end if

      call data(i)%init(0, c = 256)
    end do

    ! first data line
    line_number = 1

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
        valid = .false.
        err = "load_csv("//filename//"): IO Error"
        return
      end if
      line_number = line_number + 1
      write (snum, "(I0)") line_number

      ! separate line at ',' characters
      call sep%resize(0)
      call sep%push(0)
      do i = 1, len_trim(line)
        if (line(i:i) == ',') call sep%push(i)
      end do
      call sep%push(len_trim(line)+1)
      if (sep%n < nvar + 1) then
        close(funit)
        valid = .false.
        err = "load_csv("//filename//"): Error in line "//trim(snum)//": Not enough values separated by ,"
        return
      end if
      if (sep%n > nvar + 1) then
        close(funit)
        valid = .false.
        err = "load_csv("//filename//"): Error in line "//trim(snum)//": Too many values"
        return
      end if

      do i = 1, nvar
        if (sep%d(i)+1 > sep%d(i+1)-1) then
          close(funit)
          valid = .false.
          err = "load_csv("//filename//"): Error in line "//trim(snum)//": Value empty"
          return
        end if
        tmp = adjustl(trim(line(sep%d(i)+1:sep%d(i+1)-1)))
        if (len_trim(tmp) < 1) then
          close(funit)
          valid = .false.
          err = "load_csv("//filename//"): Error in line "//trim(snum)//": Value empty"
          return
        end if

        ! convert to real
        read (tmp, *, iostat = status) val
        if (status /= 0) then
          close(funit)
          valid = .false.
          err   = "load_csv("//filename//"): Error in line "//trim(snum)//": Invalid real data"
          return
        end if
        call data(i)%push(val)
      end do
    end do

    if (data(1)%n < 1) then
      valid = .false.
      err   = "load_csv("//filename//"): No data"
      return
    end if

    ! if no sections, add one without name
    if (sections%n == 0) then
      call sections%resize(1)
      allocate (character(0) :: tmp)
      call sections%d(1)%init(tmp)
    end if

    ! add variables to current section
    v%size = data(1)%n
    allocate (v%real_data(data(1)%n))
    do i = 1, nvar
      v%name      = names(i)%s
      v%type      = INPUT_TYPE_REAL
      v%real_data = data(i)%d(1:data(i)%n)
      v%unit      = units(i)%s
      call sections%d(sections%n)%add_var(v)
    end do
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
    integer :: i, status_

    status_ = -1
    do i = 1, this%sections%n
      if (this%sections%d(i)%name /= section_name) cycle
      section_id = i
      status_    = status_ + 1
    end do

    if (present(status)) then
      status = status_
    else
      if (status_ < 0) call program_error("section name '"//section_name//"' not found")
      if (status_ > 0) call program_error("found multiple sections with name '"//section_name//"', use get_sections instead")
    end if
  end subroutine

  subroutine input_file_get_sections(this, section_name, section_ids)
    !! get all section indices with this name

    class(input_file),    intent(in)  :: this
    character(*),         intent(in)  :: section_name
      !! name of section to search
    integer, allocatable, intent(out) :: section_ids(:)
      !! output indices

    integer          :: i
    type(vector_int) :: sv

    ! collect sections with this name
    call sv%init(0, c = this%sections%n)
    do i = 1, this%sections%n
      if (this%sections%d(i)%name == section_name) call sv%push(i)
    end do

    ! output section ids
    allocate (section_ids(sv%n), source = sv%d(1:sv%n))
  end subroutine

  subroutine input_file_get_var(this, section_id, name, type, scalar, var_idx, status)
    !! search for variable by name
    class(input_file)              :: this
    integer,           intent(in)  :: section_id
      !! section index
    character(*),      intent(in)  :: name
      !! variable name
    integer,           intent(in)  :: type
      !! variable type
    logical,           intent(in)  :: scalar
      !! scalar or array
    integer,           intent(out) :: var_idx
      !! output variable index
    logical, optional, intent(out) :: status
      !! optional output if appropriate var was found (if not present: error if not found)

    if (present(status)) status = .false.

    associate (sec => this%sections%d(section_id))
      do var_idx = 1, sec%vars%n
        associate (var => sec%vars%d(var_idx))
          if (var%name /= name) cycle
          if (scalar .and. (var%size /= 1)) then
            if (present(status)) return
            call program_error("variable '"//name//"' is not scalar")
          end if
          if (var%type /= type) then
            if (present(status)) return
            call program_error("variable '"//name//"' is not of type "//trim(INPUT_TYPE_NAMES(type)))
          end if
          if (present(status)) status = .true.
          return
        end associate
      end do
    end associate

    if (.not. present(status)) call program_error("variable '"//name//"' not found")
  end subroutine

  m4_define({m4_X},{subroutine input_file_get_$1(this, section_id, name, value, status{}m4_ifelse($1,real,{, normalize, norm_object},))
    !! get scalar value
    class(input_file), intent(in)  :: this
    integer,           intent(in)  :: section_id
      !! section index
    character(*),      intent(in)  :: name
      !! variable name
    m4_type($1),       intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)
    m4_ifelse($1,real,{
    logical,             optional, intent(in) :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in) :: norm_object
      !! optional normalization object (default: use global normconst)
    },)

    integer :: var_idx
    m4_ifelse($1,real,{logical :: normalize_},)

    call this%get_var(section_id, name, m4_get_input_type($1), .true., var_idx, status = status)
    if (present(status)) then
      if (.not. status) return
    end if

    associate (var => this%sections%d(section_id)%vars%d(var_idx))
      m4_ifelse($1,int32,{
        value = int(var%int_data(1), kind = int32)
      },{
        value = var%{}m4_get_data($1){}(1)
      })

      m4_ifelse($1,real,{
        normalize_ = .true.
        if (present(normalize)) normalize_ = normalize
        if (normalize_) value = norm(value, var%unit, n = norm_object)
      })
    end associate
  end subroutine})
  m4_list

  m4_define({m4_X},{subroutine input_file_get_$1_n(this, section_id, name, values, status{}m4_ifelse($1,real,{, normalize, norm_object},))
    !! get array of values
    class(input_file),             intent(in)  :: this
    integer,                       intent(in)  :: section_id
      !! section index
    character(*),                  intent(in)  :: name
      !! variable name
    m4_type($1), allocatable,      intent(out) :: values(:)
      !! output values
    logical,             optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)
    m4_ifelse($1,real,{
    logical,             optional, intent(in) :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in) :: norm_object
      !! optional normalization object (default: use global normconst)
    },)

    integer :: var_idx
    m4_ifelse($1,real,{logical :: normalize_},)

    call this%get_var(section_id, name, m4_get_input_type($1), .false., var_idx, status = status)
    if (present(status)) then
      if (.not. status) return
    end if

    associate (var => this%sections%d(section_id)%vars%d(var_idx))
      m4_ifelse($1,int32,{
        allocate (values(size(var%int_data)), source = int(var%int_data, kind = int32))
      },{
        values = var%{}m4_get_data($1)
      })

      m4_ifelse($1,real,{
        normalize_ = .true.
        if (present(normalize)) normalize_ = normalize
        if (normalize_) values = norm(values, var%unit, n = norm_object)
      })
    end associate
  end subroutine})
  m4_list

  m4_define({m4_X},{subroutine input_file_get_name_$1(this, section_name, name, value, status{}m4_ifelse($1,real,{, normalize, norm_object},))
    !! get scalar value, provide section name instead of index
    class(input_file), intent(in)  :: this
    character(*),      intent(in)  :: section_name
      !! section name
    character(*),      intent(in)  :: name
      !! variable name
    m4_type($1),       intent(out) :: value
      !! output value
    logical, optional, intent(out) :: status
      !! optional output if name was found (if not present: error if not found)
    m4_ifelse($1,real,{
    logical,             optional, intent(in) :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in) :: norm_object
      !! optional normalization object (default: use global normconst)
    },)

    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      status = (st == 0)
      if (.not. status) return
    else
      call this%get_section(section_name, section_id)
    end if

    ! get value
    call this%get(section_id, name, value, status = status{}m4_ifelse($1,real,{, normalize = normalize, norm_object = norm_object},))
  end subroutine})
  m4_list

  m4_define({m4_X},{subroutine input_file_get_name_$1_n(this, section_name, name, values, status{}m4_ifelse($1,real,{, normalize, norm_object},))
    !! get array of values, provide section name instead of index
    class(input_file),        intent(in)  :: this
    character(*),             intent(in)  :: section_name
      !! section name
    character(*),             intent(in)  :: name
      !! variable name
    m4_type($1), allocatable, intent(out) :: values(:)
      !! output values
    logical, optional,           intent(out) :: status
      !! optional output if name was found (if not present: error if not found)
    m4_ifelse($1,real,{
    logical,             optional, intent(in) :: normalize
      !! normalize value by using unit token (default: true)
    type(normalization), optional, intent(in) :: norm_object
      !! optional normalization object (default: use global normconst)
    },)

    integer :: section_id, st

    ! get section
    if (present(status)) then
      call this%get_section(section_name, section_id, status = st)
      status = (st == 0)
      if (.not. status) return
    else
      call this%get_section(section_name, section_id)
    end if

    ! get values
    call this%get(section_id, name, values, status = status{}m4_ifelse($1,real,{, normalize = normalize, norm_object = norm_object},))
  end subroutine})
  m4_list

end module
