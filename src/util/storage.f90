m4_include(macro.f90.inc)

module storage_m
  
  use iso_fortran_env, only: int8, int16, int32, int64, real64

  use error_m
  use normalization_m, only: norm, denorm
  use string_m
  use vector_m

  implicit none
  
  private
  public storage, variable
  public STORAGE_NEW, STORAGE_REPLACE, STORAGE_EDIT, STORAGE_UNPACK
  public DYNAMIC_NO, DYNAMIC_EXT, DYNAMIC_APP
  public ERR_INVALID_ARGUMENTS, ERR_NOT_FOUND, ERR_ALREADY_EXISTS, ERR_FILE_OP, ERR_INTERNAL, ERR_MSGS

  ! Storage general argument flags
  integer, parameter :: STORAGE_UNPACK  = 1
  integer, parameter :: STORAGE_NEW     = 2
  integer, parameter :: STORAGE_EDIT    = 3
  integer, parameter :: STORAGE_REPLACE = 4

  ! Dynamic variable flags
  integer, parameter :: DYNAMIC_NO  = 0 ! No dynamic variable 
  integer, parameter :: DYNAMIC_APP = 1 ! Append to the highest dimension 
  integer, parameter :: DYNAMIC_EXT = 2 ! Dynamic dimension will be on top

  ! Errors
  integer, parameter :: ERR_INVALID_ARGUMENTS = -1
  integer, parameter :: ERR_NOT_FOUND         = -2
  integer, parameter :: ERR_ALREADY_EXISTS    = -3
  integer, parameter :: ERR_FILE_OP           = -4
  integer, parameter :: ERR_INTERNAL          = -5

  character(*), parameter :: ERR_MSGS(5) = [ &
    & "Procedure call with invalid arguments     ", & 
    & "The variable has not been found           ", &
    & "The variable already exists in the storage", &
    & "File operation failed                     ", &
    & "An internal error occured                 "]
  ! Cannot use strings here: Function 'new_string' in initialization expression must be an intrinsic function

  ! Take care when using this macro as commas are not allowed, they will destroy everything (the macro arguments actually)
  m4_define({m4_error}, {
  if (present(stat)) then;
    stat = $1;
    if (present(err_msg)) err_msg = new_string($2); 
    return;
  else;
    call program_error(trim(ERR_MSGS(-$1)) // ": " // $2);
  end if;})

  ! Data types
  integer, parameter :: DT_BIN      = 1
  integer, parameter :: DT_INT32    = 2
  integer, parameter :: DT_INT64    = 3
  integer, parameter :: DT_LOG      = 4
  integer, parameter :: DT_REAL64   = 5
  integer, parameter :: DT_CMPLX128 = 6
  integer, parameter :: DT_CHAR     = 7
  integer, parameter :: DT_STRING   = 8

  character, parameter :: DT_NAMES(8) = ["b", "i", "j", "l", "r", "c", "a", "s"]
  integer,   parameter :: DT_SIZES(7) = [ 1,   4,   8,   1,   8,   16,  1]

  ! All available output types
  m4_define({m4_typelist},{
    m4_X(BIN)
    m4_X(INT32)
    m4_X(INT64)
    m4_X(LOG)
    m4_X(REAL64)
    m4_X(CMPLX128)
    m4_X(STRING)
  })

  ! get expanded typename (e.g. int => integer; my_type => type(my_type))
  m4_define({m4_type},{m4_dnl
  m4_ifelse($1,BIN,integer(kind=int8),{m4_dnl
  m4_ifelse($1,INT32,integer(kind=int32),{m4_dnl
  m4_ifelse($1,INT64,integer(kind=int64),{m4_dnl
  m4_ifelse($1,LOG,logical,{m4_dnl
  m4_ifelse($1,REAL64,real(kind=real64),{m4_dnl
  m4_ifelse($1,CMPLX128,complex(kind=real64),{m4_dnl
  m4_ifelse($1,CHAR,character,{m4_dnl
  type($1)m4_dnl
  })})})})})})})})

  ! dimensions (0..max_dim)
  m4_define({m4_max_dim},{8})
  m4_define({m4_dimlist},{
    m4_ifelse($1,0,,{m4_dimlist(m4_decr($1),$2)})
    m4_Y($1,$2)
  })

  ! array allocate (e.g. m4_shape(4) => :,:,:,:)
  m4_define({m4_allocate},{m4_ifelse($1,1,lbounds(1):ubounds(1),{m4_allocate(m4_decr($1)),lbounds($1):ubounds($1)})})
  m4_define({m4_pallocate},{m4_ifelse($1,0,,{m4_allocate($1)})})

  ! combine type-list with dimension-list to create full list
  m4_define({m4_X},{m4_dimlist(m4_max_dim,$1)})
  m4_define({m4_list},m4_typelist)

  type variable
    !! Corresponds to the journal entry for the variable
    type(string)                 :: name
    integer                      :: type
      !! DT_<type>
    integer(kind=int64), allocatable :: shape(:)
      !! Needed for dynamic arrays as all entries must match in shape
    integer(kind=int64), allocatable :: lbounds(:)
      !! Lower array bounds
    integer(kind=int64)              :: count
      !! Number of blobs
    integer(kind=int64), allocatable :: addr(:)
      !! "Addresses" of the blobs belonging to this variable
  contains
    procedure :: serialize   => variable_serialize
    procedure :: deserialize => variable_deserialize
  end type

  m4_define({T},{string})
  m4_define({U},{variable})
  m4_include(../util/map_def.f90.inc)

  type storage
    !! Interface to Fortran Basic
    integer :: funit = -1
    integer :: mode  = 0
    integer :: wal   = -1
      !! Write-along-log file unit
    
    type(map_string_variable) :: variables
  contains
    procedure :: open     => storage_open 
    procedure :: storage_open
    procedure :: close    => storage_close
    procedure :: contains => storage_contains
    procedure :: goto     => storage_goto

    m4_define({m4_Y},{
    generic :: write => storage_write_$2_$1, storage_writec_$2_$1 
    generic :: read => storage_read_$2_$1, storage_readc_$2_$1
    procedure, private :: storage_write_$2_$1, storage_writec_$2_$1
    procedure, private :: storage_read_$2_$1, storage_readc_$2_$1
    })
    m4_list
  end type

  integer, parameter :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2

  interface
    function fsync (fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int) :: fsync
    end function fsync
  end interface

contains

  function get_dtype(c) result(i)
    character, intent(in) :: c
    integer :: i

    do i = 1, size(DT_names)
      if (DT_names(i) == c) return
    end do
    i = -1
  end function

  subroutine write_header(funit, name, type, lbounds, sizes, unit, compression, data_length, stat, err_msg)
    integer,                       intent(in)  :: funit
    type(string),                  intent(in)  :: name
    integer,                       intent(in)  :: type
    integer(kind=int64),           intent(in)  :: lbounds(:)
    integer(kind=int64),           intent(in)  :: sizes(:)
    character(*),        optional, intent(in)  :: unit
    character,           optional, intent(in)  :: compression
    integer(kind=int64), optional, intent(in)  :: data_length
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg

    character                       :: compr_ = 'n'
    integer                         :: i, name_l, unit_l
    integer(kind=int64)             :: length, data_l
    integer(kind=int8), allocatable :: blob(:)
    if (present(compression)) compr_ = compression

    if (size(lbounds) /= size(sizes)) then
      m4_error(ERR_INVALID_ARGUMENTS, "Sizes and lower bounds do not account for the same number of dimensions")
    end if

    ! The binary blob has the structure
    ! int16:        name length
    ! x * char:     name of the variable
    ! int64:        length of the whole blob
    ! char:         data type
    ! char:         type of compression
    ! int8 + str    unit
    ! int8:         number of dimensions
    ! ndim*2*int64: lower bound + size for each dimension

    name_l = len(name%s) + 2
    if (name_l > 1026) then; m4_error(ERR_INVALID_ARGUMENTS, "Name too long (exceeds 1024 characters)"); end if;
    
    unit_l = 1
    if (present(unit)) then
      unit_l = len(unit) + 1
      if (unit_l > 128) then; m4_error(ERR_INVALID_ARGUMENTS, "Unit too long (exceeds 127 characters)"); end if;
    end if

    length = name_l + 8 + 1 + 1 + unit_l + 1 + 16*size(lbounds)

    if (type == DT_STRING) then
      if (.not. present(data_length)) then; m4_error(ERR_INVALID_ARGUMENTS, "Data length must be provided to write header for string"); end if;
      data_l = data_length
    else
      if (present(data_length)) then 
        data_l = data_length 
      else
        data_l = DT_SIZES(type)
        if (size(lbounds) > 0) data_l = data_l*product(sizes)
      end if
    end if

    allocate(blob(length))

    blob(1:2) = transfer(name_l-2, blob(1), 2)
    do i = 1, int(name_l - 2)
      blob(i+2) = ichar(name%s(i:i), kind=int8)
    end do
    blob(name_l+1: name_l + 8) = transfer(length + data_l, blob(1), 8)
    blob(name_l+9)  = ichar(DT_NAMES(type), kind=int8)
    blob(name_l+10) = ichar(compr_,         kind=int8)
    blob(name_l+11) = int(unit_l-1, kind=int8)
    do i = 1, int(unit_l - 1)
      blob(name_l+11+i) = ichar(unit(i:i), kind=int8)
    end do
    blob(name_l+unit_l+11:name_l+unit_l+11) = transfer(int(size(sizes), kind=int8), blob(1), 1)
    do i = 1, size(sizes)
      blob(name_l+unit_l+12+(i-1)*16 : name_l+unit_l+3 +i*16) = transfer(lbounds(i), blob(1), 8)
      blob(name_l+unit_l+20+(i-1)*16 : name_l+unit_l+11+i*16) = transfer(sizes(i),   blob(1), 8)
    end do

    write(funit) blob
  end subroutine

  subroutine write_log_entry(funit, name, dynamic, dtype, lbounds, sizes)
    integer,             intent(in) :: funit
    type(string),        intent(in) :: name
    integer,             intent(in) :: dynamic
    integer,             intent(in) :: dtype
    integer(kind=int64), intent(in) :: lbounds(:)
    integer(kind=int64), intent(in) :: sizes(:)

    integer      :: i
    character(7) :: type ! "dynamic" or "static"
    if (dynamic == DYNAMIC_NO) then
      type = "dynamic"
    else
      type = "static "
    end if

    write(funit, "(A)", advance="no") '{{"op": "new", "type": "' // trim(type) // '", "name": "' // name%s //'", "dtype": "' // DT_NAMES(dtype) // '", "shape": ['
      do i = 1, size(sizes)
        if (i == 1) then
          write(funit, "(A, I0, A, I0, A)", advance="no") "[",   lbounds(i), ", ", sizes(i), "]"
        else
          write(funit, "(A, I0, A, I0, A)", advance="no") ", [", lbounds(i), ", ", sizes(i), "]"
        end if
      end do
      write(funit, "(A)") ']}}'
  end subroutine

  function get_data_length(var, stat, err_msg) result(length)
    type(string), intent(in) :: var(:)
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg

    integer(kind=int64)      :: i, length

    length = 0
    do i = 1, size(var)
      if (len(var(i)%s) > huge(int(0, kind=int32))) then; m4_error(ERR_INVALID_ARGUMENTS, "String exceeds maximum number of characters"); end if;
      length = length + 4 + len(var(i)%s)
    end do

  end function
  subroutine write_string(funit, var, stat, err_msg)
    integer,      intent(in) :: funit
    type(string), intent(in) :: var(:)
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg

    integer(kind=int64) :: i
    ! Alternatively we could provide custom write and read procedures for string

    do i = 1, size(var)
      if (len(var(i)%s) > 32767) then; m4_error(ERR_INVALID_ARGUMENTS, "String is too long (max 32767 characters)"); end if;
      write(funit) len(var(i)%s, kind=int32)
      write(funit) var(i)%s
    end do

  end subroutine

  subroutine read_header(funit, name, type, unit, lbounds, sizes, compression, data_length)
    integer,                                    intent(in)   :: funit
    type(string),                     optional, intent(out)  :: name
    integer,                          optional, intent(out)  :: type
    character(len=:),    allocatable, optional, intent(out)  :: unit
    integer(kind=int64), allocatable, optional, intent(out)  :: lbounds(:)
    integer(kind=int64), allocatable, optional, intent(out)  :: sizes(:)
    character,                        optional, intent(out)  :: compression
    integer(kind=int64),              optional, intent(out)  :: data_length

    type(string)        :: name_
    integer(kind=int8)  :: i, ndim, type_, compr_, unit_l
    integer(kind=int16) :: name_l
    integer(kind=int64) :: length
    integer(kind=int64), allocatable :: lbounds_(:), sizes_(:)
    character(len=:),    allocatable :: unit_

    length = ftell(funit)

    read (funit) name_l
    allocate(character(name_l) :: name_%s)
    read (funit) name_%s
    if (present(name)) name = name_

    read (funit) length

    read (funit) type_
    if (present(type)) type = get_dtype(achar(type_))    
    read (funit) compr_
    if (present(compression)) compression = achar(compr_)

    read(funit) unit_l
    allocate(character(unit_l) :: unit_)
    read (funit) unit_
    if (present(unit)) call move_alloc(unit_, unit)

    read (funit) ndim
    allocate(lbounds_(ndim), sizes_(ndim))
    do i = 1, ndim
      read (funit) lbounds_(i)
      read (funit) sizes_(i)
    end do
    if (present(lbounds)) call move_alloc(lbounds_, lbounds)
    if (present(sizes))   call move_alloc(sizes_, sizes)

    if (present(data_length)) data_length = length - name_l - 2 - 8 - 1 - 1 - 1 - 16*ndim
  end subroutine

  subroutine read_string(funit, var)
    integer,      intent(in)    :: funit
    type(string), intent(inout) :: var(:)

    integer(kind=int32) :: str_len
    integer(kind=int64) :: i
    ! Alternatively we could provide custom write and read procedures for string

    do i = 1, size(var)
      read (funit) str_len
      allocate(character(str_len) :: var(i)%s)
      read (funit) var(i)%s
    end do
  end subroutine

  subroutine variable_serialize(this, funit)
    class(variable), intent(inout) :: this
    integer,         intent(in)    :: funit
      ! Needs to be opened in binary mode

    call write_header(funit, this%name, this%type, this%lbounds, this%shape, data_length=8+this%count*8)

    write(funit) this%count
    write(funit) this%addr(1:this%count) ! Addresses is not necessarily the size of count
  end subroutine

  subroutine variable_deserialize(this, funit)
    !! Read in the variable blob from the current position of the file descriptor
    class(variable), intent(inout) :: this
    integer,         intent(in)    :: funit

    call read_header(funit, name=this%name, type=this%type, lbounds=this%lbounds, sizes=this%shape)

    read(funit) this%count
    allocate(this%addr(this%count))
    read(funit) this%addr
  end subroutine

  subroutine storage_open(this, file, flag, stat, err_msg)
    !! Initialize the main storgae object which can be used to 
    class(storage), intent(inout) :: this
    character(*),   intent(in)    :: file
    integer,        intent(in)    :: flag
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg
    
    integer         :: stat_
    integer(kind=int64) :: i, addr, n
    character(64)   :: msg
    character(7)    :: journal_start
    
    type(variable), allocatable :: vars(:)


    if (flag == STORAGE_UNPACK) then
      open (newunit=this%funit, file=file, access="stream", form="unformatted", status="old", action="read", iostat=stat_, iomsg=msg)

      call fseek(this%funit, -8, SEEK_END, stat_)
      read(this%funit) addr
      call fseek(this%funit, addr, SEEK_SET, stat_)
      read(this%funit) journal_start

      if ("JOURNAL" /= journal_start) then; m4_error(ERR_INTERNAL, "Invalid journal"); end if;
      read(this%funit) n

      call this%variables%init()
      allocate(vars(n))
      do i = 1, n
        call vars(i)%deserialize(this%funit)
        call this%variables%set(vars(i)%name, vars(i))
      end do

    else if (flag /= STORAGE_NEW) then
      m4_error(ERR_INVALID_ARGUMENTS, "Invalid flag for creating new storages")
    else
      open(newunit=this%funit, file=file, access="stream", form="unformatted", status="new", action="write", iostat=stat_)
      open(newunit=this%wal, file=file//".log", status="new", action="write", iostat=stat_)
    end if
    
    if (stat_ /= 0) then
      m4_error(ERR_FILE_OP, "Could not open storage: " // msg)
    end if

    this%mode = flag

    if (flag == STORAGE_NEW) then
      ! Create the header
      write (this%funit, iostat=stat_) "FBS1"
      call this%variables%init()
    end if

  end subroutine

  subroutine storage_goto(this, name, stat)
    !! Find a variable in the storage and leaves the file position AFTER THE LENGTH OF THE BLOB
    class(storage), intent(inout) :: this
    type(string),   intent(in)    :: name
    integer,        intent(out)   :: stat

    integer(kind=int64) :: pos
    type(variable) :: var
    
    if (.not. this%contains(name)) then
      stat = ERR_NOT_FOUND
      return
    end if

    var = this%variables%get(name)
    if (.not. allocated(var%addr)) then
      stat = ERR_INTERNAL
      return
    end if

    pos = var%addr(1)
    call fseek(this%funit, pos, SEEK_SET, stat)
    if (stat /= 0) stat = ERR_FILE_OP
  end subroutine

  function storage_contains(this, name) result(found)
    class(storage), intent(inout) :: this
    type(string),   intent(in)    :: name
    logical :: found

    type(mapnode_string_variable), pointer :: p => null()

    p => this%variables%find(name)
    found = associated(p)
  end function

  subroutine storage_close(this, delete_log, stat, err_msg)
    class(storage), intent(inout) :: this
    logical, optional, intent(in) :: delete_log
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg
    
    integer :: stat_
    logical :: delete_log_
    integer(kind=int64)              :: i, start, n
    type(variable), allocatable  :: variables(:)        

    delete_log_ = .false.
    if (present(delete_log)) delete_log_ = delete_log

    if (this%mode /= STORAGE_UNPACK) then
      ! Write the journal to the end of the file
      n = this%variables%n
      allocate(variables(n))
      call this%variables%to_array(values=variables)

      call fseek(this%funit, 0, SEEK_END)
      call ftell(this%funit, start)
      write(this%funit) "JOURNAL"
      write(this%funit) n

      do i = 1, n
        call variables(i)%serialize(this%funit)
      end do

      write(this%funit) start

      flush(this%funit)
      stat_ = fsync(fnum(this%funit))
      if (stat_ /= 0) then
        m4_error(ERR_FILE_OP, "Error calling fsync")
      end if

      if (delete_log_) then
        close(this%wal, status="delete")
      else
        close(this%wal)
      end if
    end if

    call this%variables%destruct()
    close(this%funit)
    this%funit = -1
    this%mode = 0
  end subroutine

  m4_define({T},{string})
  m4_define({U},{variable})
  m4_include(../util/map_imp.f90.inc)

  m4_define({m4_Y},{subroutine storage_write_$2_$1(this, name, var, unit, dynamic, stat, err_msg)
    !! 
    class(storage),         intent(inout) :: this
    type(string),           intent(in)    :: name
    m4_type($2),            intent(in)    :: var{}m4_pshape($1)
    character(*), optional, intent(in)    :: unit
    integer,      optional, intent(in)    :: dynamic
    integer,      optional, intent(out)   :: stat
    type(string), optional, intent(out)   :: err_msg

    integer(kind=int64)              :: offset
    integer(kind=int64), allocatable :: sizes(:), lbounds(:), tmp(:)
    integer                          :: stat_
    integer                          :: dflag
    
    type(variable) :: new_var
    type(mapnode_string_variable), pointer :: p => null()

    m4_ifelse($2,STRING,{
    integer(kind=int64)       :: data_l
    type(string), allocatable :: bvar(:)
    },)

    dflag = DYNAMIC_NO
    if (present(dynamic)) dflag = dynamic
    if (dflag < DYNAMIC_NO .or. dflag > DYNAMIC_EXT) then; m4_error(ERR_INVALID_ARGUMENTS, "Unknown dynamic flag"); end if;

    if (this%mode <= STORAGE_UNPACK) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write in this storage") end if;


    if (DT_$2 /= DT_REAL64 .and. DT_$2 /= DT_CMPLX128 .and. present(unit)) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write in this storage") end if;

    call fseek(this%funit, 0, SEEK_END, stat_)
    if (stat_ /= 0) then
      m4_error(ERR_FILE_OP, "Could not go to end of file")
    end if
    call ftell(this%funit, offset)

    sizes = shape(var, kind=int64)
    allocate(lbounds(size(sizes)))
    m4_ifelse($1,0,,{
      lbounds = lbound(var, kind=int64)
    })

    if (size(sizes) /= $1) then
      m4_error(ERR_INVALID_ARGUMENTS, "Shape does not match subroutine definition")
    end if
    
    p => this%variables%find(name)
    if (dflag > DYNAMIC_NO) then
      if (associated(p)) then
        ! Do not allow writing if shape does not match, we simply ignore the lower bounds
                
        if (dflag == DYNAMIC_APP) then
          if (size(sizes) /= size(p%value%shape)) then;                                m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the dimension does not match"); end if; 
          if (any(sizes(1:size(sizes) - 1) /= p%value%shape(1:size(sizes) - 1))) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the shapes do not match"); end if;
          ! Update the lbounds and sizes for writing the binary blob and add the entries to the variable
          lbounds(size(lbounds)) = p%value%lbounds(size(sizes)) - p%value%shape(size(sizes))
          sizes(size(sizes))     = -sizes(size(sizes)) 
          p%value%shape(size(sizes)) = p%value%shape(size(sizes)) + sizes(size(sizes))
        else if (dflag == DYNAMIC_EXT) then
          if (size(sizes) /= size(p%value%shape) - 1) then;     m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the dimension does not match"); end if;
          if (any(sizes /= p%value%shape(1:size(sizes)))) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the shape do not match"); end if;
          
          allocate(tmp(size(p%value%lbounds)))
          tmp(1:size(lbounds)) = lbounds
          tmp(size(tmp))       = p%value%lbounds(size(sizes) + 1) - p%value%shape(size(sizes) + 1)
          call move_alloc(tmp, lbounds)

          allocate(tmp(size(p%value%shape)))
          tmp(1:size(sizes)) = sizes
          tmp(size(tmp))     = int(-1, kind=int64)
          call move_alloc(tmp, sizes)
          
          p%value%shape(size(p%value%shape)) = p%value%shape(size(p%value%shape)) - 1 ! Add one entry to the last dimension
        end if
      else
        new_var%name  = name
        new_var%type  = DT_$2
        new_var%count = int(0,kind=int64)
        new_var%addr  = [offset]
        
        if (dflag == DYNAMIC_APP) then
          sizes(size(sizes)) = -sizes(size(sizes)) ! Indicate that this is a dynamic variable
          new_var%shape              = sizes
          new_var%lbounds            = lbounds
        else if (dflag == DYNAMIC_EXT) then
          allocate(new_var%shape(size(sizes) + 1))
          new_var%shape(1:size(sizes)) = sizes
          new_var%shape(size(sizes) + 1) = int(-1, kind=int64)

          allocate(new_var%lbounds(size(lbounds) + 1))
          new_var%lbounds(1:size(lbounds)) = lbounds
          new_var%lbounds(size(lbounds) + 1) = int(1, kind=int64)
        end if

        lbounds = new_var%lbounds
        sizes   = new_var%shape

        call this%variables%set(name, new_var)
        p => this%variables%find(name)
        call write_log_entry(this%wal, new_var%name, dflag, new_var%type, lbounds, sizes)
      end if

    else if (associated(p)) then
      m4_error(ERR_ALREADY_EXISTS, "The variable already exists")
    else
      new_var%name    = name
      new_var%type    = DT_$2
      new_var%shape   = sizes
      new_var%lbounds = lbounds
      new_var%count   = int(1,kind=int64)
      new_var%addr    = [offset]
      call this%variables%set(name, new_var)
      p => this%variables%find(name) 
    end if

    m4_ifelse($2,STRING,{
      m4_ifelse($1,0,{
        bvar = [var]},{
        bvar = reshape(var, [product(sizes)])
      })
      data_l = get_data_length(bvar)
      call write_header(this%funit, name, DT_$2, lbounds, sizes, data_length=data_l)
      call write_string(this%funit, bvar)  
    },{
      call write_header(this%funit, name, DT_$2, lbounds, sizes)
      m4_ifelse($2,REAL64,{
      if (present(unit)) then
        write (this%funit) denorm(var, unit)
      else 
        write (this%funit) var
      end if},m4_ifelse($2,CMPLX128,{
      if (present(unit)) then
        write (this%funit) denorm(var, unit)
      else 
        write (this%funit) var
      end if
      },{write (this%funit) var}))
    })

    flush(this%funit)
    stat_ = fsync(fnum(this%funit))
    if (stat_ /= 0) then
      m4_error(ERR_FILE_OP, "Error calling fsync")
    end if

    if (all(p%value%shape > 0)) then
      call write_log_entry(this%wal, p%value%name, dflag, p%value%type, p%value%lbounds, p%value%shape)
      return
    end if
    ! Must be a dynamic variable now        
    
    p%value%count = p%value%count + 1
    write(this%wal, "(3A, I0, A)") '{{"op": "add", "name": "', p%value%name%s, '", "idx": ', p%value%count, '}}'

    ! Resize the arrays of the variable structure
    if (p%value%count > size(p%value%addr)) then
      block
        integer(kind=int64), allocatable :: tmp_addr(:)
        allocate (tmp_addr(size(p%value%addr)*2))  
        tmp_addr(1:size(p%value%addr)) = p%value%addr
        call move_alloc(tmp_addr, p%value%addr)
      end block
    end if

    p%value%addr(p%value%count) = offset

    flush(this%wal)
  end subroutine

  subroutine storage_writec_$2_$1(this, name, var, unit, dynamic, stat, err_msg)
    class(storage),         intent(inout) :: this
    character(*),           intent(in)    :: name
    m4_type($2),            intent(in)    :: var{}m4_pshape($1)
    character(*), optional, intent(in)    :: unit
    integer,      optional, intent(in)    :: dynamic
    integer,      optional, intent(out)   :: stat
    type(string), optional, intent(out)   :: err_msg

    call this%write(new_string(name), var, unit, dynamic, stat, err_msg)
  end subroutine})
  m4_list

  m4_define({m4_Y},{subroutine storage_read_$2_$1(this, name, var, stat, err_msg)
    class(storage), intent(inout) :: this
    type(string),   intent(in)    :: name
    m4_type($2) m4_ifelse($1,0,,{, allocatable}), intent(out) :: var{}m4_pshape($1)
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg

    integer             :: stat_
    integer             :: dtype    
    integer(kind=int64) :: sizes_($1)
    integer(kind=int8)  :: ndim
    m4_ifelse($1,0,,{integer(kind=int64) :: i})
    character(len=:),    allocatable :: unit
    integer(kind=int64), allocatable :: sizes(:), lbounds(:), ubounds(:)
    m4_ifelse($2,STRING,{
      type(string), allocatable :: bvar(:)
      ! Better solution would be to give bvar the pointer attribute and use array bounds remapping
      ! to avoid additional memory usage of reshape. However, testing this resulted in the compiler
      ! complaining about 'var' not being SIMPLY CONTIGUOUS, but it should be. 
    },)

    type(mapnode_string_variable), pointer :: p => null()
    
    ! Implies that the variable is registered
    call this%goto(name, stat_)
    if (stat_ /= 0) then
      m4_error(ERR_NOT_FOUND, "Could not find the variable " // name%s)
    end if
    p => this%variables%find(name)

    call read_header(this%funit, type=dtype, unit=unit, lbounds=lbounds, sizes=sizes)
    ubounds = lbounds + sizes - 1
    if (dtype /= DT_$2) then
      m4_error(ERR_INVALID_ARGUMENTS, "The variable type '" // DT_NAMES(dtype) // "' does not match the called procedure ($2)")
    end if
    ndim = size(sizes, kind=int8)
    if (ndim /= $1) then
      m4_error(ERR_INVALID_ARGUMENTS, "Variable declaration does not match the stored variable")
    end if
    sizes_ = sizes

    m4_ifelse($1,0,{
    m4_ifelse($2,STRING,{
      bvar = [var]
      call read_string(this%funit, bvar)
      var = bvar(1) 
    },{
      read (this%funit) var

      m4_ifelse($2,REAL64,{
      if (len(unit) > 0) then
        var = norm(var, unit)
      end if},m4_ifelse($2,CMPLX128,{
      if (len(unit) > 0) then
        var = norm(var, unit)
      end if}))
      
    })
    },{
    if (any(sizes <= 0)) then ! Dynamic array
      if (all(p%value%shape > 0)) then; m4_error(ERR_INTERNAL, "Dynamic disagree between blob and journal"); end if;

      ! Allocate the array with the complete shape
      sizes(size(sizes)) = abs(p%value%shape(size(p%value%shape)))
      lbounds(size(lbounds)) = abs(p%value%lbounds(size(p%value%lbounds)))
      ubounds = lbounds + sizes - 1
      allocate(var({}m4_pallocate($1)))

      do i = 1, p%value%count
        ! TODO: Check validity of each binary blob
        call fseek(this%funit, p%value%addr(i), SEEK_SET)
        deallocate(lbounds, sizes, ubounds)
        call read_header(this%funit, lbounds=lbounds, sizes=sizes)
        sizes(size(sizes)) = abs(sizes(size(sizes)))
        ubounds = lbounds + sizes - 1

        m4_ifelse($2,STRING,{
        bvar = reshape(var({}m4_pallocate($1)), [product(sizes)])
        call read_string(this%funit, bvar)
        var({}m4_pallocate($1)) = reshape(bvar, sizes_)
        },{
        read (this%funit) var({}m4_pallocate($1))
        m4_ifelse($2,REAL64,{
        if (len(unit) > 0) then
          var = norm(var, unit)
        end if},m4_ifelse($2,CMPLX128,{
        if (len(unit) > 0) then
          var = norm(var, unit)
        end if}))
        })
      end do
    else ! Static array
      allocate(var({}m4_pallocate($1)))
      m4_ifelse($2,STRING,{
        bvar = reshape(var, [product(sizes)])
        call read_string(this%funit, bvar)  
        var = reshape(bvar, sizes_)
      },{
        read(this%funit) var
        m4_ifelse($2,REAL64,{
        if (len(unit) > 0) then
          var = norm(var, unit)
        end if},m4_ifelse($2,CMPLX128,{
        if (len(unit) > 0) then
          var = norm(var, unit)
        end if}))
      })
    end if
    })
    
  end subroutine
  subroutine storage_readc_$2_$1(this, name, var, stat, err_msg)
    class(storage), intent(inout) :: this
    character(*),   intent(in)    :: name
    m4_type($2) m4_ifelse($1,0,,{, allocatable}), intent(out) :: var{}m4_pshape($1)
    integer,             optional, intent(out) :: stat
    type(string),        optional, intent(out) :: err_msg

    call this%read(new_string(name), var, stat, err_msg)
  end subroutine})
  m4_list
end module