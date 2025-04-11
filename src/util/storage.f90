m4_include(macro.f90.inc)

module storage_m
  
  use iso_fortran_env, only: int8, int16, int32, int64, real64, logical_kinds
  use, intrinsic :: iso_c_binding

  use error_m
  use normalization_m, only: norm, denorm
  use string_m
  use vector_m
  use util_m, only: int2str
  m4_ifdef({m4_zlib},use zlib_m)
  m4_ifdef({m4_blosc},use blosc_m)

  implicit none

  private
  public storage, variable

  ! Storage general argument flags
  integer, parameter, public :: STORAGE_READ  = 1
    !! Read-only flag for storage usage
  integer, parameter, public :: STORAGE_WRITE = 2
    !! Use the storage object for write operations.

  ! Dynamic variable flags
  integer, parameter, public :: DYNAMIC_NO  = 0 ! No dynamic variable 
  integer, parameter, public :: DYNAMIC_APP = 1 ! Append to the highest dimension 
  integer, parameter, public :: DYNAMIC_EXT = 2 ! Dynamic dimension will be on top

  ! Compression flags
  character, parameter, public :: COMPR_DEFAULT = m4_ifdef({m4_blosc},'b',m4_ifdef({m4_zlib}, 'z','n'))
  character, parameter, public :: COMPR_NONE  = 'n' ! No 
  character, parameter, public :: COMPR_ZLIB  = 'z' ! Zlib deflate best compression
  character, parameter, public :: COMPR_BLOSC = 'b' ! Blosc compression 9 with filter

  ! Errors
  integer, parameter, public :: ERR_INVALID_ARGUMENTS = -1
  integer, parameter, public :: ERR_NOT_FOUND         = -2
  integer, parameter, public :: ERR_ALREADY_EXISTS    = -3
  integer, parameter, public :: ERR_FILE_OP           = -4
  integer, parameter, public :: ERR_INTERNAL          = -5

  character(*), parameter, public :: ERR_MSGS(5) = [ &
    & "Invalid arguments    ", & 
    & "Not found            ", &
    & "Already exists       ", &
    & "File operation failed", &
    & "Internal error       "]
  ! Cannot use strings here: Function 'new_string' in initialization expression must be an intrinsic function

  ! Take care when using this macro as commas are not allowed, they will destroy everything (the macro arguments actually)
  m4_define({m4_error}, {
  if (present(stat)) then;
    stat = $1;
    if (present(err_msg)) err_msg = new_string($2); 
    return;
  else;
    call program_error($2);
  end if;})

  ! Data types
  integer, parameter :: DT_BIN      = 1
  integer, parameter :: DT_INT32    = 2
  integer, parameter :: DT_INT64    = 3
  integer, parameter :: DT_LOG32    = 4
  integer, parameter :: DT_LOG64    = 5
  integer, parameter :: DT_REAL64   = 6
  integer, parameter :: DT_CMPLX128 = 7
  integer, parameter :: DT_CHAR     = 8
  integer, parameter :: DT_STRING   = 9

  character, parameter :: DT_NAMES(9) = ["b", "i", "j", "l", "m", "r", "c", "a", "s"]
  integer,   parameter :: DT_SIZES(8) = [ 1,   4,   8,   4,   8,   8,   16,  1]

  ! All available output types
  m4_define({m4_typelist},{
    m4_X(BIN)
    m4_X(INT32)
    m4_X(INT64)
    m4_X(LOG32)
    m4_X(LOG64)
    m4_X(REAL64)
    m4_X(CMPLX128)
    m4_X(STRING)
  })

  ! get expanded typename (e.g. int => integer; my_type => type(my_type))
  m4_define({m4_type},{m4_dnl
  m4_ifelse($1,BIN,integer(kind=int8),{m4_dnl
  m4_ifelse($1,INT32,integer(kind=int32),{m4_dnl
  m4_ifelse($1,INT64,integer(kind=int64),{m4_dnl
  m4_ifelse($1,LOG32,logical(kind=4),{m4_dnl
  m4_ifelse($1,LOG64,logical(kind=8),{m4_dnl
  m4_ifelse($1,REAL64,real(kind=real64),{m4_dnl
  m4_ifelse($1,CMPLX128,complex(kind=real64),{m4_dnl
  m4_ifelse($1,CHAR,character,{m4_dnl
  type($1)m4_dnl
  })})})})})})})})})

  m4_define({m4_denormable},{m4_dnl
    m4_ifelse($1,REAL64,$2,{m4_ifelse($1,CMPLX128,$2)})m4_dnl
  })

  ! dimensions (0..max_dim)
  m4_define({m4_max_dim},{8})
  m4_define({m4_dimlist},{
    m4_ifelse($1,0,,{m4_dimlist(m4_decr($1),$2)})
    m4_Y($1,$2)
  })

  ! array allocate (e.g. m4_pallocate(2) => lbounds(1):ubounds(1),lbounds(2):ubounds(2)
  m4_define({m4_allocate},{m4_ifelse($1,1,lbounds(1):ubounds(1),{m4_allocate(m4_decr($1)),lbounds($1):ubounds($1)})})
  m4_define({m4_pallocate},{m4_ifelse($1,0,,{m4_allocate($1)})})

  m4_define({m4_dimm1},{m4_ifelse($1,1,,{m4_dimm1(m4_decr($1)):,})})
  m4_define({m4_pdimm1},{m4_ifelse($1,0,,{m4_dimm1($1)})})

  ! combine type-list with dimension-list to create full list
  m4_define({m4_X},{m4_dimlist(m4_max_dim,$1)})
  m4_define({m4_list},m4_typelist)

  type variable
    !! Corresponds to the journal entry for the variable
    type(string)                 :: name
    integer                      :: type
      !! DT_<type>
    integer(kind=int64), allocatable :: sizes(:)
      !! Needed for dynamic arrays as all entries must match in shape
    integer(kind=int64), allocatable :: lbounds(:)
      !! Lower array bounds
    integer(kind=int64)              :: count
      !! Number of blobs
    integer(kind=int64), allocatable :: span(:)
    integer(kind=int64), allocatable :: addr(:)
      !! "Addresses" of the blobs belonging to this variable
  contains
    procedure :: serialize   => variable_serialize
    procedure :: deserialize => variable_deserialize
    generic   :: write(formatted) => variable_write_formatted ! Formatted output corresponds to log entries

    procedure, private :: variable_write_formatted
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

  integer, parameter :: SEEK_SET = int(0, kind=4), SEEK_CUR = int(1, kind=4), SEEK_END = int(2, kind=4)

  interface
    function fsync (fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int) :: fsync
    end function fsync
  end interface

contains

  subroutine init_error(stat, err_msg)
    integer(kind=4), optional, intent(inout) :: stat
    type(string),    optional, intent(inout) :: err_msg

    if (present(stat)) then
      stat = 0
      if (present(err_msg)) err_msg = new_string("")
    end if

  end subroutine

  function get_dtype(c) result(i)
    character, intent(in) :: c
    integer :: i

    do i = 1, size(DT_names)
      if (DT_names(i) == c) return
    end do
    i = -1
  end function

  function write_header(funit, name, type, lbounds, sizes, unit, compression, data_length, stat, err_msg) result(length)
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
  end function

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

    if (present(data_length)) data_length = length - 2 - name_l - 8 - 1 - 1 - 1 - unit_l - 1 - 16*ndim

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

    integer(kind=int64) :: i, header_l

    header_l =  write_header(funit, this%name, this%type, this%lbounds, this%sizes, data_length=8+this%count*8)

    write(funit) this%count
    do i = 1, this%count ! Addresses and indices are not necessarily the size of count
      write(funit) this%addr(i) 
      write(funit) this%span(i)
    end do
  end subroutine

  subroutine variable_deserialize(this, funit)
    !! Read in the variable blob from the current position of the file descriptor
    class(variable), intent(inout) :: this
    integer,         intent(in)    :: funit

    integer(kind=int64) :: i

    call read_header(funit, name=this%name, type=this%type, lbounds=this%lbounds, sizes=this%sizes)

    read(funit) this%count
    allocate(this%addr(this%count), this%span(this%count))
    
    do i = 1, this%count
      read(funit) this%addr(i)
      read(funit) this%span(i)
    end do
  end subroutine

  subroutine variable_write_formatted(this, unit, iotype, v_list, iostat, iomsg)
    CLASS(variable), intent(in)    :: this
    integer,         intent(in)    :: unit
    character(*),    intent(in)    :: iotype
    integer,         intent(in)    :: v_list(:)
    integer,         intent(out)   :: iostat
    character(*),    intent(inout) :: iomsg

    integer :: i

    m4_ignore(iotype)
    m4_ignore(v_list)

    write(unit, "(A)", advance="no", iostat=iostat, iomsg=iomsg) '{{"op": "...", "type": "...", "name": "' // this%name%s //'", "dtype": "' // DT_NAMES(this%type) // '", "shape": ['
      do i = 1, size(this%sizes)
        if (i == 1) then
          write(unit, "(A, I0, A, I0, A)", advance="no", iostat=iostat, iomsg=iomsg) "[",   this%lbounds(i), ", ", this%sizes(i), "]"
        else
          write(unit, "(A, I0, A, I0, A)", advance="no", iostat=iostat, iomsg=iomsg) ", [", this%lbounds(i), ", ", this%sizes(i), "]"
        end if
      end do
      write(unit, "(A)", iostat=iostat, iomsg=iomsg) ']}}'
  end subroutine

  subroutine storage_open(this, file, flag, stat, err_msg)
    !! Initialize the main storgae object which can be used to 
    class(storage),            intent(inout) :: this
    character(*),              intent(in)    :: file
    integer,         optional, intent(in)    :: flag
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg
    
    integer             :: flag_
    integer(kind=4)     :: stat_
    integer(kind=int64) :: i, addr, n
    character(64)       :: msg
    logical             :: exists
    
    type(variable), allocatable :: vars(:)

    flag_ = STORAGE_READ
    if (present(flag)) flag_ = flag

    call init_error(stat, err_msg)

    select case (flag_)
    case (STORAGE_READ)
      open (newunit=this%funit, file=file, access="stream", form="unformatted", status="old", action="read", iostat=stat_, iomsg=msg)
      if (stat_ /= 0) then; m4_error(ERR_INTERNAL, trim(msg)); return; end if;

      call read_journal(stat_)
      if (stat_ /= 0) then; m4_error(ERR_INTERNAL, "Reading while writing is not supported"); return; end if;

    case (STORAGE_WRITE)
      inquire (file=file, exist=exists)
      
      ! New storage
      if (.not. exists) then
        open(newunit=this%funit, file=file, access="stream", form="unformatted", status="new", action="readwrite", iostat=stat_, iomsg=msg)
        open(newunit=this%wal, file=file//".log", status="new", action="write", iostat=stat_, iomsg=msg)
        
        if (stat_ /= 0) then
          m4_error(ERR_FILE_OP, "Could not open storage: " // msg)
        end if
  
        write (this%funit, iostat=stat_) "FBS1"
        call this%variables%init()
      
      ! Extend old storage
      else
        open (newunit=this%funit, file=file, access="stream", form="unformatted", status="old", action="readwrite", iostat=stat_, iomsg=msg)
        if (stat_ /= 0) then; m4_error(ERR_FILE_OP, trim(msg)); return; end if;

        inquire (file=file//".log", exist=exists)
        if (exists) then; m4_error(ERR_ALREADY_EXISTS, "The storage is currently openend for writing and no more than one writer allowed"); return; end if;
        open(newunit=this%wal, file=file//".log", status="new", action="write", iostat=stat_, iomsg=msg)
        if (stat_ /= 0) then; m4_error(ERR_FILE_OP, trim(msg)); return; end if;
        
        call read_journal(stat_)
        if (stat_ /= 0) then
          m4_error(ERR_FILE_OP, "Could not read the journal")
        end if

      end if
      
    case default
      m4_error(ERR_INVALID_ARGUMENTS, "Invalid flag for creating new storages")
    end select

    this%mode = flag_

    ! block
    !   integer :: m
    !   type(variable), allocatable :: vars(:)
    !   allocate(vars(this%variables%n))
    !   call this%variables%to_array(values=vars)
    !   write(*,*)
    !   write(*,*) "Current state when opening"
    !   do m = 1, this%variables%n
    !     write(*,*) "    ", vars(m)%name%s, "(", vars(m)%type, "): ", vars(m)%sizes, " | ", vars(m)%lbounds, " at: ", vars(m)%addr, " | ", vars(m)%span
    !   end do
    !   write(*,*)
    ! end block

    contains

    subroutine read_journal(jstat)
      integer(kind=4), intent(out)   :: jstat
      
      character(7) :: journal_start

      call fseek(this%funit, -8, SEEK_END, jstat)
      read(this%funit) addr
      call fseek(this%funit, addr, SEEK_SET, jstat)
      read(this%funit) journal_start

      if ("JOURNAL" /= journal_start) then; jstat = ERR_NOT_FOUND; return; end if;
      read(this%funit) n

      call this%variables%init()
      allocate(vars(n))
      do i = 1, n
        call vars(i)%deserialize(this%funit)
        call this%variables%set(vars(i)%name, vars(i))
      end do
    end subroutine
  end subroutine

  subroutine storage_goto(this, name, stat)
    !! Find a variable in the storage and leaves the file position AFTER THE LENGTH OF THE BLOB
    class(storage),  intent(inout) :: this
    type(string),    intent(in)    :: name
    integer(kind=4), intent(out)   :: stat

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

  subroutine storage_close(this, stat, err_msg)
    class(storage),            intent(inout) :: this
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg
    
    integer(kind=c_int)         :: stat_
    integer(kind=int64)         :: i, start, n
    type(variable), allocatable :: variables(:)

    call init_error(stat, err_msg)

    ! block
    !   integer :: m
    !   type(variable), allocatable :: vars(:)
    !   allocate(vars(this%variables%n))
    !   call this%variables%to_array(values=vars)
    !   write(*,*)
    !   write(*,*) "Current state when closing"
    !   do m = 1, this%variables%n
    !     write(*,*) "    ", vars(m)%name%s, "(", vars(m)%type, "): ", vars(m)%sizes, " | ", vars(m)%lbounds, " at: ", vars(m)%addr, " | ", vars(m)%span
    !   end do
    !   write(*,*)
    ! end block

    if (this%mode > STORAGE_READ) then
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
      stat_ = fsync(int(fnum(this%funit), kind=c_int))
      if (stat_ /= 0) then
        m4_error(ERR_FILE_OP, "Error calling fsync")
      end if

      close(this%wal, status="delete")
    end if

    call this%variables%destruct()
    close(this%funit)
    this%funit = -1
    this%mode = 0
  end subroutine

  m4_define({T},{string})
  m4_define({U},{variable})
  m4_include(../util/map_imp.f90.inc)

  m4_define({m4_Y},{subroutine storage_write_$2_$1(this, name, var, unit, dynamic, compression, stat, err_msg)
    !! 
    class(storage),            intent(inout) :: this
    type(string),              intent(in)    :: name
    m4_type($2), target,       intent(in)    :: var{}m4_pshape($1)
    character(*),    optional, intent(in)    :: unit
    integer,         optional, intent(in)    :: dynamic
    character,       optional, intent(in)    :: compression
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    integer(kind=int64)              :: offset, data_l, header_l
    integer(kind=int64), allocatable :: sizes(:), lbounds(:), tmp(:)
    integer(kind=4)                  :: stat_
    integer                          :: dflag, ddim = -1
    character                        :: compr_

    type(variable) :: new_var
    type(mapnode_string_variable), pointer :: p => null()

    m4_ifelse($2,STRING,{
    type(string), allocatable :: bvar(:)
    },)

    call init_error(stat, err_msg)

    dflag = DYNAMIC_NO
    if (present(dynamic)) dflag = dynamic
    if (dflag < DYNAMIC_NO .or. dflag > DYNAMIC_EXT) then; m4_error(ERR_INVALID_ARGUMENTS, "Unknown dynamic flag"); end if;

    if (this%mode <= STORAGE_READ) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write in this storage") end if;

    compr_ = COMPR_NONE
    if (present(compression)) compr_ = compression
    if (compr_ == COMPR_ZLIB) then; m4_ifdef({m4_zlib},,{m4_error(ERR_INVALID_ARGUMENTS, "Cannot use zlib compression if the library is not included")}) end if;
    if (compr_ == COMPR_BLOSC) then; m4_ifdef({m4_blosc},,{m4_error(ERR_INVALID_ARGUMENTS, "Cannot use blosc compression if the library is not included")}) end if;

    if (DT_$2 /= DT_REAL64 .and. DT_$2 /= DT_CMPLX128 .and. present(unit)) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot denormalize variable of type other than real or complex") end if;

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
          ddim = size(sizes)
          if (size(sizes) /= size(p%value%sizes)) then;                                m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the dimension does not match"); end if; 
          if (any(sizes(1:ddim - 1) /= p%value%sizes(1:ddim - 1))) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the shapes do not match"); end if;
          ! Update the lbounds and sizes for writing the binary blob and add the entries to the variable
          lbounds(size(lbounds)) = p%value%lbounds(size(sizes)) - p%value%sizes(size(sizes))
          p%value%sizes(size(sizes)) = p%value%sizes(size(sizes)) - sizes(ddim)
        else if (dflag == DYNAMIC_EXT) then
          ddim = size(sizes) + 1
          if (size(sizes) /= size(p%value%sizes) - 1) then;     m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the dimension does not match"); end if;
          if (any(sizes /= p%value%sizes(1:size(sizes)))) then; m4_error(ERR_INVALID_ARGUMENTS, "Cannot write to active variable; the shape do not match"); end if;
          
          allocate(tmp(size(p%value%lbounds)))
          tmp(1:ddim-1) = lbounds
          tmp(ddim)     = p%value%lbounds(ddim) - p%value%sizes(ddim)
          call move_alloc(tmp, lbounds)

          allocate(tmp(size(p%value%sizes)))
          tmp(1:ddim-1) = sizes
          tmp(ddim)     = int(1, kind=int64)
          call move_alloc(tmp, sizes)
          
          p%value%sizes(ddim) = p%value%sizes(ddim) - 1 ! Add one entry to the last dimension
        end if
      else
        new_var%name  = name
        new_var%type  = DT_$2
        new_var%count = int(0,kind=int64)
        new_var%addr  = [offset]
        
        if (dflag == DYNAMIC_APP) then
          ddim = size(sizes)
          new_var%sizes       = sizes
          new_var%sizes(ddim) = -sizes(ddim) ! Indicate that this is a dynamic variable
          new_var%lbounds     = lbounds
          new_var%span        = [abs(sizes(size(sizes)))]
        else if (dflag == DYNAMIC_EXT) then
          ddim = size(sizes) + 1
          allocate(new_var%sizes(ddim))
          new_var%sizes(1:ddim-1) = sizes
          new_var%sizes(ddim) = int(-1, kind=int64)
          
          allocate(new_var%lbounds(ddim))
          new_var%lbounds(1:ddim-1) = lbounds
          new_var%lbounds(ddim) = int(1, kind=int64)
          
          new_var%span = [1]

          lbounds = new_var%lbounds
          sizes   = new_var%sizes
          sizes(ddim) = 1
        end if
        

        call this%variables%set(name, new_var)
        p => this%variables%find(name)
        write(this%wal, "(dt'all')") new_var
      end if

    else if (associated(p)) then
      m4_error(ERR_ALREADY_EXISTS, "The variable '"// p%value%name%s //"' already exists")
    else
      new_var%name    = name
      new_var%type    = DT_$2
      new_var%sizes   = sizes
      new_var%lbounds = lbounds
      new_var%count   = int(1,kind=int64)
      new_var%addr    = [offset]
      new_var%span    = [0]
      call this%variables%set(name, new_var)
      p => this%variables%find(name) 
    end if

    call write_blob()

    flush(this%funit)
    stat_ = fsync(int(fnum(this%funit), kind=c_int))
    if (stat_ /= 0) then
      m4_error(ERR_FILE_OP, "Error calling fsync")
    end if

    if (all(p%value%sizes > 0)) then
      write(this%wal, *) p%value
      return
    end if
    ! Must be a dynamic variable now        
    
    p%value%count = p%value%count + 1
    write(this%wal, *) p%value

    ! Resize the arrays of the variable structure
    if (p%value%count > size(p%value%addr)) then
      block
        integer(kind=int64), allocatable :: tmp(:)
        allocate (tmp(size(p%value%addr)*2))  
        tmp(1:size(p%value%addr)) = p%value%addr
        call move_alloc(tmp, p%value%addr)
        allocate (tmp(size(p%value%span)*2))  
        tmp(1:size(p%value%span)) = p%value%span
        call move_alloc(tmp, p%value%span)
      end block
    end if

    p%value%addr(p%value%count) = offset
    if (dflag == DYNAMIC_EXT) then
      p%value%span(p%value%count) = 1
    else
      p%value%span(p%value%count) = sizes(size(sizes))
    end if

    flush(this%wal)
  
  contains
    subroutine write_blob()
      character, allocatable :: data(:)
      m4_ignore(data)

      m4_ifelse($2,STRING,{
        m4_ifelse($1,0,{
          bvar = [var]},{
          bvar = reshape(var, [product(sizes)])
        })
        data_l = get_data_length(bvar)
        header_l = write_header(this%funit, name, DT_$2, lbounds, sizes, compression=COMPR_NONE, data_length=data_l)
        call write_string(this%funit, bvar)  
      },{
        data_l = DT_SIZES(DT_$2)
        if (size(lbounds) > 0) data_l = data_l*product(sizes)
        header_l = write_header(this%funit, name, DT_$2, lbounds, sizes, compression=compr_, data_length=data_l)

        if (.not. present(unit)) then
          data = transfer(var, data)
        m4_denormable($2,{else
          data = transfer(denorm(var, unit), data)  
        })
        end if

        if (compr_ == COMPR_NONE) then          
          write (this%funit) data
        m4_ifdef({m4_zlib},{elseif (compr_ == COMPR_ZLIB) then
          data_l = compress_zlib(this%funit, data)
          call fseek(this%funit, offset + len(name%s, kind=int64)+2, SEEK_SET)
          write (this%funit) header_l + data_l
        })
        m4_ifdef({m4_blosc}, {elseif (compr_ == COMPR_BLOSC) then
          data_l = compress_blosc(this%funit, data, int(DT_SIZES(DT_$2), kind=c_size_t))
          call fseek(this%funit, offset + len(name%s, kind=int64)+2, SEEK_SET)
          write (this%funit) header_l + data_l
        })
        else
          m4_error(ERR_INVALID_ARGUMENTS, "Undefined compression argument " // compr_)
        end if
      })
    end subroutine
  end subroutine

  subroutine storage_writec_$2_$1(this, name, var, unit, dynamic, compression, stat, err_msg)
    class(storage),            intent(inout) :: this
    character(*),              intent(in)    :: name
    m4_type($2),               intent(in)    :: var{}m4_pshape($1)
    character(*),    optional, intent(in)    :: unit
    integer,         optional, intent(in)    :: dynamic
    character,       optional, intent(in)    :: compression
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    call this%write(new_string(name), var, unit, dynamic, compression, stat, err_msg)
  end subroutine})
  m4_list

  m4_define({m4_Y},{subroutine storage_read_$2_$1(this, name, var, index, stat, err_msg)
    class(storage),            intent(inout) :: this
    type(string),              intent(in)    :: name
    m4_type($2) m4_ifelse($1,0,,{, allocatable}), target, intent(out) :: var{}m4_pshape($1)
    integer,         optional, intent(in)    :: index
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    integer(kind=4)                  :: stat_
    m4_type($2), contiguous, pointer :: pvar({}m4_pdimm1($1):)
    m4_ifelse($1,0,,{
    integer(kind=int64) :: i, idx
    integer(kind=int64), allocatable :: sizes(:), lbounds(:), ubounds(:)
    })
    type(mapnode_string_variable), pointer :: p => null()

    call init_error(stat, err_msg)
    
    ! Implies that the variable is registered
    call this%goto(name, stat_)
    if (stat_ /= 0) then
      m4_error(ERR_NOT_FOUND, "Could not find the variable " // name%s)
    end if
    p => this%variables%find(name)

    m4_ifelse($1,0,{
    allocate(pvar(1:1))
    call read_blob_$2_$1(pvar)
    var = pvar(1)
    
    },{
    ! Allocate the array with the complete shape
    sizes   = p%value%sizes
    lbounds = p%value%lbounds
    sizes(size(sizes)) = abs(sizes(size(sizes)))
    ubounds = lbounds + sizes - 1

    if (present(index)) then
      ! Search for index in blobs
      idx = lbounds(size(lbounds))
      do i = 1, p%value%count
        if (idx <= index .and. idx + p%value%span(i) - 1 >= index) exit
        idx = idx + p%value%span(i)
      end do
      lbounds(size(lbounds)) = idx
      ubounds(size(ubounds)) = idx + p%value%span(i) - 1
      sizes(size(sizes))     = p%value%span(i)
      allocate(var({}m4_pallocate($1)))

      call fseek(this%funit, p%value%addr(i), SEEK_SET)
      pvar => var
      call read_blob_$2_$1(pvar)
      return
    end if

    allocate(var({}m4_pallocate($1)))

    if (all(p%value%sizes > 0)) then ! Static array
      pvar => var
      call read_blob_$2_$1(pvar)
      return
    end if

    idx = lbounds(size(lbounds))
    do i = 1, p%value%count
      call fseek(this%funit, p%value%addr(i), SEEK_SET)
      pvar({}m4_pallocate(m4_decr($1)) m4_ifelse($1,1,,{,}) idx:idx+p%value%span(i)-1) => var({}m4_pdimm1($1) idx:idx+p%value%span(i)-1)
      call read_blob_$2_$1(pvar)
      idx = idx + p%value%span(i)
    end do
    })
  contains
    subroutine read_blob_$2_$1(data)
      m4_type($2), contiguous, pointer, intent(inout) :: data({}m4_pdimm1($1):)

      integer                          :: dtype
      integer(kind=int64)              :: data_length, pos
      m4_ifelse($2,STRING,,{integer(kind=int64) :: var_length})
      character(len=:),    allocatable :: bunit
      character                        :: compr
      integer(kind=int64), allocatable :: form(:), lower_bounds(:)

      call read_header(this%funit, type=dtype, unit=bunit, lbounds=lower_bounds, sizes=form, compression=compr, data_length=data_length)
      call ftell(this%funit, pos)
      
      if (dtype /= DT_$2) then
        m4_error(ERR_INVALID_ARGUMENTS, "The variable type '" // DT_NAMES(dtype) // "' does not match the called procedure ($2)")
      end if

      if (compr == COMPR_ZLIB) then; m4_ifdef({m4_zlib},,{m4_error(ERR_INVALID_ARGUMENTS, "Cannot decompress from zlib compression if the library is not included")}) end if;
      if (compr == COMPR_BLOSC) then; m4_ifdef({m4_blosc},,{m4_error(ERR_INVALID_ARGUMENTS, "Cannot decompress from blosc compression if the library is not included")}) end if;

      m4_ifelse($2,STRING,{
        block
          type(string), contiguous, pointer :: pstr(:)
          pstr(1:product(form)) => data
          call read_string(this%funit, pstr)
        end block
        m4_ignore(index)
      },{
        m4_ifelse($1,0,,{
        if (.not. all(lbound(data, kind=int64) == lower_bounds)) then
          m4_error(ERR_INVALID_ARGUMENTS, "Allocated lower bounds for " // name%s // " do not match")
        end if

        if (.not. all(shape(data, kind=int64) == form)) then
          m4_error(ERR_INVALID_ARGUMENTS, "Allocated shape for " // name%s // " does not match")
        end if
        })
      
        var_length = DT_SIZES(DT_$2)
        if (size(form) > 0) var_length = var_length * product(form)

        m4_ifdef({m4_zlib},{if (compr == COMPR_ZLIB) then
          data = reshape(transfer(decompress_zlib(this%funit, data_length, var_length), data), shape(data))
        end if})
        m4_ifdef({m4_blosc},{if (compr == COMPR_BLOSC) then
          data = reshape(transfer(decompress_blosc(this%funit, data_length, var_length), data), shape(data))
        end if})
        if (compr == COMPR_NONE) then
          read(this%funit) data
        end if
        m4_ignore(index)
      })

      m4_ifelse($2,REAL64,{
      if (len(bunit) > 0) then
        data = norm(data, bunit)
      end if},{m4_ifelse($2,CMPLX128,{
      if (len(bunit) > 0) then
        data = norm(data, bunit)
      end if})})

      data_length = data_length + pos
      call ftell(this%funit, pos)
      data_length = pos - data_length
      if (data_length /= 0) then
        m4_error(ERR_FILE_OP, "Read too many bytes: " // int2str(int(data_length)))
      end if 
    end subroutine
  end subroutine
  subroutine storage_readc_$2_$1(this, name, var, index, stat, err_msg)
    class(storage),            intent(inout) :: this
    character(*),              intent(in)    :: name
    m4_type($2) m4_ifelse($1,0,,{, allocatable}), intent(out) :: var{}m4_pshape($1)
    integer,         optional, intent(in)    :: index
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    call this%read(new_string(name), var, index, stat, err_msg)
  end subroutine})
  m4_list

  m4_ifdef({m4_zlib},{function compress_zlib(funit, data, stat, err_msg) result(data_length)
    integer,                     intent(in)  :: funit
    character, target,           intent(in)  :: data(:)
    integer,           optional, intent(out) :: stat
    type(string),      optional, intent(out) :: err_msg
    integer(kind=int64)                      :: data_length

    character, target, allocatable :: out(:)
    type(z_stream)                 :: zstr
    integer(kind=c_int)            :: ret

    zstr%avail_in = int(sizeof(data), kind=z_uint)
    zstr%next_in = c_loc(data)
    ret = deflate_init(zstr, Z_BEST_COMPRESSION)
    if (ret /= Z_OK) then; m4_error(ERR_INTERNAL, "Cannot initialize deflate"); end if;
    
    data_length = deflate_bound(zstr, sizeof(data))
    
    allocate(out(data_length))
    zstr%avail_out = int(data_length, kind=z_uint)
    zstr%next_out  = c_loc(out)

    ret = deflate(zstr, Z_FINISH)
    if (ret /= Z_STREAM_END) then; m4_error(ERR_INTERNAL, "Not all input has been consumed"); end if;
    
    data_length = zstr%total_out
    ret = deflate_end(zstr)

    write (funit) out(:data_length)
  end function
  
  function decompress_zlib(funit, length, var_length, stat, err_msg) result(out)
    integer,                intent(in)  :: funit
    integer(kind=int64),    intent(in)  :: length
    integer(kind=int64),    intent(in)  :: var_length
    integer,      optional, intent(out) :: stat
    type(string), optional, intent(out) :: err_msg

    character, target :: out(var_length)
    character, target :: in(length)
    type(z_stream)    :: zstr
    integer           :: ret

    read (funit) in

    zstr%avail_in  = int(length, kind=z_uint)
    zstr%next_in   = c_loc(in)
    zstr%avail_out = int(var_length, kind=z_uint)
    zstr%next_out  = c_loc(out)

    ret = inflate_init(zstr)
    if (ret /= Z_OK) then; m4_error(ERR_INTERNAL, "Cannot initialize inflate"); end if;

    ret = inflate(zstr, Z_FINISH)
    if (ret /= Z_STREAM_END) then; m4_error(ERR_INTERNAL, "Could not inflate the data: " // int2str(ret)); end if;
    
    if (var_length /= zstr%total_out) then
      m4_error(ERR_INTERNAL, "Inflated variable size not as expected")
    end if

  end function})

  m4_ifdef({m4_blosc},{function compress_blosc(funit, data, typesize, stat, err_msg) result(data_length)
    integer,                intent(in)  :: funit
    character, target,      intent(in)  :: data(:)
    integer(kind=c_size_t), intent(in)  :: typesize
    integer,      optional, intent(out) :: stat
    type(string), optional, intent(out) :: err_msg
    integer(kind=int64)                 :: data_length
    
    character(kind=int8), target, allocatable :: out(:)
    type(c_ptr)            :: data_ptr
    integer(kind=c_size_t) :: isize

    call blosc_init()

    isize = int(sizeof(data), kind=c_size_t)
    data_ptr = c_loc(data)

    allocate(out(isize + BLOSC_MIN_HEADER_LENGTH))
    
    data_length = blosc_compress(int(9, kind=c_int), BLOSC_BITSHUFFLE, typesize, isize, data_ptr, c_loc(out), isize + BLOSC_MIN_HEADER_LENGTH)
    if (data_length == 0) then; m4_error(ERR_INTERNAL, "The variable is not compressible"); end if;
    if (data_length < 0) then; m4_error(ERR_INTERNAL, "Got error on blosc compression: " // int2str(int(data_length))); end if;

    write (funit) out(:data_length)
    call blosc_destroy()
  end function

  function decompress_blosc(funit, length, var_length, stat, err_msg) result(out)
    integer,                intent(in)  :: funit
    integer(kind=int64),    intent(in)  :: length
    integer(kind=int64),    intent(in)  :: var_length
    integer,      optional, intent(out) :: stat
    type(string), optional, intent(out) :: err_msg
    character, target                   :: out(var_length)

    character, target   :: in(length)
    integer(kind=c_int) :: ret

    read (funit) in

    call blosc_init()

    ret = blosc_decompress(c_loc(in), c_loc(out), var_length)
    if (ret <= 0) then; m4_error(ERR_INTERNAL, "Cannot decompress: " // int2str(ret)); end if;

    if (ret /= var_length) then
      m4_error(ERR_INTERNAL, "Decompressed variable size not as expected")
    end if

    call blosc_destroy()

  end function})

end module