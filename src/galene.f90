m4_include(util/macro.f90.inc)

module galene_m

  use, intrinsic :: iso_fortran_env, only: int32, int64

  use error_m,         only: program_error
  use map_m,           only: map_string_int, map_string_string
  use normalization_m, only: norm
  use string_m,        only: string

  implicit none

  private
  public gal_block, gal_file
  public GALDATA_INT, GALDATA_REAL, GALDATA_CHAR

  integer, parameter :: GALDATA_INT  = 1
  integer, parameter :: GALDATA_REAL = 2
  integer, parameter :: GALDATA_CHAR = 3

  type gal_block
    !! GALENE III data block

    type(string) :: name
      !! identifier
    type(string) :: unit
      !! physical unit

    integer :: dtype
      !! type of data elements
    integer :: ndata
      !! number of data elements

    integer,      allocatable :: idata(:)
      !! integer data
    real,         allocatable :: rdata(:)
      !! real data
    character(:), allocatable :: cdata
      !! character data
  contains
    procedure :: read => gal_block_read
  end type

  m4_define({T}, {gal_block})
  m4_include(util/vector_def.f90.inc)

  type gal_file
    !! GALENE III output file

    type(vector_gal_block) :: blocks
    type(map_string_int)   :: blmap

  contains
    procedure :: init      => gal_file_init
    procedure :: get_block => gal_file_get_block
  end type

  type(map_string_string) :: units

contains

  m4_define({T}, {gal_block})
  m4_include(util/vector_imp.f90.inc)

  subroutine gal_block_read(this, funit, dtype)
    class(gal_block), intent(out) :: this
    integer,          intent(in)  :: funit
      !! file unit
    integer,          intent(in)  :: dtype
      !! data type

    integer(int32) :: n1, n2
    integer        :: nch
    ! m4_ifelse(m4_intsize,64,{
    integer(int32), allocatable :: tmp(:)
    ! },{})
    logical :: status

    this%dtype = dtype

    ! read number of data elements + size of name
    read (funit) n1, n2
    this%ndata = int(n1)
    nch        = int(n2)

    ! read name
    allocate (character(nch) :: this%name%s)
    read (funit) this%name%s

    ! read data
    select case (dtype)
    case (GALDATA_INT)
      allocate (this%idata(this%ndata))
      if (this%ndata > 0) then
      ! m4_ifelse(m4_intsize,64,{
        allocate (tmp(this%ndata))
        read (funit) tmp
        this%idata = int(tmp, int64)
      ! },{
        read (funit) this%idata
      ! })
      end if

    case (GALDATA_REAL)
      allocate (this%rdata(this%ndata))
      if (this%ndata > 0) then
        read (funit) this%rdata
        call units%get(this%name, this%unit, status = status)
        if (status) then
          this%rdata = norm(this%rdata, this%unit%s)
        else
          this%unit%s = "1"
        end if
      end if

    case (GALDATA_CHAR)
      allocate (character(this%ndata) :: this%cdata)
      if (this%ndata > 0) then
        read (funit) this%cdata
      end if

    end select
  end subroutine

  subroutine gal_file_init(this, file)
    class(gal_file), intent(out) :: this
    character(*),    intent(in)  :: file
      !! filename

    integer         :: funit, iostat, dtype
    type(gal_block) :: b

    if (units%root <= 0) then
      ! add units for SOME quantities
      call units%init()
      call units%set(string("n-temp"), string("K"))
      call units%set(string("p-temp"), string("K"))
      call units%set(string("n-density"), string("1/cm^3"))
      call units%set(string("p-density"), string("1/cm^3"))
      call units%set(string("n-imref"), string("V"))
      call units%set(string("p-imref"), string("V"))
      call units%set(string("n-velocity"), string("cm/s"))
      call units%set(string("p-velocity"), string("cm/s"))
      call units%set(string("n-mobility"), string("cm^2/V/s"))
      call units%set(string("p-mobility"), string("cm^2/V/s"))
      call units%set(string("electron-current"), string("A/cm^2"))
      call units%set(string("hole-current"), string("A/cm^2"))
      call units%set(string("electron-energy"), string("J/(cm^2*s)"))
      call units%set(string("hole-energy"), string("J/(cm^2*s)"))
      call units%set(string("potential"), string("V"))
      call units%set(string("avalanche-generation"), string("1/(cm^3*s)"))
      call units%set(string("srh-recombination"), string("1/(cm^3*s)"))
      call units%set(string("blank.dd"), string("1"))
      call units%set(string("blank.td"), string("1"))
      call units%set(string("ni"), string("1/cm^3"))
      call units%set(string("nv"), string("1/cm^3"))
      call units%set(string("ec"), string("eV"))
      call units%set(string("ev"), string("eV"))
      call units%set(string("delta-ec"), string("eV"))
      call units%set(string("delta-ev"), string("eV"))
      call units%set(string("eg"), string("eV"))
      call units%set(string("ge-concentration"), string("1"))
      call units%set(string("mc-pot"), string("V"))
      call units%set(string("electric-field"), string("kV/cm"))
      call units%set(string("x-coord"), string("um"))
      call units%set(string("y-coord"), string("um"))
    end if

    call this%blocks%init(0, c = 32)
    call this%blmap%init()

    open (newunit = funit, file = file, form = "unformatted", status = "old", action = "read")

    do while (.true.)
      ! read type of following data
      read (funit, iostat = iostat) dtype

      ! check for end of file
      if (iostat /= 0) exit

      ! check for end of g3-record
      if (dtype < 0) cycle

      ! check if type is valid
      if ((dtype > 3) .or. (dtype == 0)) call program_error("Error in gal3 binary file: unknown data type")

      ! read block
      call b%read(funit, dtype)

      ! save block + index
      call this%blocks%push(b)
      call this%blmap%set(b%name, this%blocks%n)
    end do

    close (funit)
  end subroutine

  function gal_file_get_block(this, name) result(b)
    class(gal_file), target, intent(in) :: this
    character(*),            intent(in) :: name
      !! block name
    type(gal_block), pointer            :: b
      !! return pointer to galene data block

    integer :: i
    logical :: status

    call this%blmap%get(string(name), i, status = status)
    nullify (b)
    if (status) b => this%blocks%d(i)
  end function

end module
