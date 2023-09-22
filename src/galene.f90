m4_include(util/macro.f90.inc)

module galene_m

  use error_m,         only: program_error
  use iso_fortran_env, only: int32, int64
  use map_m,           only: map_string_int, mapnode_string_int, map_string_string, mapnode_string_string
  use normalization_m, only: norm
  use string_m,        only: new_string

  implicit none

  private
  public gal_block, gal_file
  public GALDATA_INT, GALDATA_REAL, GALDATA_CHAR

  integer, parameter :: GALDATA_INT  = 1
  integer, parameter :: GALDATA_REAL = 2
  integer, parameter :: GALDATA_CHAR = 3

  type gal_block
    !! GALENE III data block

    character(:), allocatable :: name
      !! identifier
    character(:), allocatable :: unit
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
    m4_ifelse(m4_intsize,64,{
      integer(int32), allocatable :: tmp(:)
    },{})
    type(mapnode_string_string), pointer :: nd

    this%dtype = dtype

    ! read number of data elements + size of name
    read (funit) n1, n2
    this%ndata = int(n1)
    nch        = int(n2)

    ! read name
    allocate (character(nch) :: this%name)
    read (funit) this%name

    ! read data
    select case (dtype)
    case (GALDATA_INT)
      allocate (this%idata(this%ndata))
      if (this%ndata > 0) then
        m4_ifelse(m4_intsize,64,{
          allocate (tmp(this%ndata))
          read (funit) tmp
          this%idata = int(tmp, int64)
        },{
          read (funit) this%idata
        })
      end if

    case (GALDATA_REAL)
      allocate (this%rdata(this%ndata))
      if (this%ndata > 0) then
        read (funit) this%rdata
        nd => units%find(new_string(this%name))
        if (associated(nd)) then
          this%unit = nd%value%s
          this%rdata = norm(this%rdata, this%unit)
        else
          this%unit = "1"
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

    if (.not. associated(units%root)) then
      ! add units for SOME quantities
      call units%init()
      call units%insert(new_string("n-temp"), new_string("K"))
      call units%insert(new_string("p-temp"), new_string("K"))
      call units%insert(new_string("n-density"), new_string("1/cm^3"))
      call units%insert(new_string("p-density"), new_string("1/cm^3"))
      call units%insert(new_string("n-imref"), new_string("V"))
      call units%insert(new_string("p-imref"), new_string("V"))
      call units%insert(new_string("n-velocity"), new_string("cm/s"))
      call units%insert(new_string("p-velocity"), new_string("cm/s"))
      call units%insert(new_string("n-mobility"), new_string("cm^2/V/s"))
      call units%insert(new_string("p-mobility"), new_string("cm^2/V/s"))
      call units%insert(new_string("electron-current"), new_string("A/cm^2"))
      call units%insert(new_string("hole-current"), new_string("A/cm^2"))
      call units%insert(new_string("electron-energy"), new_string("J/(cm^2*s)"))
      call units%insert(new_string("hole-energy"), new_string("J/(cm^2*s)"))
      call units%insert(new_string("potential"), new_string("V"))
      call units%insert(new_string("avalanche-generation"), new_string("1/(cm^3*s)"))
      call units%insert(new_string("srh-recombination"), new_string("1/(cm^3*s)"))
      call units%insert(new_string("blank.dd"), new_string("1"))
      call units%insert(new_string("blank.td"), new_string("1"))
      call units%insert(new_string("ni"), new_string("1/cm^3"))
      call units%insert(new_string("nv"), new_string("1/cm^3"))
      call units%insert(new_string("ec"), new_string("eV"))
      call units%insert(new_string("ev"), new_string("eV"))
      call units%insert(new_string("delta-ec"), new_string("eV"))
      call units%insert(new_string("delta-ev"), new_string("eV"))
      call units%insert(new_string("eg"), new_string("eV"))
      call units%insert(new_string("ge-concentration"), new_string("1"))
      call units%insert(new_string("mc-pot"), new_string("V"))
      call units%insert(new_string("electric-field"), new_string("kV/cm"))
      call units%insert(new_string("x-coord"), new_string("um"))
      call units%insert(new_string("y-coord"), new_string("um"))
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
      call this%blmap%insert(new_string(b%name), this%blocks%n)
    end do

    close (funit)
  end subroutine

  function gal_file_get_block(this, name) result(b)
    class(gal_file), target, intent(in) :: this
    character(*),            intent(in) :: name
      !! block name
    type(gal_block), pointer            :: b
      !! return pointer to galene data block

    type(mapnode_string_int), pointer :: nd

    nd => this%blmap%find(new_string(name))

    nullify(b)
    if (associated(nd)) b => this%blocks%d(nd%value)
  end function

end module
