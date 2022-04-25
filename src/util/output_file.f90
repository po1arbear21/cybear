m4_include(macro.f90.inc)

module output_file_m

  use error_m,         only: assert_failed, program_error
  use iso_fortran_env, only: int32, int64
  use json_m,          only: json_file, json, json_array, json_object
  use normalization_m, only: denorm

  implicit none

  private
  public output_file

  ! output types
  m4_define({m4_typelist},{
    m4_X(int32)
    m4_X(int64)
    m4_X(log)
    m4_X(real)
    m4_X(cmplx)
  })

  ! real or complex (use physical unit)
  m4_define({m4_real_or_cmplx},{m4_ifelse($1,real,{$2},{m4_ifelse($1,cmplx,{$2},{$3})})})

  ! typename
  m4_define({m4_typename},{m4_ifelse($1,int32,i32,{m4_ifelse($1,int64,i64,{m4_ifelse($1,log,l32,{m4_ifelse($1,real,f64,{m4_ifelse($1,cmplx,c128)})})})})})

  ! size in bytes
  m4_define({m4_sizeof},{m4_ifelse($1,int32,4,{m4_ifelse($1,int64,8,{m4_ifelse($1,log,4,{m4_ifelse($1,real,8,{m4_ifelse($1,cmplx,16)})})})})})

  ! dimensions (0..max_dim)
  m4_define({m4_max_dim},{8})
  m4_define({m4_dimlist},{
    m4_ifelse($1,0,,{m4_dimlist(m4_decr($1),$2)})
    m4_Y($1,$2)
  })

  ! combine type-list with dimension-list to create full list
  m4_define({m4_X},{m4_dimlist(m4_max_dim,$1)})
  m4_define({m4_list},m4_typelist)

  type output_file
    !! output file handle (represents two files, json and binary)

    character(:), allocatable :: name
      !! file name (without extension)
    character(:), allocatable :: folder
      !! folder name

    character(:), allocatable :: fname_json
      !! full name of json file
    character(:), allocatable :: fname_bin
      !! full name of binary file

    integer :: funit_bin = 0
      !! file unit of binary file

    integer :: index_bin
      !! index in binary file

    type(json_file)            :: jsfile
      !! json file
    type(json_object), pointer :: obj => null()
      !! pointer to main json object
  contains
    procedure :: init       => output_file_init
    procedure :: save       => output_file_save
    procedure :: close      => output_file_close
    procedure :: is_open    => output_file_is_open
    procedure :: new_object => output_file_new_object

    m4_define({m4_Y},{generic :: write => output_file_write_$2_$1})
    m4_list

    m4_define({m4_Y},{procedure, private :: output_file_write_$2_$1})
    m4_list
  end type

contains

  subroutine output_file_init(this, name, folder)
    !! initialize output file
    class(output_file),     intent(out) :: this
    character(*),           intent(in)  :: name
      !! file name without extension (e.g. "my_output")
    character(*), optional, intent(in)  :: folder
      !! folder name including trailing slash; default = "./"

    character(80) :: iomsg
    integer       :: iostat

    ! set name and folder
    this%name   = name
    this%folder = "./"
    if (present(folder)) this%folder = folder

    ! actual (full) filenames
    this%fname_json = this%folder // name // ".json"
    this%fname_bin  = this%folder // name // ".bin"

    ! open binary file for writing
    open (newunit = this%funit_bin, file = this%fname_bin, access = "stream", form = "unformatted", status = "replace", action = "write", iostat = iostat, iomsg = iomsg)
    if (iostat /= 0) call program_error("Error while opening file '" // this%fname_bin // "': " // iomsg)
    this%index_bin = 0

    ! initialize json file
    call this%jsfile%init()
    this%obj => this%jsfile%js%cast_object()
    call this%obj%add_string("DataFile", name // ".bin")
  end subroutine

  subroutine output_file_save(this)
    !! save json file
    class(output_file), intent(in) :: this

    call this%jsfile%save(this%fname_json)
  end subroutine

  subroutine output_file_close(this)
    !! save json file and close binary file
    class(output_file), intent(inout) :: this

    m4_assert(this%is_open())

    call this%save()
    call this%jsfile%destruct()

    close (this%funit_bin)
    this%funit_bin = -1
  end subroutine

  function output_file_is_open(this) result(o)
    !! determine if output file is open
    class(output_file), intent(in) :: this
    logical                        :: o
      !! return true if file is open, otherwise false

    inquire (unit = this%funit_bin, opened = o)
  end function

  function output_file_new_object(this, category) result(obj)
    !! create new json object in specified category
    class(output_file), intent(inout) :: this
    character(*),       intent(in)    :: category
      !! category name
    type(json_object), pointer        :: obj
      !! pointer to newly created json object

    class(json),      pointer :: js
    type(json_array), pointer :: js_category

    js => this%obj%get_json(category)
    if (associated(js)) then
      js_category => js%cast_array()
    else
      call this%obj%add_array(category, p = js_category)
    end if

    call js_category%add_object(p = obj)
  end function

  m4_define({m4_Y},{subroutine output_file_write_$2_$1(this, obj, name, values{}m4_real_or_cmplx($2,{, unit}))
    !! write data to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    m4_type($2),                intent(in)    :: values{}m4_pshape($1)
      !! data values
    m4_real_or_cmplx($2,{
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"
    })

    integer                    :: i
    integer, allocatable       :: tmp(:)
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    m4_real_or_cmplx($2,{
    character(:), allocatable :: unit_

    unit_ = "1"
    if (present(unit)) unit_ = unit
    })

    call obj%add_object(name, p = dat)
    call dat%add_string("Type", "m4_typename($2)")
    call dat%add_array("Shape", p = sh)
    tmp = shape(values)
    do i = 1, size(tmp)
      call sh%add_int(tmp(i))
    end do
    call dat%add_int("Index", this%index_bin)

    ! write data to binary file
    write (this%funit_bin) m4_real_or_cmplx($2,{denorm(values, unit_)},{values})
    this%index_bin = this%index_bin + m4_ifelse($1,0,{1},{size(values)}) * m4_sizeof($2)
  end subroutine})
  m4_list

end module
