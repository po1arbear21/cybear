m4_include(macro.f90.inc)

module output_file_m

  use error_m,         only: assert_failed, program_error
  use iso_fortran_env, only: int32, int64
  use json_m,          only: json, json_array, json_object, json_save
  use normalization_m, only: denorm

  implicit none

  private
  public output_file

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

    type(json_object) :: obj
      !! main json object
  contains
    procedure :: init       => output_file_init
    procedure :: save       => output_file_save
    procedure :: close      => output_file_close
    procedure :: is_open    => output_file_is_open
    procedure :: new_object => output_file_new_object
    generic   :: write      => output_file_write_i32_1,   output_file_write_i32_2,   output_file_write_i32_3, &
      &                        output_file_write_i64_1,   output_file_write_i64_2,   output_file_write_i64_3, &
      &                        output_file_write_real_1,  output_file_write_real_2,  output_file_write_real_3, &
      &                        output_file_write_cmplx_1, output_file_write_cmplx_2, output_file_write_cmplx_3

    procedure, private :: output_file_write_i32_1,   output_file_write_i32_2,   output_file_write_i32_3
    procedure, private :: output_file_write_i64_1,   output_file_write_i64_2,   output_file_write_i64_3
    procedure, private :: output_file_write_real_1,  output_file_write_real_2,  output_file_write_real_3
    procedure, private :: output_file_write_cmplx_1, output_file_write_cmplx_2, output_file_write_cmplx_3
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

    ! initialize json object
    call this%obj%init()
    call this%obj%add("DataFile", name // ".bin")
  end subroutine

  subroutine output_file_save(this)
    !! save json object to file
    class(output_file), intent(in) :: this

    call json_save(this%obj, this%fname_json)
  end subroutine

  subroutine output_file_close(this)
    !! save json object and close binary file
    class(output_file), intent(inout) :: this

    m4_assert(this%is_open())

    call this%save()
    call this%obj%destruct()

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

    call this%obj%get_json(category, js)
    if (associated(js)) then
      select type (js)
      type is (json_array)
        js_category => js
      end select
    else
      allocate (js_category)
      call js_category%init()
      call this%obj%add(category, js_category)
    end if

    allocate (obj)
    call obj%init()
    call js_category%add(obj)
  end function

  subroutine output_file_write_i32_1(this, obj, name, values)
    !! write 1D 32-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int32),        intent(in)    :: values(:)
      !! 1D 32-bit integer array

    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i32")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 4
  end subroutine

  subroutine output_file_write_i32_2(this, obj, name, values)
    !! write 2D 32-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int32),        intent(in)    :: values(:,:)
      !! 2D 32-bit integer array

    integer                    :: i
    integer, allocatable       :: tmp(:)
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    tmp = shape(values)
    do i = 1, size(tmp)
      call sh%add(tmp(i))
    end do

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i32")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 4
  end subroutine

  subroutine output_file_write_i32_3(this, obj, name, values)
    !! write 3D 32-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int32),        intent(in)    :: values(:,:,:)
      !! 3D 32-bit integer array

    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))
    call sh%add(size(values, 3))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i32")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 4
  end subroutine

  subroutine output_file_write_i64_1(this, obj, name, values)
    !! write 1D 64-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int64),        intent(in)    :: values(:)
      !! 1D 64-bit integer array

    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i64")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_i64_2(this, obj, name, values)
    !! write 2D 64-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int64),        intent(in)    :: values(:,:)
      !! 2D 64-bit integer array

    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i64")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_i64_3(this, obj, name, values)
    !! write 3D 64-bit integer array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    integer(kind=int64),        intent(in)    :: values(:,:,:)
      !! 3D 64-bit integer array

    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))
    call sh%add(size(values, 3))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "i64")
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) values
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_real_1(this, obj, name, values, unit)
    !! write 1D real array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    real,                       intent(in)    :: values(:)
      !! 1D real array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "f64")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) denorm(values, unit_)
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_real_2(this, obj, name, values, unit)
    !! write 2D real array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    real,                       intent(in)    :: values(:,:)
      !! 2D real array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "f64")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) denorm(values, unit_)
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_real_3(this, obj, name, values, unit)
    !! write 3D real array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    real,                       intent(in)    :: values(:,:,:)
      !! 3D real array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    integer                    :: i
    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))
    call sh%add(size(values, 3))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "f64")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    do i = 1, size(values, 3)
      write (this%funit_bin) denorm(values(:,:,i), unit_)
    end do
    this%index_bin = this%index_bin + size(values) * 8
  end subroutine

  subroutine output_file_write_cmplx_1(this, obj, name, values, unit)
    !! write 1D complex array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    complex,                    intent(in)    :: values(:)
      !! 1D complex array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "c128")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) denorm(values, unit_)
    this%index_bin = this%index_bin + size(values) * 16
  end subroutine

  subroutine output_file_write_cmplx_2(this, obj, name, values, unit)
    !! write 2D complex array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    complex,                    intent(in)    :: values(:,:)
      !! 2D complex array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "c128")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    write (this%funit_bin) denorm(values, unit_)
    this%index_bin = this%index_bin + size(values) * 16
  end subroutine

  subroutine output_file_write_cmplx_3(this, obj, name, values, unit)
    !! write 3D complex array to binary file
    class(output_file),         intent(inout) :: this
    type(json_object), pointer, intent(inout) :: obj
      !! parent json object
    character(*),               intent(in)    :: name
      !! data name
    complex,                    intent(in)    :: values(:,:,:)
      !! 3D complex array
    character(*), optional,     intent(in)    :: unit
      !! physical unit token; default = "1"

    integer                    :: i
    character(:), allocatable  :: unit_
    type(json_array),  pointer :: sh
    type(json_object), pointer :: dat

    unit_ = "1"
    if (present(unit)) unit_ = unit

    ! create json shape array
    allocate (sh)
    call sh%init()
    call sh%add(size(values, 1))
    call sh%add(size(values, 2))
    call sh%add(size(values, 3))

    ! create json data descriptor and add it to parent object
    allocate (dat)
    call dat%init()
    call dat%add("Type", "c128")
    call dat%add("Unit", unit_)
    call dat%add("Shape", sh)
    call dat%add("Index", this%index_bin)
    call obj%add(name, dat)

    ! write data to binary file
    do i = 1, size(values, 3)
      write (this%funit_bin) denorm(values(:,:,i), unit_)
    end do
    this%index_bin = this%index_bin + size(values) * 16
  end subroutine

end module
