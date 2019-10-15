module src_file_m
  use, intrinsic :: iso_fortran_env
  implicit none

  ! string length
  integer, parameter :: SLEN = 256

#define T string
#define TT character(len=SLEN)
#include "vector_def.f90.inc"

#define T int
#define TT integer
#include "vector_def.f90.inc"

  type src_file
    character(len=SLEN) :: folder
    character(len=SLEN) :: filename
    character(len=SLEN) :: objname
    character(len=SLEN) :: ancname

    character(len=SLEN) :: program_name
    type(vector_string) :: modules
    type(vector_string) :: submodules
    type(vector_string) :: use_modules
    type(vector_string) :: includes

    ! indices for dependencies
    integer, allocatable :: dep_modules(:)
    type(vector_int)     :: dep_anchor
  contains
    procedure :: init => src_file_init
  end type

  type src_item
    character(len=SLEN) :: name
    integer             :: file_index
  end type

#define T src_item
#define TT type(src_item)
#include "vector_def.f90.inc"

contains

#define T string
#define TT character(len=SLEN)
#include "vector_imp.f90.inc"

#define T int
#define TT integer
#include "vector_imp.f90.inc"

#define T src_item
#define TT type(src_item)
#include "vector_imp.f90.inc"

  subroutine src_file_init(this, folder, filename)
    class(src_file),     intent(out) :: this
    character(len=SLEN), intent(in)  :: folder
    character(len=SLEN), intent(in)  :: filename

    ! local variables
    integer             :: funit, line_number, status, i
    character(len=SLEN) :: line, lfmt

    ! init members
    this%folder       = folder
    this%filename     = filename
    this%program_name = ""
    call this%modules%init(    0, c = 8)
    call this%submodules%init( 0, c = 8)
    call this%use_modules%init(0, c = 8)
    call this%dep_anchor%init( 0, c = 8)

    ! object and anchor name
    i = scan(filename, '.', back = .true.)
    if (i == 0) call error("filename must be *.f90")
    if (filename(i:len_trim(filename)) /= ".f90") call error("filename must be *.f90")
    this%objname = trim(filename(1:i))//"o"
    this%ancname = trim(filename(1:i))//"anc"

    ! line format
    write(lfmt, "(A,I0,A)") "(A", SLEN, ")"

    ! open file
    open(newunit = funit, file = trim(folder)//trim(filename), status = "old", action = "read")

    ! start at the top
    line_number = 0

    ! read and parse lines to get includes
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
        call error("IO Error")
      end if
      line_number = line_number + 1

      ! remove leading and trailing spaces
      line = trim(adjustl(line))

      ! remove comments
      block
        character :: quote

        quote = ""

        ! got through line and search for comment character
        do i = 1, len_trim(line)
          ! check if inside of quoted string
          if (quote == "") then ! not inside quoted string
            ! check for beginning of quoted string or comment
            if ((line(i:i) == "'") .or. (line(i:i) == '"')) then
              quote = line(i:i)
            elseif (line(i:i) == "!") then
              ! trim line up to comment character
              line = trim(line(1:i-1))
              exit
            end if
          elseif (line(i:i) == quote) then
            ! end of quoted string
            quote = ""
          end if
        end do
      end block

      ! check for program statement
      if (to_lower(line(1:8)) == "program ") then
        block
          character(len=SLEN) :: program_name

          if (len_trim(line) < 9) call error("program statement error")
          program_name = adjustl(line(9:len_trim(line)))

          if (this%program_name /= "") call error("multiple program statements in one file")
          this%program_name = to_lower(program_name)
        end block
      end if

      ! check for module statements
      if (to_lower(line(1:7)) == "module ") then
        block
          character(len=SLEN) :: module_name

          if (len_trim(line) < 8) call error("module specification error")
          module_name = trim(adjustl(line(8:len_trim(line))))

          ! make sure this is not a "module procedure ..." etc. statement
          if (scan(module_name(1:len_trim(module_name)), " ") == 0) then
            call push_unique(this%modules, to_lower(module_name))
            cycle
          end if
        end block
      end if

      ! check for submodule statement
      if ((to_lower(line(1:10)) == "submodule ") .or. &
          (to_lower(line(1:10)) == "submodule(")) then
        block
          character(len=SLEN) :: parent, submodule_name

          submodule_name = adjustl(line(10:len_trim(line)))

          ! get parent name
          if (submodule_name(1:1) /= '(') call error("submodule specification error")
          i = scan(submodule_name, ')')
          if (i <= 2) call error("submodule specification error")
          parent = trim(adjustl(submodule_name(2:i-1)))
          if (parent == "") call error("submodule specification error")

          ! get submodule name
          if (i == len_trim(submodule_name)) call error("submodule specification error")
          submodule_name = adjustl(submodule_name(i+1:len_trim(submodule_name)))
          if (submodule_name == "") call error("submodule specification error")

          ! use parent module
          call push_unique(this%use_modules, to_lower(parent))

          ! add submodule
          call push_unique(this%submodules, trim(to_lower(parent))//"@"//trim(to_lower(submodule_name)))
        end block
      end if

      ! check for use statement
      if (to_lower(line(1:3)) == "use") then
        block
          character(len=SLEN) :: module_name
          integer             :: i0

          if (len_trim(line) < 5) call error("use statement error")
          module_name = adjustl(line(4:len_trim(line)))
          if (module_name == "") call error("use statement error")

          ! search for ::
          i0 = scan(module_name, ':')
          if (i0 /= 0) then
            if (module_name(i0:i0+1) == "::") then
              if (len_trim(module_name) < i0+2) call error("use statement error")
              module_name = adjustl(module_name(i0+2:len_trim(module_name)))
            end if
          end if

          ! remove stuff after ,
          i0 = scan(module_name, ',')
          if (i0 /= 0) then
            if (i0 < 2) call error("use statement error")
            module_name = trim(module_name(1:i0-1))
          end if

          if (module_name == "") call error("use statement error")
          call push_unique(this%use_modules, to_lower(module_name))
        end block
      end if

      ! check for includes
      block
        character(len=SLEN) :: inc
        character(len=1)    :: quote

        if (to_lower(line(1:7)) == "include") then
          inc = adjustl(line(8:len_trim(line)))
        elseif (to_lower(line(1:8)) == "#include") then
          inc = adjustl(line(9:len_trim(line)))
        else
          inc = ""
        end if
        if (inc /= "") then
          if (len_trim(inc) < 3) call error("include parse error")
          quote = inc(1:1)
          if ((quote /= "'") .and. (quote /= '"')) call error("include parse error")
          if (inc(len_trim(inc):len_trim(inc)) /= quote) call error("include parse error")
          inc = trim(adjustl(inc(2:len_trim(inc)-1)))
          if (inc == "") call error("include parse error")
          call push_unique(this%includes, trim(inc))
        end if
      end block
    end do

    allocate (this%dep_modules( this%use_modules%n), source = 0)
  end subroutine

  function to_lower(str) result(lstr)
    !! covert string to lower case

    character(len=*), intent(in) :: str
    character(len=len(str))      :: lstr

    ! local variables
    integer :: i

    lstr = ""
    do i = 1, len_trim(str)
      select case (str(i:i))
      case ("A":"Z")
        lstr(i:i) = achar(iachar(str(i:i))+32)
      case default
        lstr(i:i) = str(i:i)
      end select
    end do
  end function

  function to_upper(str) result(ustr)
    !! covert string to lower case

    character(len=*), intent(in) :: str
    character(len=len(str))      :: ustr

    ! local variables
    integer :: i

    ustr = ""
    do i = 1, len_trim(str)
      select case (str(i:i))
      case ("a":"z")
        ustr(i:i) = achar(iachar(str(i:i))-32)
      case default
        ustr(i:i) = str(i:i)
      end select
    end do
  end function

  subroutine push_unique(str_list, str, idx)
    !! add string to vector if not already included

    type(vector_string), intent(inout) :: str_list
    character(len=*),    intent(in)    :: str
    integer, optional,   intent(out)   :: idx

    ! local variables
    integer             :: i
    character(len=SLEN) :: str_tmp

    str_tmp = str

    do i = 1, str_list%n
      if (str_list%d(i) == str) then
        if (present(idx)) idx = i
        return
      end if
    end do
    call str_list%push(str_tmp)
    if (present(idx)) idx = str_list%n
  end subroutine

  subroutine error(msg)
    character(len=*), intent(in) :: msg

    write(error_unit,*) msg
    error stop
  end subroutine

end module
