m4_include(macro.f90.inc)

module cl_options_m

  use error_m,  only: assert_failed, program_error
  use map_m,    only: map_string_int, mapnode_string_int
  use string_m, only: new_string, string
  use util_m,   only: int2str
  use vector_m, only: vector_int

  implicit none

  private
  public cl_option_descriptor, new_cl_option_descriptor
  public cl_option, get_cl_options, get_cl_options_simple

  integer, parameter :: MAX_LONG_NAME = 256

  type cl_option_descriptor
    character(1)              :: short
      !! short name (-m)
    character(MAX_LONG_NAME)  :: long
      !! long name (--my_argument)
    logical                   :: req
      !! option required ?
    logical                   :: multi
      !! multiple options with this name allowed ?
    logical                   :: arg
      !! option has argument ?
    logical                   :: arg_req
      !! argument is required ?
  end type

  type cl_option
    character(1)              :: short
      !! short name
    character(:), allocatable :: arg
      !! argument
  end type

contains

  function new_cl_option_descriptor(short, long, arg, arg_req) result(desc)
    !! create new command line option descriptor
    character(1)               :: short
      !! short name (-m)
    character(:), allocatable  :: long
      !! long name (--my_argument)
    logical                    :: arg
      !! option has argument
    logical                    :: arg_req
      !! argument is required
    type(cl_option_descriptor) :: desc
      !! return descriptor

    desc%short   = short
    desc%long    = long
    desc%arg     = arg
    desc%arg_req = arg_req
  end function

  subroutine get_cl_options(desc, opt, iopt, jopt)
    type(cl_option_descriptor),   intent(in)  :: desc(:)
      !! list of command line descriptors
    type(cl_option), allocatable, intent(out) :: opt(:)
      !! output command line options
    integer,         allocatable, intent(out) :: iopt(:)
      !! output indices of command line options sorted by descriptors (size(iopt) = size(desc) + 1)
      !! options for idesc-th descriptor: opt(jopt(iopt(idesc):iopt(idesc+1)-1))
    integer,         allocatable, intent(out) :: jopt(:)
      !! output indices of command line options sorted by descriptors (size(iopt) = size(desc) + 1)
      !! options for idesc-th descriptor: opt( jopt( iopt(idesc):iopt(idesc+1)-1 ) )

    integer                           :: i, idesc, iarg, length, nargs, nopt
    integer,      allocatable         :: kopt(:)
    type(map_string_int)              :: map_long, map_short
    type(mapnode_string_int), pointer :: node
    type(string)                      :: str
    type(string), allocatable         :: args(:)
    type(vector_int)                  :: vparse

    ! create maps
    call map_short%init()
    call map_long%init()
    do idesc = 1, size(desc)
      str = new_string(desc(idesc)%short)
      if (associated(map_short%find(str))) call program_error("Found multiple descriptors with the same short option name")
      call map_short%set(str, idesc)

      str = new_string(trim(desc(idesc)%long))
      if (associated(map_long%find(str))) call program_error("Found multiple descriptors with the same long option name")
      call map_long%set(str, idesc)
    end do

    ! get raw command line arguments
    args  = get_cl_args()
    nargs = size(args)

    ! parse options
    idesc = 0
    nopt = 0
    call vparse%init(0, c = nargs)
    do iarg = 1, nargs
      length = len(args(iarg)%s)

      ! handle optional arguments
      if (idesc > 0) then
        if (args(iarg)%s(1:1) == "-") then
          if (desc(idesc)%arg_req) call program_error("Command line argument '" // trim(desc(idesc)%long) // "' has required option")
          idesc = 0
        end if
      end if

      if (idesc == 0) then ! expect option name
        if (length < 2) call program_error("Command line argument must have length >= 2")
        if (args(iarg)%s(1:2) == "--") then ! long name
          if (length < 3) call program_error("Command line argument with long name expected")
          str = new_string(args(iarg)%s(3:length))
          node => map_long%find(str)
          if (.not. associated(node)) call program_error("Unknown command line argument with long name '" // str%s // "'")
        elseif (args(iarg)%s(1:1) == "-") then ! short name
          if (length > 2) call program_error("Command line argument with short name must have length 2 (including hyphen)")
          str = new_string(args(iarg)%s(2:2))
          node => map_short%find(str)
          if (.not. associated(node)) call program_error("Unknown command line argument with short name '" // str%s // "'")
        else
          call program_error("Command line argument must start with - or --")
        end if

        ! new command line option
        nopt = nopt + 1

        ! two integers per option: descriptor index and iarg for option (actual value is set later)
        idesc = node%value
        call vparse%push(idesc)
        call vparse%push(0)

        ! reset descriptor index if no argument expected
        if (.not. desc(idesc)%arg) idesc = 0
      else ! option argument
        vparse%d(vparse%n) = iarg
        idesc = 0
      end if
    end do
    if (idesc > 0) then
      if (desc(idesc)%arg_req) call program_error("Command line argument '" // trim(desc(idesc)%long) // "' has required option")
    end if

    ! return
    allocate (opt(nopt))
    allocate (iopt(size(desc) + 1), jopt(nopt), kopt(size(desc)), source = 0)
    do i = 1, nopt
      idesc = vparse%d(2 * i - 1)
      iarg  = vparse%d(2 * i    )
      opt(i)%short = desc(idesc)%short
      if (desc(idesc)%arg) then
        if (iarg > 0) opt(i)%arg = args(iarg)%s
      end if
      iopt(idesc + 1) = iopt(idesc + 1) + 1
    end do
    iopt(1) = 1
    do idesc = 2, size(desc)+1
      if (desc(idesc-1)%req .and. (iopt(idesc) < 1)) call program_error("Missing required command line option '" // trim(desc(idesc-1)%long) // "'")
      if ((.not. desc(idesc-1)%multi) .and. (iopt(idesc) > 1)) call program_error("Found multiple instances of singular command line option '" // trim(desc(idesc-1)%long) //"'")
      iopt(idesc) = iopt(idesc) + iopt(idesc - 1)
    end do
    do i = 1, nopt
      idesc = vparse%d(2 * i - 1)
      jopt(iopt(idesc) + kopt(idesc)) = i
      kopt(idesc) = kopt(idesc) + 1
    end do
  end subroutine

  subroutine get_cl_options_simple(names, values)
    !! provide long names, return corresponding values, without need for descriptors
    type(string), intent(in)  :: names(:)
      !! long names (short = first letter, req = true, multi = false, arg = true, arg_req = true)
    type(string), intent(out) :: values(:)
      !! return values

    character(1)                            :: short
    type(string)                            :: s
    integer                                 :: i, ii, j, k, n
    type(cl_option_descriptor), allocatable :: desc(:)
    type(cl_option),            allocatable :: opt(:)
    type(map_string_int)                    :: short_map
    type(mapnode_string_int),   pointer     :: node
    integer,                    allocatable :: iopt(:), jopt(:)

    n = size(names)
    m4_assert(size(values) == n)

    ! create descriptors
    allocate (desc(n))
    call short_map%init()
    do i = 1, n
      short = names(i)%s(1:1)
      s = new_string(short)
      node => short_map%find(s)
      if (associated(node)) then
        call program_error("Found multiple option names starting with letter '" // short // "'")
      end if
      desc(i) = cl_option_descriptor(names(i)%s(1:1), names(i)%s, .true., .false., .true., .true.)
      call short_map%insert(s, i)
    end do

    ! get command line options and collect corresponding values
    call get_cl_options(desc, opt, iopt, jopt)
    do i = 1, n
      do ii = iopt(i), iopt(i + 1) - 1
        j = jopt(ii)
        k = short_map%get(new_string(opt(j)%short))
        values(k)%s = opt(j)%arg
      end do
    end do
  end subroutine

  function get_cl_arg(iarg) result(arg)
    !! retrieve command line argument
    integer, intent(in)       :: iarg
      !! get iarg-th argument
    character(:), allocatable :: arg
      !! return argument

    integer :: len, stat

    call get_command_argument(number = iarg, length = len, status = stat)
    if (stat > 0) call program_error("Can not retrieve argument number " // int2str(iarg))
    allocate (character(len) :: arg)
    call get_command_argument(number = iarg, value = arg, status = stat)
  end function

  function get_cl_args() result(args)
    !! get all command line arguments
    type(string), allocatable :: args(:)

    integer :: iarg, nargs

    nargs = command_argument_count()

    allocate (args(nargs))
    do iarg = 1, nargs
      args(iarg)%s = get_cl_arg(iarg)
    end do
  end function

end module
