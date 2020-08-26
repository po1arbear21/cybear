module test_plotmtv_m

  use plotmtv_m
  use test_case_m
  use string_m

  implicit none

contains

  subroutine test_plotmtv()
    type(test_case)    :: tc

    print "(A)", "test_plotmtv"
    call tc%init("plotmtv")

    ! test1: simple write (no options altered)
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)

      allocate (character(0) :: fname)      ! removes gfortran warning
      fname = '/tmp/test_plotmtv.asc'

      x = [-1,   0,   2,  4]
      y = [-5, -10, -15, 20]

      call write_plotmtv(fname, x, y)

      ! compare data
      block
        type(string) :: header(2)

        header(1)%s = '$ DATA=CURVE2D'
        header(2)%s = ''

        call compare_file(tc, "test1", fname, header, x, y)
      end block
    end block

    ! test2: write and plotset options altered
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)
      type(plotmtv_opts)        :: opts

      allocate (character(0) :: fname)      ! removes gfortran warning
      fname = '/tmp/test_plotmtv.asc'

      x = [-1,   0,   2,  4]
      y = [-5, -10, -15, 20]

      opts%ps%toplabel = "my plot"
      opts%ps%xlabel   = "my x data"
      opts%ps%ylabel   = "my y data"
      opts%ps%grid     = .true.

      call write_plotmtv(fname, x, y, opts=opts)

      ! compare data
      block
        type(string) :: header(6)

        header(1)%s = '$ DATA=CURVE2D'
        header(2)%s = ''
        header(3)%s = '% xlabel = "my x data"'
        header(4)%s = '% ylabel = "my y data"'
        header(5)%s = '% toplabel = "my plot"'
        header(6)%s = '% grid = True'

        call compare_file(tc, "test2", fname, header, x, y)
      end block
    end block

    ! test3: write and plotset+curve options altered
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)
      type(plotmtv_opts)        :: opts

      allocate (character(0) :: fname)      ! removes gfortran warning
      fname = '/tmp/test_plotmtv.asc'

      x = [-1,   0,   2,  4]
      y = [-5, -10, -15, 20]

      ! plotset options
      opts%ps%toplabel = "my plot"
      opts%ps%xlabel   = "my x data"
      opts%ps%ylabel   = "my y data"
      opts%ps%grid     = .true.

      ! curve options
      opts%c%linetype   = 4
      opts%c%markertype = 5
      opts%c%markersize = 3
      opts%c%filltype   = 6

      call write_plotmtv(fname, x, y, opts=opts)

      ! compare data
      block
        character(10) :: c
        type(string)  :: header(10)

        header( 1)%s = '$ DATA=CURVE2D'
        header( 2)%s = ''
        header( 3)%s = '% xlabel = "my x data"'
        header( 4)%s = '% ylabel = "my y data"'
        header( 5)%s = '% toplabel = "my plot"'
        header( 6)%s = '% grid = True'
        write (c, '(I10)') 4
        header( 7)%s = '% linetype = '   // c
        write (c, '(I10)') 5
        header( 8)%s = '% markertype = ' // c
        write (c, '(I10)') 3
        header( 9)%s = '% markersize = ' // c
        write (c, '(I10)') 6
        header(10)%s = '% filltype = '   // c

        call compare_file(tc, "test3", fname, header, x, y)
      end block
    end block

    ! test4: write3d and plotset+curve+view3d options altered
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:), z(:)
      type(plotmtv_opts)        :: opts

      allocate (character(0) :: fname)      ! removes gfortran warning
      fname = '/tmp/test_plotmtv.asc'

      x = [-1,   0,   2,  4]
      y = [-5, -10, -15, 20]
      z = [-6, -11, -16, 21]

      ! plotset options
      opts%ps%toplabel = "my plot"
      opts%ps%xlabel   = "my x data"
      opts%ps%ylabel   = "my y data"
      opts%ps%grid     = .true.

      ! curve options
      opts%c%linetype   = 4
      opts%c%markertype = 5
      opts%c%markersize = 3
      opts%c%filltype   = 6

      ! view3d options
      opts%v3%eyepos_x   = 2.0
      opts%v3%yaxisscale = 3.0
      opts%v3%hiddenline = .true.

      call write_plotmtv(fname, x, y, z=z, opts=opts)

      ! compare data
      block
        character(10) :: c_int
        character(15) :: c_real
        type(string)  :: header(13)

        header( 1)%s = '$ DATA=CURVE3D'
        header( 2)%s = ''
        header( 3)%s = '% xlabel = "my x data"'
        header( 4)%s = '% ylabel = "my y data"'
        header( 5)%s = '% toplabel = "my plot"'
        header( 6)%s = '% grid = True'
        write (c_int, '(I10)') 4
        header( 7)%s = '% linetype = '   // c_int
        write (c_int, '(I10)') 5
        header( 8)%s = '% markertype = ' // c_int
        write (c_int, '(I10)') 3
        header( 9)%s = '% markersize = ' // c_int
        write (c_int, '(I10)') 6
        header(10)%s = '% filltype = '   // c_int
        write (c_real, '(E15.5)') 2.0
        header(11)%s = '% eyepos.x = '   // c_real
        write (c_real, '(E15.5)') 3.0
        header(12)%s = '% yaxisscale = ' // c_real
        header(13)%s = '% hiddenline = True'

        call compare_file(tc, "test4", fname, header, x, y, z=z)
      end block
    end block

    call tc%finish
  end subroutine

  subroutine compare_file(tc, tc_msg, fname, header, x, y, z)
    !! checks if in file "fname" the following holds:
    !!    1) first lines given by header
    !!    2) after that lines represent [x, y] matrix (or [x, y, z])

    type(test_case), intent(inout)        :: tc
    character(*),    intent(in)           :: tc_msg
    character(*),    intent(in)           :: fname
    type(string),    intent(in)           :: header(:)
    real,            intent(in)           :: x(:), y(:)
    real,            intent(in), optional :: z(:)

    character(4096) :: char_line
    integer         :: iounit, ios, i
    real            :: xi, yi, zi
    type(string)    :: str_line

    open (newunit=iounit, file=fname, iostat=ios, action='read')
    call tc%assert_eq(0, ios, "opening file")

    do i = 1, size(header)
      read (iounit, "(A)") char_line
      str_line%s = trim(char_line)
      call tc%assert_eq(header(i), str_line, tc_msg // ": header wrong")
    end do

    do i = 1, size(x)
      if (present(z)) then
        read (iounit, *) xi, yi, zi
      else
        read (iounit, *) xi, yi
      end if

      call tc%assert_eq(                x(i), xi, epsilon(1.0), tc_msg // ": 1st col wrong")
      call tc%assert_eq(                y(i), yi, epsilon(1.0), tc_msg // ": 2nd col wrong")
      if (present(z)) call tc%assert_eq(z(i), zi, epsilon(1.0), tc_msg // ": 3rd col wrong")
    end do

    close (unit=iounit, iostat=ios)
    call tc%assert_eq(0, ios, "closing file")
  end subroutine

end module
