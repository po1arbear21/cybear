module test_plotmtv_m

  use plotmtv_m,    only: curve_options, plotmtv, plotmtv_write, plotset_options, view3d_options
  use string_m,     only: string
  use test_case_m,  only: test_case

  implicit none

  private
  public test_plotmtv

contains

  subroutine test_plotmtv()
    type(test_case) :: tc

    call tc%init("plotmtv")

    ! test1: simple 2D write (no options altered)
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)

      allocate (fname, source = "/tmp/plotmtv_test.asc")
      ! fname = "/tmp/plotmtv_test.asc"
      x     = [-1,   0,   2,  4]
      y     = [-5, -10, -15, 20]

      call plotmtv_write(fname, x, y)

      ! compare data
      block
        type(string) :: header(2)

        header(1)%s = "$ DATA=CURVE2D"
        header(2)%s = ""

        call compare_file(tc, "test1", fname, header, x, y)
      end block

      call execute_command_line("rm -f " // fname)
    end block

    ! test2: simple 2D write and plotset options altered
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)
      type(plotset_options)     :: opts

      allocate (fname, source = "/tmp/plotmtv_test.asc")
      x     = [-1,   0,   2,  4]
      y     = [-5, -10, -15, 20]

      opts%toplabel = "my plot"
      opts%xlabel   = "my x data"
      opts%ylabel   = "my y data"
      opts%grid     = .true.

      call plotmtv_write(fname, x, y, plotset_opts=opts)

      ! compare data
      block
        type(string) :: header(6)

        header(1)%s = "$ DATA=CURVE2D"
        header(2)%s = ""
        header(3)%s = '% xlabel = "my x data"'
        header(4)%s = '% ylabel = "my y data"'
        header(5)%s = '% toplabel = "my plot"'
        header(6)%s = "% grid = True"

        call compare_file(tc, "test2", fname, header, x, y)
      end block

      call execute_command_line("rm -f " // fname)
    end block

    ! test3: write and plotset+curve options altered
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:)
      type(plotset_options)     :: ps_opts
      type(curve_options)       :: c_opts

      allocate (fname, source = "/tmp/plotmtv_test.asc")
      x     = [-1,   0,   2,  4]
      y     = [-5, -10, -15, 20]

      ! plotset options
      ps_opts%toplabel = "my plot"
      ps_opts%xlabel   = "my x data"
      ps_opts%ylabel   = "my y data"
      ps_opts%grid     = .true.

      ! curve options
      c_opts%linetype   = 4
      c_opts%markertype = 5
      c_opts%markersize = 3
      c_opts%filltype   = 6

      call plotmtv_write(fname, x, y, plotset_opts=ps_opts, curve_opts=c_opts)

      ! compare data
      block
        character(10) :: c
        type(string)  :: header(11)

        header( 1)%s = "$ DATA=CURVE2D"
        header( 2)%s = ""
        header( 3)%s = '% xlabel = "my x data"'
        header( 4)%s = '% ylabel = "my y data"'
        header( 5)%s = '% toplabel = "my plot"'
        header( 6)%s = "% grid = True"
        header( 7)%s = ""
        write (c, "(I0)") 4
        header( 8)%s = "% linetype = "   // c
        write (c, "(I0)") 5
        header( 9)%s = "% markertype = " // c
        write (c, "(I0)") 3
        header(10)%s = "% markersize = " // c
        write (c, "(I0)") 6
        header(11)%s = "% filltype = "   // c

        call compare_file(tc, "test3", fname, header, x, y)
      end block

      call execute_command_line("rm -f " // fname)
    end block

    ! test4: write3d. create pyramid (possible for grid visualization)
    block
      character(:), allocatable :: fname
      real,         allocatable :: x(:), y(:), z(:)
      real,         allocatable :: xyz_tot(:,:,:)
      type(plotmtv)             :: pmtv
      type(plotset_options)     :: ps_opts
      type(curve_options)       :: c_opts, gl_c_opts, c_opts_tot(4)
      type(view3d_options)      :: v3_opts

      ! save all points in global xyz_tot array
      allocate (xyz_tot(4,3,4))

      ! plotset options
      ps_opts%toplabel = "my plot"
      ps_opts%xlabel   = "my x data"
      ps_opts%ylabel   = "my y data"
      ps_opts%grid     = .true.

      ! global curve options
      gl_c_opts%linetype   = 4
      gl_c_opts%markertype = 5
      gl_c_opts%markersize = 3
      gl_c_opts%filltype   = 6

      ! view3d options
      v3_opts%eyepos_x   = 2.0
      v3_opts%yaxisscale = 3.0
      v3_opts%hiddenline = .true.

      allocate (fname, source = "/tmp/plotmtv_test.asc")

      call pmtv%init(fname)
      call pmtv%write_header(three_dim=.true., plotset_opts=ps_opts, view3d_opts=v3_opts, gl_curve_opts=gl_c_opts)

      ! write pyramid side 1
      c_opts%fillcolor = 1
      x = [0, 1, 2, 0]
      y = [0, 1, 0, 0]
      z = [0, 1, 0, 0]
      call pmtv%write_curve(x, y, z=z, opts=c_opts)

      xyz_tot(:,1,1) = x; xyz_tot(:,2,1) = y; xyz_tot(:,3,1) = z
      c_opts_tot(1) = c_opts

      ! write pyramid side 2
      c_opts%fillcolor = 2
      x = [2, 1, 2, 2]
      y = [0, 1, 2, 0]
      z = [0, 1, 0, 0]
      call pmtv%write_curve(x, y, z=z, opts=c_opts)

      xyz_tot(:,1,2) = x; xyz_tot(:,2,2) = y; xyz_tot(:,3,2) = z
      c_opts_tot(2) = c_opts

      ! write pyramid side 3
      c_opts%fillcolor = 3
      x = [2, 1, 0, 2]
      y = [2, 1, 2, 2]
      z = [0, 1, 0, 0]
      call pmtv%write_curve(x, y, z=z, opts=c_opts)

      xyz_tot(:,1,3) = x; xyz_tot(:,2,3) = y; xyz_tot(:,3,3) = z
      c_opts_tot(3) = c_opts

      ! write pyramid side 4
      c_opts%fillcolor = 4
      x = [0, 1, 0, 0]
      y = [2, 1, 0, 2]
      z = [0, 1, 0, 0]
      call pmtv%write_curve(x, y, z=z, opts=c_opts)

      xyz_tot(:,1,4) = x; xyz_tot(:,2,4) = y; xyz_tot(:,3,4) = z
      c_opts_tot(4) = c_opts

      call pmtv%close

      ! compare data
      block
        character(10) :: c_int
        character(32) :: c_real
        type(string)  :: header(15)

        header( 1)%s = "$ DATA=CURVE3D"
        header( 2)%s = ""
        header( 3)%s = '% xlabel = "my x data"'
        header( 4)%s = '% ylabel = "my y data"'
        header( 5)%s = '% toplabel = "my plot"'
        header( 6)%s = "% grid = True"
        header( 7)%s = ""
        write (c_real, "(ES25.16E3)") 2.0
        header( 8)%s = "% eyepos.x = "   // c_real
        write (c_real, "(ES25.16E3)") 3.0
        header( 9)%s = "% yaxisscale = " // c_real
        header(10)%s = "% hiddenline = True"
        header(11)%s = ""
        write (c_int, "(I0)") 4
        header(12)%s = "% dlinetype = "   // c_int
        write (c_int, "(I0)") 5
        header(13)%s = "% dmarkertype = " // c_int
        write (c_int, "(I0)") 3
        header(14)%s = "% dmarkersize = " // c_int
        write (c_int, "(I0)") 6
        header(15)%s = "% dfilltype = "   // c_int

        call compare_file_3d(tc, "test4", fname, header, xyz_tot, c_opts_tot)
      end block

      call execute_command_line("rm -f " // fname)
    end block

    call tc%finish
  end subroutine

  subroutine compare_file(tc, tc_msg, fname, header, x, y)
    !! checks if in file "fname" the following holds:
    !!    1) first lines given by header
    !!    2) after that lines represent [x, y]

    type(test_case), intent(inout) :: tc
    character(*),    intent(in)    :: tc_msg
    character(*),    intent(in)    :: fname
    type(string),    intent(in)    :: header(:)
    real,            intent(in)    :: x(:), y(:)

    character(4096) :: char_line
    integer         :: iounit, ios, i
    real            :: xi, yi
    type(string)    :: str_line

    open (newunit=iounit, file=fname, iostat=ios, status="old", action="read")
    call tc%assert_eq(0, ios, "opening file")

    do i = 1, size(header)
      read (iounit, "(A)") char_line

      str_line%s = trim(char_line)
      call tc%assert_eq(header(i), str_line, tc_msg // ": header wrong")
    end do

    do i = 1, size(x)
      read (iounit, *) xi, yi

      call tc%assert_eq(x(i), xi, 1e-14, epsilon(1.0), tc_msg // ": 1st col wrong")
      call tc%assert_eq(y(i), yi, 1e-14, epsilon(1.0), tc_msg // ": 2nd col wrong")
    end do

    close (unit=iounit, iostat=ios)
    call tc%assert_eq(0, ios, "closing file")
  end subroutine

  subroutine compare_file_3d(tc, tc_msg, fname, header, xyz, c_opts)
    !! checks if in file "fname" the following holds:
    !!    1) first lines given by header
    !!    2) after that lines represent [x, y, z]

    type(test_case),     intent(inout) :: tc
    character(*),        intent(in)    :: tc_msg
    character(*),        intent(in)    :: fname
    type(string),        intent(in)    :: header(:)
    real,                intent(in)    :: xyz(:,:,:)
      !! indices(i,ix,ic)
      !! values represent:
      !!    point i
      !!    coordinate value x(ix=1), y(ix=2), z(ix=3)
      !!    curve ic
    type(curve_options), intent(in)    :: c_opts(:)
      !! index(ic)
      !! curve option for curve ic

    character(4096) :: char_line
    character(10)   :: c_int
    integer         :: iounit, ios, i, ic
    real            :: xi, yi, zi
    type(string)    :: str_line, new_line, fillcolor_line_exp, fillcolor_line_val

    open (newunit=iounit, file=fname, iostat=ios, action="read")
    call tc%assert_eq(0, ios, "opening file")

    do i = 1, size(header)
      read (iounit, "(A)") char_line
      str_line%s = trim(char_line)
      call tc%assert_eq(header(i), str_line, tc_msg // ": header wrong")
    end do

    call tc%assert_eq(size(c_opts), size(xyz, dim=3), tc_msg // ": wrong number of curves")

    new_line%s = ""

    do ic = 1, size(c_opts)
      ! always starts with a newline
      read (iounit, "(A)") char_line
      str_line%s = trim(char_line)
      call tc%assert_eq(new_line, str_line, tc_msg // ": curve's newline missing")

      ! check that curve option (fillcolor) is correct
      write (c_int, "(I0)") ic
      fillcolor_line_exp%s = "% fillcolor = "   // c_int
      read (iounit, "(A)") char_line
      fillcolor_line_val%s = trim(char_line)
      call tc%assert_eq(fillcolor_line_exp, fillcolor_line_val, tc_msg // ": curve's fillcolor")

      ! check node values
      do i = 1, size(xyz, dim=1)
        read (iounit, *) xi, yi, zi

        call tc%assert_eq(xyz(i,1,ic), xi, 1e-14, epsilon(1.0), tc_msg // ": 1st col wrong")
        call tc%assert_eq(xyz(i,2,ic), yi, 1e-14, epsilon(1.0), tc_msg // ": 2nd col wrong")
        call tc%assert_eq(xyz(i,3,ic), zi, 1e-14, epsilon(1.0), tc_msg // ": 3rd col wrong")
      end do
    end do

    close (unit=iounit, iostat=ios)
    call tc%assert_eq(0, ios, "closing file")
  end subroutine

end module
