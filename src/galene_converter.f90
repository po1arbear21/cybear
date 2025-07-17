program galene_converter
  !! GALENE data reader and converter
  !! Converts binary .OSV and .GEO files to readable formats
  !! Supports ASCII output, MATLAB .dat format, and TikZ-compatible data

  use galene_m
  use normalization_m, only: init_normconst, denorm
  use iso_fortran_env, only: output_unit

  implicit none

  type(gal_file)           :: gat, sil
  type(gal_block), pointer :: b
  integer                  :: i, j, choice
  character(256)           :: filename, variable_name
  character(256)           :: osv_filename, geo_filename
  logical                  :: file_exists

  ! Initialize normalization constants
  call init_normconst(300.0)

  ! Default file names (can be changed here)
  osv_filename = "DD.OSV"
  geo_filename = "DGMOS.GEO"

  ! Check if input files exist
  inquire(file=trim(osv_filename), exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "Error: ", trim(osv_filename), " file not found!"
    stop
  end if

  inquire(file=trim(geo_filename), exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "Error: ", trim(geo_filename), " file not found!"
    stop
  end if

  ! Read binary files
  write(*,*) "Reading binary data files..."
  call gat%init(trim(osv_filename))
  call sil%init(trim(geo_filename))
  write(*,*) "Files loaded successfully!"

  ! Main menu
  do while (.true.)
    write(*,*)
    write(*,*) "===== GALENE Data Converter ====="
    write(*,*) "1. List all variables in ", trim(osv_filename)
    write(*,*) "2. List all variables in ", trim(geo_filename)
    write(*,*) "3. Export variable to ASCII file"
    write(*,*) "4. Export variable to MATLAB .m format"
    write(*,*) "5. Export variable to MATLAB .mat format"
    write(*,*) "6. Export all variables to .mat database"
    write(*,*) "7. Export to standalone TikZ figure"
    write(*,*) "8. Export all variables summary"
    write(*,*) "9. Quick view variable"
    write(*,*) "0. Exit"
    write(*,*) "Enter your choice: "
    read(*,*) choice

    select case(choice)
    case(0)
      exit
    case(1)
      call list_variables(gat, trim(osv_filename))
    case(2)
      call list_variables(sil, trim(geo_filename))
    case(3)
      call export_to_ascii(gat, sil)
    case(4)
      call export_to_matlab(gat, sil)
    case(5)
      call export_to_mat_file(gat, sil)
    case(6)
      call export_to_mat_database(gat, sil)
    case(7)
      call export_to_tikz_standalone(gat, sil)
    case(8)
      call export_summary(gat, sil, osv_filename, geo_filename)
    case(9)
      call quick_view(gat, sil)
    case default
      write(*,*) "Invalid choice!"
    end select
  end do

contains

  subroutine list_variables(gf, filename)
    !! List all variables in a GALENE file
    type(gal_file), intent(in) :: gf
    character(*), intent(in) :: filename
    integer :: i

    write(*,*)
    write(*,*) "Variables in ", trim(filename), ":"
    write(*,*) "==============================="

    do i = 1, gf%blocks%n
      associate(block => gf%blocks%d(i))
        write(*,'(I3,". ",A," [",A,"] (",I0," elements)")') &
              i, trim(block%name), trim(block%unit), block%ndata
      end associate
    end do
  end subroutine

  subroutine export_to_ascii(gat, sil)
    !! Export a variable to ASCII file
    type(gal_file), intent(in) :: gat, sil
    character(256) :: varname, outfile
    type(gal_block), pointer :: block
    integer :: unit, i

    write(*,*) "Enter variable name: "
    read(*,*) varname

    ! Try to find in simulation data first
    block => gat%get_block(trim(varname))
    if (.not. associated(block)) then
      ! Try geometry data
      block => sil%get_block(trim(varname))
    end if

    if (.not. associated(block)) then
      write(*,*) "Variable '", trim(varname), "' not found!"
      return
    end if

    outfile = trim(varname) // ".dat"
    open(newunit=unit, file=outfile, status='replace')

    write(unit,'("# Variable: ",A)') trim(block%name)
    write(unit,'("# Unit: ",A)') trim(block%unit)
    write(unit,'("# Number of elements: ",I0)') block%ndata
    write(unit,'("# Data type: ",I0)') block%dtype
    write(unit,'("#")')

    select case(block%dtype)
    case(GALDATA_INT)
      do i = 1, block%ndata
        write(unit,'(I0)') block%idata(i)
      end do
    case(GALDATA_REAL)
      do i = 1, block%ndata
        write(unit,'(ES15.8)') denorm(block%rdata(i), trim(block%unit))
      end do
    case(GALDATA_CHAR)
      write(unit,'(A)') trim(block%cdata)
    end select

    close(unit)
    write(*,*) "Data exported to: ", trim(outfile)
  end subroutine

  subroutine export_to_matlab(gat, sil)
    !! Export a variable to MATLAB-compatible format
    type(gal_file), intent(in) :: gat, sil
    character(256) :: varname, outfile, matlab_varname
    type(gal_block), pointer :: block
    integer :: unit, i

    write(*,*) "Enter variable name: "
    read(*,*) varname

    ! Try to find in simulation data first
    block => gat%get_block(trim(varname))
    if (.not. associated(block)) then
      ! Try geometry data
      block => sil%get_block(trim(varname))
    end if

    if (.not. associated(block)) then
      write(*,*) "Variable '", trim(varname), "' not found!"
      return
    end if

    ! Create MATLAB-compatible variable name
    matlab_varname = trim(varname)
    ! Replace hyphens with underscores for MATLAB compatibility
    do i = 1, len_trim(matlab_varname)
      if (matlab_varname(i:i) == '-') matlab_varname(i:i) = '_'
    end do

    outfile = trim(matlab_varname) // ".m"
    open(newunit=unit, file=outfile, status='replace')

    write(unit,'("% Variable: ",A)') trim(block%name)
    write(unit,'("% Unit: ",A)') trim(block%unit)
    write(unit,'("% Number of elements: ",I0)') block%ndata
    write(unit,'("%")')

    select case(block%dtype)
    case(GALDATA_INT)
      write(unit,'(A," = [")') trim(matlab_varname)
      do i = 1, block%ndata
        if (i < block%ndata) then
          write(unit,'(I0,";")') block%idata(i)
        else
          write(unit,'(I0)') block%idata(i)
        end if
      end do
      write(unit,'("];")')
    case(GALDATA_REAL)
      write(unit,'(A," = [")') trim(matlab_varname)
      do i = 1, block%ndata
        if (i < block%ndata) then
          write(unit,'(ES15.8,";")') denorm(block%rdata(i), trim(block%unit))
        else
          write(unit,'(ES15.8)') denorm(block%rdata(i), trim(block%unit))
        end if
      end do
      write(unit,'("];")')
    case(GALDATA_CHAR)
      write(unit,'(A," = ''",A,"'';")') trim(matlab_varname), trim(block%cdata)
    end select

    close(unit)
    write(*,*) "MATLAB data exported to: ", trim(outfile)
  end subroutine

  subroutine export_to_tikz(sil)
    !! Export geometry data to TikZ-compatible format
    type(gal_file), intent(in) :: sil
    type(gal_block), pointer :: x_block, y_block
    integer :: unit, i

    x_block => sil%get_block("x-coord")
    y_block => sil%get_block("y-coord")

    if (.not. associated(x_block) .or. .not. associated(y_block)) then
      write(*,*) "Error: x-coord or y-coord not found in geometry file!"
      return
    end if

    if (x_block%ndata /= y_block%ndata) then
      write(*,*) "Error: x and y coordinates have different sizes!"
      return
    end if

    open(newunit=unit, file="geometry.tikz", status='replace')

    write(unit,'("% TikZ geometry data")')
    write(unit,'("% X-coordinates in ",A)') trim(x_block%unit)
    write(unit,'("% Y-coordinates in ",A)') trim(y_block%unit)
    write(unit,'("% Number of points: ",I0)') x_block%ndata
    write(unit,'("%")')
    write(unit,'("\\coordinate (device) at (0,0);")')
    write(unit,'("\\begin{scope}[shift={(device)}]")')

    do i = 1, x_block%ndata
      write(unit,'("\\coordinate (p",I0,") at (",F12.6,",",F12.6,");")') &
            i, x_block%rdata(i), y_block%rdata(i)
    end do

    write(unit,'("\\end{scope}")')
    close(unit)
    write(*,*) "TikZ geometry exported to: geometry.tikz"
  end subroutine

  subroutine export_summary(gat, sil, osv_filename, geo_filename)
    !! Export summary of all variables
    type(gal_file), intent(in) :: gat, sil
    character(*), intent(in) :: osv_filename, geo_filename
    integer :: unit, i

    open(newunit=unit, file="variables_summary.txt", status='replace')

    write(unit,'("GALENE Data Summary")')
    write(unit,'("===================")')
    write(unit,'()')
    write(unit,'("Simulation Data (",A,"):")') trim(osv_filename)
    write(unit,'("-------------------------")')

    do i = 1, gat%blocks%n
      associate(block => gat%blocks%d(i))
        write(unit,'("Variable: ",A)') trim(block%name)
        write(unit,'("  Unit: ",A)') trim(block%unit)
        write(unit,'("  Type: ",I0," (1=int, 2=real, 3=char)")') block%dtype
        write(unit,'("  Elements: ",I0)') block%ndata

        if (block%dtype == GALDATA_REAL .and. block%ndata > 0) then
          write(unit,'("  Min: ",ES12.5)') denorm(minval(block%rdata), trim(block%unit))
          write(unit,'("  Max: ",ES12.5)') denorm(maxval(block%rdata), trim(block%unit))
        end if
        write(unit,'()')
      end associate
    end do

    write(unit,'("Geometry Data (",A,"):")') trim(geo_filename)
    write(unit,'("--------------------------")')

    do i = 1, sil%blocks%n
      associate(block => sil%blocks%d(i))
        write(unit,'("Variable: ",A)') trim(block%name)
        write(unit,'("  Unit: ",A)') trim(block%unit)
        write(unit,'("  Type: ",I0," (1=int, 2=real, 3=char)")') block%dtype
        write(unit,'("  Elements: ",I0)') block%ndata

        if (block%dtype == GALDATA_REAL .and. block%ndata > 0) then
          write(unit,'("  Min: ",ES12.5)') denorm(minval(block%rdata), trim(block%unit))
          write(unit,'("  Max: ",ES12.5)') denorm(maxval(block%rdata), trim(block%unit))
        end if
        write(unit,'()')
      end associate
    end do

    close(unit)
    write(*,*) "Summary exported to: variables_summary.txt"
  end subroutine

  subroutine quick_view(gat, sil)
    !! Quick view of a variable's data
    type(gal_file), intent(in) :: gat, sil
    character(256) :: varname
    type(gal_block), pointer :: block
    integer :: i, max_display

    write(*,*) "Enter variable name: "
    read(*,*) varname

    ! Try to find in simulation data first
    block => gat%get_block(trim(varname))
    if (.not. associated(block)) then
      ! Try geometry data
      block => sil%get_block(trim(varname))
    end if

    if (.not. associated(block)) then
      write(*,*) "Variable '", trim(varname), "' not found!"
      return
    end if


    write(*,*)
    write(*,*) "Variable: ", trim(block%name)
    write(*,*) "Unit: ", trim(block%unit)
    write(*,*) "Type: ", block%dtype, " (1=int, 2=real, 3=char)"
    write(*,*) "Elements: ", block%ndata
    write(*,*)

    max_display = min(10, block%ndata)

    select case(block%dtype)
    case(GALDATA_INT)
      write(*,*) "First", max_display, "values:"
      do i = 1, max_display
        write(*,'("  ",I0,": ",I0)') i, block%idata(i)
      end do
      if (block%ndata > max_display) then
        write(*,*) "  ... (", block%ndata - max_display, " more values)"
      end if
    case(GALDATA_REAL)
      write(*,*) "First", max_display, "values:"
      do i = 1, max_display
        write(*,'("  ",I0,": ",ES15.8)') i, denorm(block%rdata(i), trim(block%unit))
      end do
      if (block%ndata > max_display) then
        write(*,*) "  ... (", block%ndata - max_display, " more values)"
      end if
      if (block%ndata > 1) then
        write(*,*)
        write(*,*) "Statistics:"
        write(*,*) "  Min: ", denorm(minval(block%rdata), trim(block%unit))
        write(*,*) "  Max: ", denorm(maxval(block%rdata), trim(block%unit))
        write(*,*) "  Mean: ", denorm(sum(block%rdata) / block%ndata, trim(block%unit))
      end if
    case(GALDATA_CHAR)
      write(*,*) "Content: ", trim(block%cdata)
    end select
  end subroutine

  subroutine export_to_mat_file(gat, sil)
    !! Export a variable to binary MATLAB .mat format using Python
    type(gal_file), intent(in) :: gat, sil
    character(256) :: varname, outfile, matlab_varname
    type(gal_block), pointer :: block
    integer :: unit, i
    character(1024) :: temp_pyfile

    write(*,*) "Enter variable name: "
    read(*,*) varname

    ! Try to find in simulation data first
    block => gat%get_block(trim(varname))
    if (.not. associated(block)) then
      ! Try geometry data
      block => sil%get_block(trim(varname))
    end if

    if (.not. associated(block)) then
      write(*,*) "Variable '", trim(varname), "' not found!"
      return
    end if

    ! Create MATLAB-compatible variable name
    matlab_varname = trim(varname)
    do i = 1, len_trim(matlab_varname)
      if (matlab_varname(i:i) == '-') matlab_varname(i:i) = '_'
      if (matlab_varname(i:i) == '#') matlab_varname(i:i) = '_'
      if (matlab_varname(i:i) == '.') matlab_varname(i:i) = '_'
      if (matlab_varname(i:i) == ' ') matlab_varname(i:i) = '_'
      if (matlab_varname(i:i) == '(') matlab_varname(i:i) = '_'
      if (matlab_varname(i:i) == ')') matlab_varname(i:i) = '_'
    end do

    ! Ensure name doesn't start with underscore or digit
    if (len_trim(matlab_varname) > 0) then
      if (matlab_varname(1:1) == '_' .or. &
            (matlab_varname(1:1) >= '0' .and. matlab_varname(1:1) <= '9')) then
        matlab_varname = 'var_' // trim(matlab_varname)
      end if
    end if

    ! Truncate to 31 characters (leaving room for _info suffix)
    if (len_trim(matlab_varname) > 26) then
      matlab_varname = matlab_varname(1:26)
    end if

    ! Create temporary Python file
    temp_pyfile = 'temp_' // trim(matlab_varname) // '.py'
    open(newunit=unit, file=temp_pyfile, status='replace')

    ! Write Python script to create .mat file
    write(unit,'("import numpy as np")')
    write(unit,'("from scipy.io import savemat")')
    write(unit,'("")')

    ! Write variable data
    select case(block%dtype)
    case(GALDATA_INT)
      write(unit,'(A," = np.array([")') trim(matlab_varname)
      do i = 1, block%ndata
        if (i < block%ndata) then
          write(unit,'(I0,",")') block%idata(i)
        else
          write(unit,'(I0)') block%idata(i)
        end if
      end do
      write(unit,'("], dtype=int)")')
    case(GALDATA_REAL)
      write(unit,'(A," = np.array([")') trim(matlab_varname)
      do i = 1, block%ndata
        if (i < block%ndata) then
          write(unit,'(ES15.8,",")') denorm(block%rdata(i), trim(block%unit))
        else
          write(unit,'(ES15.8)') denorm(block%rdata(i), trim(block%unit))
        end if
      end do
      write(unit,'("], dtype=float)")')
    case(GALDATA_CHAR)
      write(unit,'(A," = ''",A,"''")') trim(matlab_varname), trim(block%cdata)
    end select

    ! Add metadata as a dictionary
    write(unit,'("")')
    write(unit,'(A,"_info = {")') trim(matlab_varname)
    write(unit,'("    ''name'': ''",A,"'',")') trim(block%name)
    write(unit,'("    ''unit'': ''",A,"'',")') trim(block%unit)
    write(unit,'("    ''dtype'': ",I0,",")') block%dtype
    write(unit,'("    ''ndata'': ",I0)') block%ndata
    write(unit,'("}")')

    outfile = trim(matlab_varname) // '.mat'
    write(unit,'("")')
    write(unit,'("savemat(''",A,"'', {''",A,"'': ",A,", ''",A,"_info'': ",A,"_info})")') &
          trim(outfile), trim(matlab_varname), trim(matlab_varname), trim(matlab_varname), trim(matlab_varname)

    close(unit)

    ! Execute Python to create .mat file
    call system('python3 ' // trim(temp_pyfile) // ' > /dev/null 2>&1')
    call system('rm -f ' // trim(temp_pyfile))

    write(*,*) "MATLAB .mat file exported to: ", trim(outfile)
  end subroutine

  subroutine export_to_mat_database(gat, sil)
    !! Export all variables to a single .mat database file using Python
    type(gal_file), intent(in) :: gat, sil
    character(256) :: temp_pyfile, matlab_varname
    integer :: unit, i, j
    type(gal_block), pointer :: block

    temp_pyfile = 'temp_database.py'
    open(newunit=unit, file=temp_pyfile, status='replace')

    write(unit,'("import numpy as np")')
    write(unit,'("from scipy.io import savemat")')
    write(unit,'("")')
    write(unit,'("# GALENE Data Database")')
    write(unit,'("data = {}")')
    write(unit,'("data[''sim''] = {}")')
    write(unit,'("data[''geo''] = {}")')
    write(unit,'("")')

    ! Export simulation data
    write(unit,'("# Simulation Data (.OSV file)")')
    do i = 1, gat%blocks%n
      associate(block => gat%blocks%d(i))
        matlab_varname = trim(block%name)
        do j = 1, len_trim(matlab_varname)
          if (matlab_varname(j:j) == '-') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '#') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '.') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == ' ') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '(') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == ')') matlab_varname(j:j) = '_'
        end do

        ! Ensure name doesn't start with underscore or digit
        if (len_trim(matlab_varname) > 0) then
          if (matlab_varname(1:1) == '_' .or. &
                (matlab_varname(1:1) >= '0' .and. matlab_varname(1:1) <= '9')) then
            matlab_varname = 'var_' // trim(matlab_varname)
          end if
        end if

        ! Truncate to 31 characters (leaving room for _info suffix)
        if (len_trim(matlab_varname) > 26) then
          matlab_varname = matlab_varname(1:26)
        end if

        select case(block%dtype)
        case(GALDATA_INT)
          write(unit,'("data[''sim''][''",A,"''] = np.array([")') trim(matlab_varname)
          do j = 1, block%ndata
            if (j < block%ndata) then
              write(unit,'(I0,",")') block%idata(j)
            else
              write(unit,'(I0)') block%idata(j)
            end if
          end do
          write(unit,'("], dtype=int)")')
        case(GALDATA_REAL)
          write(unit,'("data[''sim''][''",A,"''] = np.array([")') trim(matlab_varname)
          do j = 1, block%ndata
            if (j < block%ndata) then
              write(unit,'(ES15.8,",")') denorm(block%rdata(j), trim(block%unit))
            else
              write(unit,'(ES15.8)') denorm(block%rdata(j), trim(block%unit))
            end if
          end do
          write(unit,'("], dtype=float)")')
        case(GALDATA_CHAR)
          write(unit,'("data[''sim''][''",A,"''] = ''",A,"''")') trim(matlab_varname), trim(block%cdata)
        end select

        write(unit,'("data[''sim''][''",A,"_info''] = {")') trim(matlab_varname)
        write(unit,'("    ''name'': ''",A,"'',")') trim(block%name)
        write(unit,'("    ''unit'': ''",A,"'',")') trim(block%unit)
        write(unit,'("    ''dtype'': ",I0,",")') block%dtype
        write(unit,'("    ''ndata'': ",I0)') block%ndata
        write(unit,'("}")')
        write(unit,'("")')
      end associate
    end do

    ! Export geometry data
    write(unit,'("# Geometry Data (.GEO file)")')
    write(unit,'("print(f''Processing {",I0,"} geometry blocks'')")') sil%blocks%n
    do i = 1, sil%blocks%n
      associate(block => sil%blocks%d(i))
        matlab_varname = trim(block%name)
        do j = 1, len_trim(matlab_varname)
          if (matlab_varname(j:j) == '-') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '#') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '.') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == ' ') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == '(') matlab_varname(j:j) = '_'
          if (matlab_varname(j:j) == ')') matlab_varname(j:j) = '_'
        end do

        ! Ensure name doesn't start with underscore or digit
        if (len_trim(matlab_varname) > 0) then
          if (matlab_varname(1:1) == '_' .or. &
                (matlab_varname(1:1) >= '0' .and. matlab_varname(1:1) <= '9')) then
            matlab_varname = 'var_' // trim(matlab_varname)
          end if
        end if

        ! Truncate to 31 characters (leaving room for _info suffix)
        if (len_trim(matlab_varname) > 26) then
          matlab_varname = matlab_varname(1:26)
        end if

        write(unit,'("print(''Block ",I0,": ",A," dtype=",I0," ndata=",I0,"'')")') &
              i, trim(matlab_varname), block%dtype, block%ndata

        select case(block%dtype)
        case(GALDATA_INT)
          write(unit,'("if ",I0," > 0:")') block%ndata
          write(unit,'("    data[''geo''][''",A,"''] = np.array([")') trim(matlab_varname)
          do j = 1, block%ndata
            if (j < block%ndata) then
              write(unit,'(I0,",")') block%idata(j)
            else
              write(unit,'(I0)') block%idata(j)
            end if
          end do
          write(unit,'("    ], dtype=int)")')
          write(unit,'("else:")')
          write(unit,'("    data[''geo''][''",A,"''] = np.array([], dtype=int)")') trim(matlab_varname)
        case(GALDATA_REAL)
          write(unit,'("if ",I0," > 0:")') block%ndata
          write(unit,'("    data[''geo''][''",A,"''] = np.array([")') trim(matlab_varname)
          do j = 1, block%ndata
            if (j < block%ndata) then
              write(unit,'(ES15.8,",")') denorm(block%rdata(j), trim(block%unit))
            else
              write(unit,'(ES15.8)') denorm(block%rdata(j), trim(block%unit))
            end if
          end do
          write(unit,'("    ], dtype=float)")')
          write(unit,'("else:")')
          write(unit,'("    data[''geo''][''",A,"''] = np.array([], dtype=float)")') trim(matlab_varname)
        case(GALDATA_CHAR)
          write(unit,'("if len(''",A,"'') > 0:")') trim(block%cdata)
          write(unit,'("    data[''geo''][''",A,"''] = ''",A,"''")') trim(matlab_varname), trim(block%cdata)
          write(unit,'("else:")')
          write(unit,'("    data[''geo''][''",A,"''] = ''''")') trim(matlab_varname)
        case default
          write(unit,'("print(f''Unknown dtype ",I0," for block ",A,"'')")') block%dtype, trim(matlab_varname)
          write(unit,'("data[''geo''][''",A,"''] = None  # Unknown dtype")') trim(matlab_varname)
        end select

        write(unit,'("data[''geo''][''",A,"_info''] = {")') trim(matlab_varname)
        write(unit,'("    ''name'': ''",A,"'',")') trim(block%name)
        write(unit,'("    ''unit'': ''",A,"'',")') trim(block%unit)
        write(unit,'("    ''dtype'': ",I0,",")') block%dtype
        write(unit,'("    ''ndata'': ",I0)') block%ndata
        write(unit,'("}")')
        write(unit,'("")')
      end associate
    end do

    write(unit,'("print(''Saving to galene_database.mat...'')")')
    write(unit,'("savemat(''galene_database.mat'', data)")')
    write(unit,'("print(''Successfully saved galene_database.mat'')")')
    close(unit)

    ! Execute Python to create .mat file
    write(*,*) "Executing Python script to create .mat file..."
    call system('python3 ' // trim(temp_pyfile))
    write(*,*) "Python execution completed."
    call system('rm -f ' // trim(temp_pyfile))

    write(*,*) "MATLAB .mat database exported to: galene_database.mat"
    write(*,*) "Access data in MATLAB with: load('galene_database.mat'); data.sim.variable_name"
  end subroutine

  subroutine export_to_tikz_standalone(gat, sil)
    !! Export data to standalone TikZ figure
    type(gal_file), intent(in) :: gat, sil
    character(256) :: var1_name, var2_name, outfile
    type(gal_block), pointer :: var1_block, var2_block
    integer :: unit, i
    real :: xmin, xmax, ymin, ymax, xrange, yrange

    write(*,*) "Enter first variable name (X-axis): "
    read(*,*) var1_name
    write(*,*) "Enter second variable name (Y-axis): "
    read(*,*) var2_name

    ! Try to find variables
    var1_block => gat%get_block(trim(var1_name))
    if (.not. associated(var1_block)) then
      var1_block => sil%get_block(trim(var1_name))
    end if

    var2_block => gat%get_block(trim(var2_name))
    if (.not. associated(var2_block)) then
      var2_block => sil%get_block(trim(var2_name))
    end if

    if (.not. associated(var1_block)) then
      write(*,*) "Variable '", trim(var1_name), "' not found!"
      return
    end if

    if (.not. associated(var2_block)) then
      write(*,*) "Variable '", trim(var2_name), "' not found!"
      return
    end if

    if (var1_block%dtype /= GALDATA_REAL .or. var2_block%dtype /= GALDATA_REAL) then
      write(*,*) "Both variables must be real-valued for plotting!"
      return
    end if

    if (var1_block%ndata /= var2_block%ndata) then
      write(*,*) "Variables must have the same number of data points!"
      return
    end if

    outfile = trim(var1_name) // '_vs_' // trim(var2_name) // '.tex'
    open(newunit=unit, file=outfile, status='replace')

    ! Write standalone LaTeX document
    write(unit,'("\\documentclass[border=5pt]{standalone}")')
    write(unit,'("\\usepackage{tikz}")')
    write(unit,'("\\usepackage{pgfplots}")')
    write(unit,'("\\pgfplotsset{compat=1.18}")')
    write(unit,'("\\begin{document}")')
    write(unit,'("\\begin{tikzpicture}")')

    ! Calculate axis ranges
    xmin = minval(var1_block%rdata)
    xmax = maxval(var1_block%rdata)
    ymin = minval(var2_block%rdata)
    ymax = maxval(var2_block%rdata)
    xrange = xmax - xmin
    yrange = ymax - ymin

    write(unit,'("\\begin{axis}[")')
    write(unit,'("  xlabel={\\textbf{",A,"} [",A,"]},")') trim(var1_block%name), trim(var1_block%unit)
    write(unit,'("  ylabel={\\textbf{",A,"} [",A,"]},")') trim(var2_block%name), trim(var2_block%unit)
    write(unit,'("  grid=major,")')
    write(unit,'("  legend pos=north west,")')
    write(unit,'("  width=10cm,")')
    write(unit,'("  height=8cm,")')
    write(unit,'("  xmin=",F12.6,",")') xmin - 0.05*xrange
    write(unit,'("  xmax=",F12.6,",")') xmax + 0.05*xrange
    write(unit,'("  ymin=",F12.6,",")') ymin - 0.05*yrange
    write(unit,'("  ymax=",F12.6,",")') ymax + 0.05*yrange
    write(unit,'("]")')

    write(unit,'("\\addplot[blue, mark=*, mark size=1pt] coordinates {")')
    do i = 1, var1_block%ndata
      write(unit,'("(",F12.6,",",F12.6,")")') var1_block%rdata(i), var2_block%rdata(i)
    end do
    write(unit,'("};")')

    write(unit,'("\\legend{",A," vs ",A,"}")') trim(var1_block%name), trim(var2_block%name)
    write(unit,'("\\end{axis}")')
    write(unit,'("\\end{tikzpicture}")')
    write(unit,'("\\end{document}")')

    close(unit)
    write(*,*) "Standalone TikZ figure exported to: ", trim(outfile)
    write(*,*) "Compile with: pdflatex ", trim(outfile)
  end subroutine

end program galene_converter
