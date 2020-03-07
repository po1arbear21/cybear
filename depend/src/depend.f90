program depend
  use src_file_m
  implicit none

  type(src_file), allocatable :: files(:)
  type(vector_src_item)       :: programs
  type(vector_src_item)       :: modules
  type(vector_src_item)       :: submodules
  type(vector_string)         :: filenames, include_folders
  integer                     :: i, j, k, l, expect
  logical                     :: file_exists
  type(src_item)              :: tmp
  character(len=SLEN)         :: arg, incfolders, folder, filename, buildfolder

  ! process command line arguments
  expect = 0
  call filenames%init(0, c = 16)
  incfolders = ""
  buildfolder = ""
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)

    if (expect > 0) then
      if (expect == 1) then
        buildfolder = trim(arg)
      elseif (expect == 2) then
        incfolders = trim(arg)
      end if
      expect = 0
    else
      if (len_trim(arg) >= 2) then
        if (arg(1:2) == "-b") then
          if (len_trim(arg) > 2) then
            buildfolder = arg(3:len_trim(arg))
          else
            expect = 1
          end if
        elseif (arg(1:2) == "-I") then
          if (len_trim(arg) > 2) then
            incfolders = arg(3:len_trim(arg))
          else
            expect = 2
          end if
        else
          filename = trim(arg)
          call filenames%push(filename)
        end if
      end if
    end if
  end do

  ! include folders, separated by colons
  call include_folders%init(0, c = 8)
  if (incfolders /= "") then
    do while (.true.)
      i = scan(trim(incfolders), ':')
      if (i == 0) then
        folder = trim(incfolders)
        call include_folders%push(folder)
        exit
      else
        if ((i == 1) .or. (i == len_trim(incfolders))) call error("include folder parse error")
        folder = trim(incfolders(1:i-1))
        call include_folders%push(folder)
        incfolders = trim(adjustl(incfolders(i+1:len_trim(incfolders))))
      end if
    end do
  end if

  ! init files
  allocate (files(filenames%n))
  do i = 1, size(files)
    arg = filenames%d(i)

    ! get slash
    j = scan(trim(arg), '/', back = .true.)
    if (j == 0) then
      folder = ""
      filename = arg
    else
      folder = trim(adjustl(arg(1:j)))
      filename = trim(adjustl(arg(j+1:len_trim(arg))))
    end if
    call files(i)%init(folder, filename)
  end do

  ! collect all programs, modules and submodules
  call programs%init(  0, c = 8)
  call modules%init(   0, c = 8)
  call submodules%init(0, c = 8)
  do i = 1, size(files)
    if (files(i)%program_name /= "") then
      do j = 1, programs%n
        if (programs%d(j)%name == files(i)%program_name) then
          call error("program name multiple times defined. name: "//files(i)%program_name)
        end if
      end do
      tmp%name = files(i)%program_name
      tmp%file_index = i
      call programs%push(tmp)
    end if

    do j = 1, files(i)%modules%n
      do k = 1, modules%n
        if (modules%d(k)%name == files(i)%modules%d(j)) then
          call error("module name defined multiple times. name: "//files(i)%modules%d(j))
        end if
      end do
      tmp%name = trim(files(i)%modules%d(j))
      tmp%file_index = i
      call modules%push(tmp)
    end do
    do j = 1, files(i)%submodules%n
      do k = 1, submodules%n
        if (submodules%d(k)%name == files(i)%submodules%d(j)) then
          call error("submodule name defined multiple times. name: "//files(i)%submodules%d(j))
        end if
      end do
      tmp%name = trim(files(i)%submodules%d(j))
      tmp%file_index = i
      call submodules%push(tmp)
    end do
  end do

  ! get dependency indices
  do i = 1, size(files)
    do j = 1, files(i)%use_modules%n
      do k = 1, modules%n
        if (modules%d(k)%name == files(i)%use_modules%d(j)) then
          files(i)%dep_modules(j) = k

          ! add anchor file of this module (if it does not exist already)
          do l = 1, files(i)%dep_anchor%n
            if (files(i)%dep_anchor%d(l) == modules%d(k)%file_index) exit
          end do
          if (l > files(i)%dep_anchor%n) call files(i)%dep_anchor%push(modules%d(k)%file_index)

          exit
        end if
      end do
    end do

    ! search include files in include folders
    do j = 1, files(i)%includes%n
      ! search file folder first
      inquire(file = trim(files(i)%folder)//trim(files(i)%includes%d(j)), exist = file_exists)
      if (file_exists) then
        files(i)%includes%d(j) = trim(files(i)%folder)//trim(files(i)%includes%d(j))
      else
        ! search include folders in order
        do k = 1, include_folders%n
          inquire(file = trim(include_folders%d(k))//trim(files(i)%includes%d(j)), exist = file_exists)
          if (file_exists) then
            files(i)%includes%d(j) = trim(include_folders%d(k))//trim(files(i)%includes%d(j))
            exit
          end if
        end do
      end if
      ! do not include file if not found in folders (e.g. includes of a library interface)
      if (.not. file_exists) then
        files(i)%includes%d(j) = ""
      end if
    end do
  end do

  ! output targets variable (programs)
  write(output_unit, "(1A)", advance = "no") "TARGETS ="
  do i = 1, programs%n
    write(output_unit, "(1A)", advance = "no") " "//trim(buildfolder)//trim(programs%d(i)%name)
  end do
  write(output_unit,*)
  write(output_unit,*)

  ! output object list for each target
  do i = 1, programs%n
    write(output_unit, "(1A)", advance = "no") "OBJECTS_"//trim(programs%d(i)%name)//" ="
    do j = 1, size(files)
      if ((files(j)%program_name == "") .or. (files(j)%program_name == programs%d(i)%name)) then
        write(output_unit, "(1A)", advance = "no") " "//trim(buildfolder)//trim(files(j)%objname)
      end if
    end do
    write(output_unit,*)
  end do
  write(output_unit,*)

  ! output anchor list for each target
  do i = 1, programs%n
    write(output_unit, "(1A)", advance = "no") "ANCHORS_"//trim(programs%d(i)%name)//" ="
    do j = 1, size(files)
      if ((files(j)%program_name == "") .or. (files(j)%program_name == programs%d(i)%name)) then
        write(output_unit, "(1A)", advance = "no") " "//trim(buildfolder)//trim(files(j)%ancname)
      end if
    end do
    write(output_unit,*)
  end do
  write(output_unit,*)

  ! output targets
  do i = 1, programs%n
    write(output_unit, "(1A)") ".PHONY: "//trim(programs%d(i)%name)
    write(output_unit, "(1A)") trim(buildfolder)//trim(programs%d(i)%name)//": $(ANCHORS_"//trim(programs%d(i)%name)//") $(OBJECTS_"//trim(programs%d(i)%name)//") $(OBJECTS_C) $(LIBS)"
    write(output_unit, "(1A)") char(9)//'@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) -I'//trim(buildfolder)//' -o $(OU_COL)'//trim(buildfolder)//trim(programs%d(i)%name)//'$(NO_COL) $(IN_COL)$(OBJECTS_'//trim(programs%d(i)%name)//') $(OBJECTS_C) $(LIBS)$(NO_COL)\n\n"'
    write(output_unit, "(1A)") char(9)//"@$(FC) $(FFLAGS) -I"//trim(buildfolder)//" -o "//trim(buildfolder)//trim(programs%d(i)%name)//" $(OBJECTS_"//trim(programs%d(i)%name)//") $(OBJECTS_C) $(LIBS)"
    write(output_unit,*)
  end do

  ! output anchor files
  do i = 1, size(files)
    write(output_unit, "(1A)", advance="no") trim(buildfolder)//trim(files(i)%ancname)//" : "//trim(files(i)%folder)//trim(files(i)%filename)

    ! depend on include files
    do k = 1, files(i)%includes%n
      if (files(i)%includes%d(k) /= "") then
        write(output_unit, "(1A)", advance="no") " "//trim(files(i)%includes%d(k))
      end if
    end do

    ! depend on anchor files
    do k = 1, files(i)%dep_anchor%n
      write(output_unit, "(1A)", advance = "no") " "//trim(buildfolder)//trim(files(files(i)%dep_anchor%d(k))%ancname)
    end do

    ! next line
    write(output_unit,*)
  end do

  ! output object files
  do i = 1, size(files)
    write(output_unit, "(1A)", advance="no") trim(buildfolder)//trim(files(i)%objname)//" : "//trim(files(i)%folder)//trim(files(i)%filename)

    ! depend on include files
    do k = 1, files(i)%includes%n
      if (files(i)%includes%d(k) /= "") then
        write(output_unit, "(1A)", advance="no") " "//trim(files(i)%includes%d(k))
      end if
    end do

    ! depend on anchor files
    do k = 1, files(i)%dep_anchor%n
      write(output_unit, "(1A)", advance = "no") " "//trim(buildfolder)//trim(files(files(i)%dep_anchor%d(k))%ancname)
    end do

    ! next line
    write (output_unit,*)
  end do
end program
