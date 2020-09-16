#include "macro.f90.inc"

module filesystem_m

  use error_m

  implicit none

  private
  public create_parent_dir
  public create_tmp_file
  public remove_file

contains

  function create_tmp_file() result(fname)
    !! Creates a temporary file, e.g. "/tmp/tmp.0zXjVs53bj"

    character(:), allocatable :: fname

    integer                   :: istat, iounit, ios
    character(:), allocatable :: tmp_fname

    ! create temp file and save its file name in tmp_fname
    allocate (character(0) :: tmp_fname)      ! remove gfortran warning
    tmp_fname = "/tmp/tmp142789.fname"

    call execute_command_line("mktemp > " // tmp_fname, cmdstat=istat)
    if (istat /= 0) call program_error('Could not create temporary file!')

    ! read temp file's fname
    open (newunit=iounit, file=tmp_fname, iostat=ios, action='read')
    if (ios /= 0) call program_error("Error opening file")

    block
      character(1024) :: fname_

      read (iounit, '(A)') fname_
      fname = trim(fname_)
    end block

    close (unit=iounit, iostat=ios)
    if (ios /= 0) call program_error("Error closing file")

    call remove_file(tmp_fname)
  end function

  subroutine create_parent_dir(fname)
    !! makes sure that the parent directory exists

    character(*), intent(in) :: fname
      !! file name, e.g. 'output/folder/my_data.asc'

    integer :: i

    i = scan(fname, '/', back=.true.)
    if (i /= 0) call execute_command_line("mkdir -p " // fname(:i))
  end subroutine

  subroutine remove_file(fname)
    !! removes a file

    character(*), intent(in) :: fname
      !! file name, e.g. 'output/folder/my_data.asc'

    integer :: istat

    call execute_command_line("rm " // fname, cmdstat=istat)
    if (istat /= 0) call program_error('Could not remove file!')
  end subroutine

end module
