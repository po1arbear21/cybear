module test_variable_m

  use grid_m,          only: IDX_VERTEX
  use grid1D_m,        only: grid1D
  use normalization_m, only: denorm, norm, init_normconst, destruct_normconst
  use test_case_m,     only: test_case
  use variable_m,      only: variable

  implicit none

  private
  public test_variable

  type, extends(variable) :: potential
    !! potential variable. just for testing purposes here.
  contains
    procedure :: init => potential_init
  end type

contains

  subroutine potential_init(this, name, g, idx_type)
    class(potential), intent(out) :: this
    character(*),     intent(in)  :: name
    type(grid1D),     intent(in)  :: g
      !! position grid
    integer,          intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX or IDX_CELL)

    ! init base
    call this%variable_init(name, 'V', g, idx_type, 0)
  end subroutine

  subroutine test_variable()
    type(test_case) :: tc
    type(grid1D)    :: g
    type(potential) :: pot

    call tc%init("variable")
    call init_normconst(300.0)

    ! init grid, variable
    call g%init([1.0, 2.0, 4.0])
    call pot%init("potential", g, IDX_VERTEX)

    ! test1: load data
    ! test2: save data
    block
      character(*), parameter :: fname_tmp = "/tmp/variable_test.csv"
      real,         parameter :: d_exp(*)  = [10.0, -2.1, 3.5]
      integer :: iounit
      real    :: d(3)

      ! test1: load data
      call pot%load_data("src/test/variable_test.csv")
      call tc%assert_eq(d_exp, denorm(pot%get(), 'V'), 0.0, "load data")

      ! test2: save data
      call pot%save_data(fname_tmp)

      open (newunit = iounit, file = fname_tmp, action = 'READ')
      read (iounit, *) d
      close (iounit)
      call execute_command_line("rm " // fname_tmp)
      call tc%assert_eq(d_exp, d, 0.0, "save data")
    end block

    call destruct_normconst()
    call tc%finish()
  end subroutine

end module
