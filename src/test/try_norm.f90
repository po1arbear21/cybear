program try_norm

  use normalization_m, only: init_normconst, lexer, norm, parser, token, token_name, tree

  implicit none

  character(len=1024) :: line
  character(:), allocatable :: input
  integer :: size
  real :: T, val

  write(*, "(A)", advance="no") "Temperature? "
  read(*,*) T
  call init_normconst(T)

  print "(A)", "Type .exit or .quit to exit the interpreter"
  do
    write (*, "(A)", advance = "no") ">> "
    read (*, "(A)") line

    line = trim(line)
    if (line == ".exit" .or. line == ".quit") then
      exit
    else if (line /= "") then
      size = len_trim(line)
      input = line(1:size)
      print "(A)", "Lexer:"
      call run_line_lexer(input)
      print "(A)", "Parser:"
      call run_line(input)
      print "(A)", "Normalizing constant:"
      val = norm(1.0, input)
      print "(ES25.16)", val
    end if
  end do

contains

  subroutine run_line_lexer(input)
    character(:), allocatable, intent(in) :: input
    type(lexer) :: l
    type(token) :: tok

    call l%init(input)

    call l%next_token(tok)
    do while (tok%type /= 0)
      write(*, "(A)", advance = "no") " " // token_name(tok%type) // " '" // tok%literal // "'"
      call l%next_token(tok)
    end do
    print *

  end subroutine

  subroutine run_line(input)
    character(:), allocatable, intent(in) :: input
    type(lexer) :: l
    type(parser) :: p

    class(tree), allocatable :: tree
    integer :: i

    call l%init(input)
    call p%init(l)
    call p%parse_expr(tree)

    if (p%errors%n > 0) then
      print "(A)", "no Apfeltasche for you!"
      do i = 1, p%errors%n
        print "(A)", p%errors%d(i)%s
      end do
    else
      print "(A)", tree%to_string()
    end if
  end subroutine run_line

end program
