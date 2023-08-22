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

  write(*, *) "Type .exit or .quit to exit the interpreter"
  do
    write(*, '(A)', advance="no") ">> "
    read(*, '(A)') line

    line = trim(line)
    if (line == ".exit" .or. line == ".quit") then
      exit
    else if (line /= "") then
      size = len_trim(line)
      input = line(1:size)
      write(*,*) "Lexer:"
      call run_line_lexer(input)
      write(*,*) "Parser:"
      call run_line(input)
      write(*,*) "Normalizing constant:"
      val = norm(1.0, input)
      write(*,"(ES25.16)") val
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
      write(*, "(A)", advance="no") " " // token_name(tok%type) // " '" // tok%literal // "'"
      call l%next_token(tok)
    end do
    write (*, "(A)") ""

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

    if (allocated(p%errors)) then
      write(*, *) "no Apfeltasche for you!"
      do i = 1, count(p%errors /= "")
        write(*, *) trim(p%errors(i))
      end do
    else
      write(*, *) tree%to_string()
    end if
  end subroutine run_line

end program
