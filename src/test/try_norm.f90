program try_norm
  use normalization_m

  implicit none

  character(len=1024) :: line
  character(:), allocatable :: input
  integer :: size
  real :: T, val

  write(*, "(A)", advance="no") "Temperature? "
  read(*,*) T
  call init_normconst(T)

  write(*, *) "Type .exit to exit the interpreter"
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
      val = eval(input)
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
    do while (tok%type /= token_eof)
      write(*, "(A)", advance="no") " " // token_name(tok%type) // " '" // tok%literal // "'"
      call l%next_token(tok)
    end do
    write (*, "(A)") ""

  end subroutine
    
  subroutine run_line(input)
    character(:), allocatable, intent(in) :: input
    type(lexer) :: l
    type(parser) :: p
    
    class(ast_t), allocatable :: ast
    integer :: i

    call l%init(input)
    call p%init(l)
    call p%parse_expr(ast)

    if (allocated(p%errors)) then
      write(*, *) "no Apfeltasche for you!"
      do i = 1, count(p%errors /= "")
        write(*, *) trim(p%errors(i))
      end do
    else
      write(*, *) ast_string(ast)
    end if
  end subroutine run_line

  recursive function ast_string(this) result(s)
    class(ast_t), intent(in)  :: this
    character(:), allocatable :: s

    select type(this)
      type is (ast_expr)
        s = ast_expr_string(this)
      type is (ast_unit)
        s = ast_unit_string(this)
      type is (ast_int)
        s = ast_int_string(this)
      type is (ast_expo)
        s = ast_expo_string(this)
      type is (ast_ratio)
        s = ast_ratio_string(this)
      class default
        m4_ifdef({m4_intel},{call tracebackqq()})
        m4_ifdef({m4_gnu},{call backtrace()})
        call exit(1)
      end select
    end function

  function ast_expr_string(expr) result(s)
    class(ast_expr), intent(in) :: expr

    character(:), allocatable :: s
    integer                   :: i
    
    if (expr%operands%n > 0) then
      s = "( " // ast_string(expr%operands%d(1)%item)
    else
      s = "( "
    end if

    do i = 2, expr%operands%n 
      select case (expr%operators%d(i-1))
        case ("*")
          s = s // " * " // ast_string(expr%operands%d(i)%item)
        case ("/")
          s = s // " / " // ast_string(expr%operands%d(i)%item)
      end select
    end do
    s = s // " )"
  end function
    
  function ast_unit_string(this) result(s)
    type(ast_unit), intent(in) :: this
    character(:), allocatable :: s
    s = "unit[" // this%prefix%s // "-" // this%unit%s // "]"
  end function
    
  function ast_int_string(this) result(s)
    type(ast_int), intent(in) :: this
    character(:), allocatable :: s
    s = this%tok%literal
  end function

  function ast_expo_string(this) result(s)
    type(ast_expo), intent(in) :: this
    character(:), allocatable :: s
    s = ast_string(this%base) // "^" // ast_string(this%expo)
  end function
  
  function ast_ratio_string(this) result(s)
    type(ast_ratio), intent(in) :: this
    character(:), allocatable :: s
    s = "(" // ast_string(this%nom) // "/" // ast_string(this%den) // ")"
  end function

end program
