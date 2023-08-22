m4_include(macro.f90.inc)

module normalization_m

  use error_m,  only: program_error
  use map_m,    only: map_string_real, mapnode_string_real
  use math_m,   only: PI
  use string_m, only: string, new_string
  use vector_m, only: 

  implicit none

  ! Uff
  public

  type :: normalization
    type(map_string_real) :: map
    type(map_string_real) :: prefixes
    type(map_string_real) :: units
  contains
    procedure :: init     => normalization_init
    procedure :: destruct => normalization_destruct
  end type

  type(normalization) :: normconst

  ! maximum supported array dimension
  m4_define({m4_max_dim},{8})

  ! get actual type
  m4_define({m4_norm_type},{m4_ifelse($1,r,real,{m4_ifelse($1,c,complex,)})})

  ! get operation
  m4_define({m4_norm_op},{m4_ifelse($1,norm,/,{m4_ifelse($1,denorm,*,)})})

  ! combine dimensions 1 to max_dim with norm/denorm and r/c to get a full list
  m4_define({m4_list_help},{
    m4_ifelse($2,0,,{m4_list_help($1,m4_decr($2),$3)})
    m4_X($1,$2,$3)
  })
  m4_define({m4_list},{
    m4_list_help(norm,m4_max_dim,r)
    m4_list_help(denorm,m4_max_dim,r)
    m4_list_help(norm,m4_max_dim,c)
    m4_list_help(denorm,m4_max_dim,c)
  })

  m4_define({m4_X},{
    interface $1
      module procedure :: $1_$2_$3
    end interface
  })
  m4_list

  ! Tokens
  enum, bind(c)
    enumerator :: token_eof = 0, token_illegal
    enumerator :: token_ident
    enumerator :: token_int

    enumerator :: token_asterisk, token_slash, token_caret
    enumerator :: token_lparen, token_rparen
  end enum

  type token
    integer :: type
    character(:), allocatable :: literal
  end type

  m4_define({m4_list_vec},{
    m4_X(ast_item)
    m4_X(character)
  })

  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(vector_def.f90.inc)
  })
  m4_list_vec

  ! Abstract Syntax Tree elements 
  type, abstract :: ast_t
  contains
    procedure :: eval => ast_eval
  end type
  
  type, extends(ast_t) :: ast_expr
    type(vector_ast_item)  :: operands
    type(vector_character) :: operators
  end type
  type, extends(ast_t) :: ast_unit
    type(token)  :: tok
    type(string) :: prefix
    type(string) :: unit
  end type
  type, extends(ast_t) :: ast_int
    type(token) :: tok
    integer     :: value
  end type
  type, extends(ast_t) :: ast_ratio
    class(ast_t), allocatable :: nom
    class(ast_t), allocatable :: den
  end type
  type, extends(ast_t) :: ast_expo
    class(ast_t), allocatable :: base
    type(ast_ratio)           :: expo
  end type
  
  type :: lexer
    character(:), allocatable :: input
    integer :: pos
    integer :: read_pos
    character :: ch
  contains
    procedure :: init        => lexer_init
    procedure :: next_token  => lexer_next_token
    procedure :: read_ch     => lexer_read_ch
    procedure :: read_ident  => lexer_read_ident
    procedure :: read_number => lexer_read_number
    procedure :: peek_char   => lexer_peek_char
    procedure :: skip_ws     => lexer_skip_ws
  end type

  type :: parser
    type(lexer) :: l
    character(:), allocatable :: errors(:)
    type(token) :: curr_token
    type(token) :: peek_token
  contains 
    procedure :: init             => parser_init
    procedure :: add_error        => parser_add_error
    procedure :: expect_peek      => parser_expect_peek
    procedure :: next_token       => parser_next_token
    procedure :: parse_expr       => parser_parse_expr
    procedure :: parse_value      => parser_parse_value
    procedure :: parse_expo       => parser_parse_expo
    procedure :: parse_group      => parser_parse_group
    procedure :: parse_identifier => parser_parse_identifier
    procedure :: parse_integer    => parser_parse_integer
    procedure :: peek_error       => parser_peek_error
    procedure :: peek_token_is    => parser_peek_token_is
  end type

  type :: ast_item
    class(ast_t), allocatable :: item
  end type

contains
    
  subroutine init_normconst(T)
    !! initialize global normalization object
    real, intent(in) :: T

    call normconst%init(T)
  end subroutine

  subroutine destruct_normconst()
    !! destruct global normalization object

    call normconst%destruct()
  end subroutine

  subroutine normalization_init(this, T)
    class(normalization), intent(out) :: this
    real,                 intent(in)  :: T
    
    ! constants
    real, parameter :: EC = 1.602176634e-19
      !! elementary charge [ As ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?e
    real, parameter :: EM = 9.1093837015e-31
      !! electron rest mass [ kg ]
      !! rel err: 3e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?me
    real, parameter :: PLANCK = 6.62607015e-34 / (2*PI*EC)
      !! reduced Planck constant [ eVs ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?h
    real, parameter :: BOLTZ = 1.380649e-23/EC
      !! Boltzmann constant [ eV/K ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?k
    real, parameter :: EPS0 = 8.8541878128e-12
      !! vacuum permittivity [ F/m ]
      !! rel err: 1.5e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?ep0
    real, parameter :: MU0 = 1.25663706212e-6
      !! vacuum permeability [ N/A² ]
      !! rel err: 1.5e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?mu0

    ! local variables
    real :: ampere, coulomb, diel, farad, henry, hertz, kelvin, gram, meter, ohm, permby, second, volt, watt

    call this%map%init()
    call this%prefixes%init()
    call this%units%init()

    ! metric prefixes: INVERSE values for normalization
    call this%prefixes%insert(new_string("P"), 1e-15)
    call this%prefixes%insert(new_string("T"), 1e-12)
    call this%prefixes%insert(new_string("G"), 1e-9)
    call this%prefixes%insert(new_string("M"), 1e-6)
    call this%prefixes%insert(new_string("k"), 1e-3)
    call this%prefixes%insert(new_string("h"), 1e-2)
    call this%prefixes%insert(new_string("1"), 1e0)
    call this%prefixes%insert(new_string("d"), 1e1)
    call this%prefixes%insert(new_string("c"), 1e2)
    call this%prefixes%insert(new_string("m"), 1e3)
    call this%prefixes%insert(new_string("u"), 1e6)
    call this%prefixes%insert(new_string("n"), 1e9)
    call this%prefixes%insert(new_string("p"), 1e12)
    call this%prefixes%insert(new_string("f"), 1e15)
    call this%prefixes%insert(new_string("a"), 1e18)
    
    ! available units
    volt    = BOLTZ * T
    meter   = PLANCK / sqrt(EM / EC * volt)
    second  = PLANCK / volt
    hertz   = 1.0 / second
    gram    = EM * 1e3 ! Use grams instead of kilograms
    ampere  = EC / second
    coulomb = ampere * second
    farad   = coulomb / volt
    ohm     = volt / ampere
    henry   = ohm * second
    watt    = volt * ampere
    kelvin  = T
    diel    = sqrt(EM * EC / volt) / (PLANCK * EPS0)
    permby  = henry / meter / MU0

    call this%units%insert(new_string("V"), volt)
    call this%units%insert(new_string("eV"), volt)
    call this%units%insert(new_string("m"), meter)
    call this%units%insert(new_string("s"), second)
    call this%units%insert(new_string("Hz"), hertz)
    call this%units%insert(new_string("g"), gram)
    call this%units%insert(new_string("A"), ampere)
    call this%units%insert(new_string("C"), coulomb)
    call this%units%insert(new_string("F"), farad)
    call this%units%insert(new_string("Ohm"), ohm)
    call this%units%insert(new_string("H"), henry)
    call this%units%insert(new_string("W"), watt)
    call this%units%insert(new_string("K"), kelvin)
    call this%units%insert(new_string("Å"), 1e10 * meter)
    call this%units%insert(new_string("eps0"), diel)
    call this%units%insert(new_string("mu0"), permby)
    call this%units%insert(new_string("deg"), 180 / PI)
    call this%units%insert(new_string("rad"), 1.0)
  end subroutine

  function get_norm_value(unit, n) result (val)
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n

    real :: val
    type(string) :: unit_
    type(mapnode_string_real), pointer :: node

    if (present(n)) call program_error("Calling normalization with custom normalization object is deprecated.")

    unit_ = new_string(unit)
    node => normconst%map%find(unit_)
    if (associated(node)) then
      val = node%value
      return
    end if

    val = eval(unit_%s)
    call normconst%map%set(unit_, val)
  end function

  function eval(unit) result(res)
    character(:), allocatable, intent(in) :: unit
    real :: res
   
    class(ast_t), allocatable :: ast
    integer :: i

    type(lexer) :: l
    type(parser) :: p

    call l%init(unit)
    call p%init(l)
    call p%parse_expr(ast)

    if (allocated(p%errors)) then
      if (count(p%errors /= "") > 0) then
        do i = 1, count(p%errors /= "")
          write(*, *) trim(p%errors(i))
        end do
        call program_error("Cannot parse unit expression '" // unit // "'")
      end if
    end if

    res = ast%eval()
  end function

  m4_define({m4_nvalues_shape},{m4_ifelse($1,1,size(values,1),{m4_nvalues_shape(m4_decr($1)),size(values,$1)})})
    m4_define({m4_nvalues_pshape},{m4_ifelse($1,0,,{(m4_nvalues_shape($1))})})
    m4_define({m4_X},{function $1_$2_$3(values, unit, n) result(nvalues)
      m4_norm_type($3),              intent(in) :: values{}m4_pshape($2)
        !! value to (de-)normalize
      character(*),                  intent(in) :: unit
        !! physical unit token
      type(normalization), optional, intent(in) :: n
        !! optional normalization object (default: use global normconst)
      m4_norm_type($3)                          :: nvalues{}m4_nvalues_pshape($2)
        !! return normalized value
  
      nvalues = values m4_norm_op($1) get_norm_value(unit, n = n)
    end function})
    m4_list

  subroutine normalization_destruct(this)
    !! destruct normalization constants (release memory)
    class(normalization), intent(inout) :: this

    call this%map%destruct()
    call this%prefixes%destruct()
    call this%units%destruct()
  end subroutine

  pure function token_name(token_type) result(name)
    integer, intent(in) :: token_type
    character(:), allocatable :: name

    select case (token_type)
      case (token_eof)
        name = "EOF"
      case (token_illegal)
        name = "ILLEGAL"
      case (token_ident)
        name = "IDENT"
      case (token_int)
        name = "INT"
      case (token_asterisk)
        name = "ASTERISK"
      case (token_slash)
        name = "SLASH"
      case (token_caret)
        name = "CARET"
      case (token_lparen)
        name = "LPAREN"
      case (token_rparen)
        name = "RPAREN"
      case default
        name = "OTHER"
    end select
  end function token_name
    
  subroutine lexer_init(this, input)
    class(lexer),              intent(inout) :: this
    character(:), allocatable, intent(in)    :: input

    this%input = input
    this%pos = 1
    this%read_pos = 1
    this%ch = "x"
    call this%read_ch()
  end subroutine
    
  subroutine lexer_next_token(this, tok)
    class(lexer), intent(inout) :: this
    type(token),  intent(out)   :: tok

    character(:), allocatable :: lit

    call this%skip_ws()

    select case(this%ch)
      case ("/")
        tok = token(token_slash, this%ch)
      case ("*")
        tok = token(token_asterisk, this%ch)
      case ("^")
        tok = token(token_caret, this%ch)
      case ("(")
        tok = token(token_lparen, this%ch)
      case (")")
        tok = token(token_rparen, this%ch)
      case (achar(0))
        tok = token(token_eof, "")
      case ("-")
        if (is_digit(this%peek_char())) then
          call this%read_number(lit)
          tok = token(token_int, lit)
          return
        else
          tok = token(token_illegal, this%ch)
        end if
      case default
        if (is_letter(this%ch)) then
          call this%read_ident(lit)
          tok = token(token_ident, lit)
          return
        else if (is_digit(this%ch)) then
          call this%read_number(lit)
          tok = token(token_int, lit)
          return
        else
          tok = token(token_illegal, this%ch)
        end if
      end select

      call this%read_ch()
  end subroutine
    
  subroutine lexer_read_ch(this)
    class(lexer), intent(inout) :: this

    if (this%read_pos > len(this%input)) then
      this%ch = achar(0)
    else
      this%ch = this%input(this%read_pos:this%read_pos)
    end if
    this%pos = this%read_pos
    this%read_pos = this%read_pos + 1
  end subroutine
    
  subroutine lexer_skip_ws(this)
    class(lexer), intent(inout) :: this

    do while (is_whitepsace(this%ch))
      call this%read_ch()
    end do
  end subroutine

  subroutine lexer_read_ident(this, ident)
    class(lexer),              intent(inout) :: this
    character(:), allocatable, intent(out)   :: ident
    integer :: position

    position = this%pos
    do while (is_letter(this%ch) .or. is_digit(this%ch))
      call this%read_ch()
    end do

    ident = this%input(position:this%pos-1)
  end subroutine

  subroutine lexer_read_number(this, number)
    class(lexer),              intent(inout) :: this
    character(:), allocatable, intent(out)   :: number
    integer :: position

    position = this%pos
    if (this%ch == "-") call this%read_ch()

    do while (is_digit(this%ch))
      call this%read_ch()
    end do

    number = this%input(position:this%pos-1)
  end subroutine

  function lexer_peek_char(this) result(ch)
    class(lexer), intent(in) :: this

    character :: ch

    if (this%read_pos > len(this%input)) then
        ch = achar(0)
    else
        ch = this%input(this%read_pos:this%read_pos)
    end if
  end function
    
  pure function is_whitepsace(ch) result(ret)
    character, intent(in) :: ch
    logical :: ret
    integer :: ord

    ord = ichar(ch)
    if (ord == 32 .or. ord == 9 .or. ord == 13 .or. ord == 10) then
      ret = .true.
    else
      ret = .false.
    end if
  end function

  pure function is_letter(ch) result(ret)
    character, intent(in) :: ch
    logical :: ret

    if ((ch >= 'a' .and. ch <= 'z') .or. (ch >= 'A' .and. ch <= 'Z') .or. ch == '_') then
      ret = .true.
    else
      ret = .false.
    end if
  end function

  pure function is_digit(ch) result(ret)
    character, intent(in) :: ch
    logical :: ret

    if (ch >= '0' .and. ch <= '9') then
      ret = .true.
    else
      ret = .false.
    end if
  end function
    
  subroutine parser_init(this, l)
    class(parser), intent(inout) :: this
    type(lexer),   intent(inout) :: l

    this%l = l
    call this%next_token()
    call this%next_token()
  end subroutine

  subroutine parser_next_token(this)
    class(parser), intent(inout) :: this

    this%curr_token = this%peek_token
    call this%l%next_token(this%peek_token)
  end subroutine
    
  function parser_peek_token_is(this, tok) result(b)
    class(parser), intent(in) :: this
    integer,       intent(in) :: tok

    logical :: b

    b = tok == this%peek_token%type
  end function
    
  function parser_expect_peek(this, tok) result(b)
    class(parser), intent(inout) :: this
    integer,         intent(in)    :: tok
    
    logical :: b

    if (this%peek_token_is(tok)) then
      call this%next_token()
      b = .true.
    else
      call this%peek_error(tok)
      b = .false.
    end if
  end function

  subroutine parser_peek_error(this, tok)
    class(parser), intent(inout) :: this
    integer,       intent(in)    :: tok

    character(:), allocatable :: msg

    msg = "expected next token to be " // token_name(tok) // ", got " // token_name(this%peek_token%type)
    call this%add_error(msg)
  end subroutine

  subroutine parser_add_error(this, msg)
    class(parser), intent(inout) :: this
    
    character(:), allocatable, intent(in) :: msg
    
    if (allocated(this%errors)) then
      this%errors = [this%errors, msg]
    else
      this%errors = [msg]
    end if
  end subroutine  

  subroutine parser_parse_expr(this, ast, nested)
    class(parser),             intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: ast
    logical, optional,         intent(in)    :: nested

    type(ast_expr)            :: expr
    type(ast_item)            :: item
    character(:), allocatable :: err
    character                 :: op
    logical                   :: nested_

    nested_ = .false.
    if (present(nested)) nested_ = nested

    ast = expr
    select type(ast)
    type is (ast_expr)

    call ast%operands%init(0, c = 4)
    call ast%operators%init(0, c = 4)
    
    do while(this%curr_token%type /= token_eof .and. ((.not. nested_) .or. (this%curr_token%type /= token_rparen)))
      
      if (this%curr_token%type == token_asterisk .or. this%curr_token%type == token_slash) then
        op = this%curr_token%literal
        call this%next_token()
      else
        op = achar(0)
      end if

      if (ast%operands%n > 0 .and. op == achar(0)) then
        err = "Got no operator between values"
        call this%add_error(err)
        exit
      end if

      call this%parse_value(item%item)
      if (allocated(item%item)) then
        call ast%operands%push(item)
        if (op /= achar(0)) call ast%operators%push(op)
        deallocate(item%item)
      end if

      call this%next_token()
    end do
    end select
  end subroutine

  subroutine parser_parse_value(this, s)
    class(parser),           intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: s
    
    character(:), allocatable :: err

    select case(this%curr_token%type)
      case (token_ident)
        call this%parse_identifier(s)
      case (token_int)
        call this%parse_integer(s)
      case (token_lparen)
        call this%parse_group(s)
      case default
        err = "expected expression got " // token_name(this%curr_token%type)
        call this%add_error(err)
        return
    end select

    ! Possible exponentiation
    if (this%peek_token_is(token_caret)) then
      call this%next_token()
      call this%parse_expo(s)
    end if
  end subroutine

  subroutine parser_parse_group(this, ast)
    class(parser),             intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: ast

    character(:), allocatable :: err

    call this%next_token()
    call this%parse_expr(ast, nested=.true.)

    if (.not. this%curr_token%type == token_rparen) then
      deallocate(ast)
      err = "expected ')' after expression got " // token_name(this%curr_token%type)
      call this%add_error(err)
      return
    end if
  end subroutine
  
  subroutine parser_parse_expo(this, ast)
    class(parser),             intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: ast

    type(ast_expo)  :: expo
    type(ast_ratio) :: ratio
    type(ast_int)   :: i

    character(:), allocatable  :: err
    logical                    :: b

    if (.not. allocated(ast)) then
      call program_error("Cannot parse exponent if no base is given")
    end if

    expo%base = ast
    deallocate(ast)
    
    call this%next_token()
    select case(this%curr_token%type)
      case (token_int) 
        call this%parse_integer(ratio%nom)
        
        i%tok = token(token_int, "1")
        i%value = 1
        ratio%den = i
      case (token_lparen)
        b = this%expect_peek(token_int)
        if (.not. b) return
        call this%parse_integer(ratio%nom)
        b = this%expect_peek(token_slash)
        if (.not. b) return
        b = this%expect_peek(token_int)
        if (.not. b) return
        call this%parse_integer(ratio%den)
        b = this%expect_peek(token_rparen)
        if (.not. b) return
      case default
        err = "Expected integer or left parenthesis, got " // token_name(this%curr_token%type)
        call this%add_error(err)
        return
    end select

    expo%expo = ratio
    ast = expo
  end subroutine
    
  subroutine parser_parse_integer(this, ast)
    class(parser),             intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: ast

    type(ast_int) :: s

    s%tok = this%curr_token
    read(this%curr_token%literal, *) s%value ! todo: handle read error
    ast = s
  end subroutine

  subroutine parser_parse_identifier(this, ast)
    class(parser),             intent(inout) :: this
    class(ast_t), allocatable, intent(inout) :: ast

    type(ast_unit)                     :: s
    type(mapnode_string_real), pointer :: node
    type(string)                       :: prefix, unit
    character(:), allocatable          :: err

    s%tok = this%curr_token
    unit = new_string(this%curr_token%literal)
    
    ! identifier = unit
    node => normconst%units%find(unit)
    if (associated(node)) then
      s%prefix = new_string("1")
      s%unit = unit
      ast = s
      return
    end if
    
    if (len(unit%s) .eq. 1) then
      err = "Cannot find unit " // unit%s
      call this%add_error(err)
      return
    end if

    ! identifier = prefix + unit?
    prefix = new_string(unit%s(1:1))
    unit = new_string(unit%s(2:))

    node => normconst%prefixes%find(prefix)
    if (.not. associated(node)) then
      err = "Undefined unit prefix " // prefix%s
      call this%add_error(err)
      return
    end if
    node => normconst%units%find(unit)
    if (.not. associated(node)) then
      err = "Undefined unit " // unit%s
      call this%add_error(err)
      return
    end if
    
    s%prefix = prefix
    s%unit = unit
    ast = s
  end subroutine

  recursive function ast_eval(this) result(val)
    !! Serves as eval and destruct
    class(ast_t), intent(inout) :: this
    
    real    :: val
    integer :: i

    select type(this)
      type is (ast_expr)
        if (this%operands%n < 1) call program_error("Empty expression given")
      
        val = this%operands%d(1)%item%eval()
        do i = 2, this%operands%n
          select case (this%operators%d(i-1))
            case ("*")
              val = val * this%operands%d(i)%item%eval()
            case ("/")
              val = val / this%operands%d(i)%item%eval()
            case default
              call program_error("Undefined operator given")
          end select
        end do

        call this%operands%destruct()
        call this%operators%destruct()
      type is (ast_int)
        val = real(this%value)
      type is (ast_unit)
        val = normconst%prefixes%get(this%prefix) * normconst%units%get(this%unit)
      type is (ast_expo)
        val = this%base%eval() ** this%expo%eval()
      type is (ast_ratio)
        val = this%nom%eval() / this%den%eval()
      class default
        call program_error("Got undefined ast type")
    end select
  end function

  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(vector_imp.f90.inc)
  })
  m4_list_vec
end module
