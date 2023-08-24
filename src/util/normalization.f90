m4_include(macro.f90.inc)

module normalization_m

  use error_m,  only: program_error
  use map_m,    only: map_string_real, mapnode_string_real
  use math_m,   only: PI
  use string_m, only: string, new_string
  use util_m,   only: is_digit, is_letter, is_whitespace
  use vector_m, only: vector_char, vector_string

  implicit none

  private
  public normalization
  public init_normconst, destruct_normconst
  public norm, denorm
  public token, token_name, lexer, parser, tree

  type normalization
    type(map_string_real) :: prefixes
    type(map_string_real) :: units
    type(map_string_real) :: cache
  contains
    procedure :: init     => normalization_init
    procedure :: destruct => normalization_destruct
    procedure :: lookup   => normalization_lookup
  end type

  type(normalization) :: normconst

  ! maximum supported array dimension
  m4_define({m4_max_dim},{8})

  ! get actual type
  m4_define({m4_norm_type},{m4_ifelse($1,r,real,{m4_ifelse($1,c,complex,)})})

  ! get operation
  m4_define({m4_norm_op},{m4_ifelse($1,norm,*,{m4_ifelse($1,denorm,/,)})})

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
  integer, parameter :: TOKEN_EOF      = 0
  integer, parameter :: TOKEN_ILLEGAL  = 1
  integer, parameter :: TOKEN_IDENT    = 2
  integer, parameter :: TOKEN_INT      = 3
  integer, parameter :: TOKEN_ASTERISK = 4
  integer, parameter :: TOKEN_SLASH    = 5
  integer, parameter :: TOKEN_CARET    = 6
  integer, parameter :: TOKEN_LPAREN   = 7
  integer, parameter :: TOKEN_RPAREN   = 8

  type token
    integer                   :: type
    character(:), allocatable :: literal
  end type

  ! Abstract Syntax Tree nodes
  type, abstract :: tree
  contains
    procedure :: eval      => tree_eval
    procedure :: to_string => tree_to_string
  end type

  type item
    class(tree), allocatable :: item
  end type
  m4_define({T},{item})
  m4_include(vector_def.f90.inc)

  type, extends(tree) :: tree_expr
    type(vector_item) :: operands
    type(vector_char) :: operators
  end type
  type, extends(tree) :: tree_unit
    type(token)  :: tok
    type(string) :: prefix
    type(string) :: unit
  end type
  type, extends(tree) :: tree_int
    type(token) :: tok
    integer     :: value
  end type
  type, extends(tree) :: tree_ratio
    class(tree), allocatable :: nom
    class(tree), allocatable :: den
  end type
  type, extends(tree) :: tree_expo
    class(tree), allocatable :: base
    type(tree_ratio)         :: expo
  end type

  type lexer
    character(:), allocatable :: input
    integer                   :: pos
    integer                   :: read_pos
    character                 :: ch
  contains
    procedure :: init        => lexer_init
    procedure :: next_token  => lexer_next_token
    procedure :: read_ch     => lexer_read_ch
    procedure :: read_ident  => lexer_read_ident
    procedure :: read_number => lexer_read_number
    procedure :: peek_char   => lexer_peek_char
    procedure :: skip_ws     => lexer_skip_ws
  end type

  type parser
    type(lexer)         :: l
    type(vector_string) :: errors
    type(token)         :: curr_token
    type(token)         :: peek_token
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

contains

  m4_define({T},{item})
  m4_include(vector_imp.f90.inc)

  subroutine init_normconst(T)
    !! initialize global normalization object
    real, intent(in) :: T
      !! Temperature in Kelvin

    call normconst%init(T)
  end subroutine

  subroutine destruct_normconst()
    !! destruct global normalization object

    call normconst%destruct()
  end subroutine

  subroutine normalization_init(this, T)
    !! initialize normalization object
    class(normalization), intent(out) :: this
    real,                 intent(in)  :: T
      !! Temperature in Kelvin

    ! constants
    real, parameter :: EC = 1.602176634e-19
      !! elementary charge [ C ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?e
    real, parameter :: EM = 9.1093837015e-31
      !! electron rest mass [ kg ]
      !! rel err: 3e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?me
    real, parameter :: PLANCK = 6.62607015e-34
      !! Planck constant [ J/Hz ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?h
    real, parameter :: BOLTZ = 1.380649e-23
      !! Boltzmann constant [ J/K ]
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

    real :: ampere, coulomb, diel, farad, henry, hertz, kelvin, gram, meter, ohm, permby, second, volt, watt

    call this%prefixes%init()
    call this%units%init()
    call this%cache%init()

    ! metric prefixes
    call this%prefixes%insert(new_string("P"), 1e15)
    call this%prefixes%insert(new_string("T"), 1e12)
    call this%prefixes%insert(new_string("G"), 1e9)
    call this%prefixes%insert(new_string("M"), 1e6)
    call this%prefixes%insert(new_string("k"), 1e3)
    call this%prefixes%insert(new_string("h"), 1e2)
    call this%prefixes%insert(new_string("1"), 1e0)
    call this%prefixes%insert(new_string("d"), 1e-1)
    call this%prefixes%insert(new_string("c"), 1e-2)
    call this%prefixes%insert(new_string("m"), 1e-3)
    call this%prefixes%insert(new_string("u"), 1e-6)
    call this%prefixes%insert(new_string("n"), 1e-9)
    call this%prefixes%insert(new_string("p"), 1e-12)
    call this%prefixes%insert(new_string("f"), 1e-15)
    call this%prefixes%insert(new_string("a"), 1e-18)

    ! basic normalization constants
    volt    = EC / (BOLTZ * T)
    meter   = 2 * PI * sqrt(EM * BOLTZ * T) / PLANCK
    second  = 2 * PI * BOLTZ * T / PLANCK
    hertz   = 1.0 / second
    gram    = 1e-3 / EM ! grams instead of kilograms
    ampere  = 1.0 / (EC * second)
    coulomb = ampere * second
    farad   = coulomb / volt
    ohm     = volt / ampere
    henry   = ohm * second
    watt    = volt * ampere
    kelvin  = 1.0 / T
    diel    = (PLANCK * EPS0) / (2 * PI * EC * sqrt(EM * EC * volt)) ! eps0 is not normalized to 1.0
    permby  = MU0 * henry / meter

    ! basic units
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
    call this%units%insert(new_string("Å"), 1e-10 * meter)
    call this%units%insert(new_string("eps0"), diel)
    call this%units%insert(new_string("mu0"), permby)
    call this%units%insert(new_string("°"), PI / 180)
    call this%units%insert(new_string("deg"), PI / 180)
    call this%units%insert(new_string("rad"), 1.0)
  end subroutine

  function normalization_lookup(this, unit) result(val)
    !! lookup unit and get normalization constant
    class(normalization), intent(inout) :: this
    character(*),         intent(in)    :: unit
      !! physical unit token
    real                                :: val
      !! return normalization constant

    integer                            :: i
    type(string)                       :: sunit
    type(mapnode_string_real), pointer :: n
    type(lexer)                        :: l
    type(parser)                       :: p
    class(tree), allocatable           :: tr

    ! lookup unit in cache
    sunit%s = unit
    n => this%cache%find(sunit)
    if (associated(n)) then
      val = n%value
      return
    end if

    ! parse unit
    call l%init(unit)
    call p%init(l)
    call p%parse_expr(tr)

    ! check for errors
    if (p%errors%n > 0) then
      do i = 1, p%errors%n
        print "(A)", p%errors%d(i)%s
      end do
      call program_error("Cannot parse unit expression '" // unit // "'")
    end if

    ! return value and save in cache
    val = tr%eval()
    call this%cache%set(sunit, val)
  end function

  m4_define({m4_nvalues_shape},{m4_ifelse($1,1,size(values,1),{m4_nvalues_shape(m4_decr($1)),size(values,$1)})})
  m4_define({m4_nvalues_pshape},{m4_ifelse($1,0,,{(m4_nvalues_shape($1))})})
  m4_define({m4_X},{function $1_$2_$3(values, unit) result(nvalues)
    m4_norm_type($3),              intent(in) :: values{}m4_pshape($2)
      !! value to (de-)normalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    m4_norm_type($3)                          :: nvalues{}m4_nvalues_pshape($2)
      !! return normalized value

    nvalues = values m4_norm_op($1) normconst%lookup(unit)
  end function})
  m4_list

  subroutine normalization_destruct(this)
    !! destruct normalization constants (release memory)
    class(normalization), intent(inout) :: this

    call this%prefixes%destruct()
    call this%units%destruct()
    call this%cache%destruct()
  end subroutine

  pure function token_name(token_type) result(name)
    integer, intent(in) :: token_type
    character(:), allocatable :: name

    select case (token_type)
    case (TOKEN_EOF)
      name = "EOF"
    case (TOKEN_ILLEGAL)
      name = "ILLEGAL"
    case (TOKEN_IDENT)
      name = "IDENT"
    case (TOKEN_INT)
      name = "INT"
    case (TOKEN_ASTERISK)
      name = "ASTERISK"
    case (TOKEN_SLASH)
      name = "SLASH"
    case (TOKEN_CARET)
      name = "CARET"
    case (TOKEN_LPAREN)
      name = "LPAREN"
    case (TOKEN_RPAREN)
      name = "RPAREN"
    case default
      name = "OTHER"
    end select
  end function token_name

  subroutine lexer_init(this, input)
    class(lexer),  intent(inout) :: this
    character(*),  intent(in)    :: input

    this%input    = input
    this%pos      = 1
    this%read_pos = 1
    this%ch       = "x"
    call this%read_ch()
  end subroutine

  subroutine lexer_next_token(this, tok)
    class(lexer), intent(inout) :: this
    type(token),  intent(out)   :: tok

    character(:), allocatable :: lit

    call this%skip_ws()

    select case(this%ch)
    case ("/")
      tok = token(TOKEN_SLASH, this%ch)
    case ("*")
      tok = token(TOKEN_ASTERISK, this%ch)
    case ("^")
      tok = token(TOKEN_CARET, this%ch)
    case ("(")
      tok = token(TOKEN_LPAREN, this%ch)
    case (")")
      tok = token(TOKEN_RPAREN, this%ch)
    case (achar(0))
      tok = token(TOKEN_EOF, "")
    case ("-")
      if (is_digit(this%peek_char())) then
        call this%read_number(lit)
        tok = token(TOKEN_INT, lit)
        return
      else
        tok = token(TOKEN_ILLEGAL, this%ch)
      end if
    case default
      if (is_letter(this%ch)) then
        call this%read_ident(lit)
        tok = token(TOKEN_IDENT, lit)
        return
      else if (is_digit(this%ch)) then
        call this%read_number(lit)
        tok = token(TOKEN_INT, lit)
        return
      else
        tok = token(TOKEN_ILLEGAL, this%ch)
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
    this%pos      = this%read_pos
    this%read_pos = this%read_pos + 1
  end subroutine

  subroutine lexer_skip_ws(this)
    class(lexer), intent(inout) :: this

    do while (is_whitespace(this%ch))
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
    logical                   :: b

    b = (tok == this%peek_token%type)
  end function

  function parser_expect_peek(this, tok) result(b)
    class(parser), intent(inout) :: this
    integer,       intent(in)    :: tok

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

    if (.not. allocated(this%errors%d)) call this%errors%init(0, c = 4)

    call this%errors%push(new_string(msg))
  end subroutine

  subroutine parser_parse_expr(this, tr, nested)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: tr
    logical, optional,        intent(in)    :: nested

    type(tree_expr)     :: expr
    type(item)                :: it
    character(:), allocatable :: err
    character                 :: op
    logical                   :: nested_

    nested_ = .false.
    if (present(nested)) nested_ = nested

    tr = expr
    select type(tr)
    type is (tree_expr)
      call tr%operands%init(0, c = 4)
      call tr%operators%init(0, c = 4)

      do while(this%curr_token%type /= TOKEN_EOF .and. ((.not. nested_) .or. (this%curr_token%type /= TOKEN_RPAREN)))

        if ((this%curr_token%type == TOKEN_ASTERISK) .or. (this%curr_token%type == TOKEN_SLASH)) then
          op = this%curr_token%literal
          call this%next_token()
        else
          op = achar(0)
        end if

        if (tr%operands%n > 0 .and. op == achar(0)) then
          err = "Got no operator between values"
          call this%add_error(err)
          exit
        end if

        call this%parse_value(it%item)
        if (allocated(it%item)) then
          call tr%operands%push(it)
          if (op /= achar(0)) call tr%operators%push(op)
          deallocate(it%item)
        end if

        call this%next_token()
      end do

    end select
  end subroutine

  subroutine parser_parse_value(this, s)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: s

    character(:), allocatable :: err

    select case(this%curr_token%type)
    case (TOKEN_IDENT)
      call this%parse_identifier(s)
    case (TOKEN_INT)
      call this%parse_integer(s)
    case (TOKEN_LPAREN)
      call this%parse_group(s)
    case default
      err = "expected expression got " // token_name(this%curr_token%type)
      call this%add_error(err)
      return
    end select

    ! Possible exponentiation
    if (this%peek_token_is(TOKEN_CARET)) then
      call this%next_token()
      call this%parse_expo(s)
    end if
  end subroutine

  subroutine parser_parse_group(this, tr)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: tr

    character(:), allocatable :: err

    call this%next_token()
    call this%parse_expr(tr, nested=.true.)

    if (.not. this%curr_token%type == TOKEN_RPAREN) then
      deallocate(tr)
      err = "expected ')' after expression got " // token_name(this%curr_token%type)
      call this%add_error(err)
      return
    end if
  end subroutine

  subroutine parser_parse_expo(this, tr)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: tr

    type(tree_expo)  :: expo
    type(tree_ratio) :: ratio
    type(tree_int)   :: i

    character(:), allocatable  :: err
    logical                    :: b

    if (.not. allocated(tr)) then
      call program_error("Cannot parse exponent if no base is given")
    end if

    expo%base = tr
    deallocate(tr)

    call this%next_token()
    select case(this%curr_token%type)
    case (TOKEN_INT)
      call this%parse_integer(ratio%nom)

      i%tok = token(TOKEN_INT, "1")
      i%value = 1
      ratio%den = i

    case (TOKEN_LPAREN)
      b = this%expect_peek(TOKEN_INT)
      if (.not. b) return
      call this%parse_integer(ratio%nom)
      b = this%expect_peek(TOKEN_SLASH)
      if (.not. b) return
      b = this%expect_peek(TOKEN_INT)
      if (.not. b) return
      call this%parse_integer(ratio%den)
      b = this%expect_peek(TOKEN_RPAREN)
      if (.not. b) return

    case default
      err = "Expected integer or left parenthesis, got " // token_name(this%curr_token%type)
      call this%add_error(err)
      return

    end select

    expo%expo = ratio
    tr = expo
  end subroutine

  subroutine parser_parse_integer(this, tr)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: tr

    type(tree_int) :: s

    s%tok = this%curr_token
    read(this%curr_token%literal, *) s%value ! todo: handle read error
    tr = s
  end subroutine

  subroutine parser_parse_identifier(this, tr)
    class(parser),            intent(inout) :: this
    class(tree), allocatable, intent(inout) :: tr

    type(tree_unit)                    :: s
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
      tr = s
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
    tr = s
  end subroutine

  recursive function tree_eval(this) result(val)
    !! Serves as eval and destruct
    class(tree), intent(inout) :: this

    real    :: val
    integer :: i

    select type(this)
    type is (tree_expr)
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

    type is (tree_int)
      val = real(this%value)

    type is (tree_unit)
      val = normconst%prefixes%get(this%prefix) * normconst%units%get(this%unit)

    type is (tree_expo)
      val = this%base%eval() ** this%expo%eval()

    type is (tree_ratio)
      val = this%nom%eval() / this%den%eval()

    class default
      call program_error("Got undefined ast type")

    end select
  end function

  recursive function tree_to_string(this) result(s)
    !! convert to string for printing
    class(tree), intent(in)   :: this
    character(:), allocatable :: s

    integer :: i

    select type (this)
    type is (tree_expr)
      if (this%operands%n == 0) then
        s = "( )"
      else
        s = "( " // this%operands%d(1)%item%to_string()
        do i = 2, this%operands%n
          select case (this%operators%d(i-1))
          case ("*")
            s = s // " * " // this%operands%d(i)%item%to_string()
          case ("/")
            s = s // " / " // this%operands%d(i)%item%to_string()
          end select
        end do
        s = s // " )"
      end if

    type is (tree_unit)
      s = "unit[" // this%prefix%s // "-" // this%unit%s // "]"

    type is (tree_int)
      s = this%tok%literal

    type is (tree_expo)
      s = this%base%to_string() // "^" // this%expo%to_string()

    type is (tree_ratio)
      s = "(" // this%nom%to_string() // "/" // this%den%to_string() // ")"

    class default
      call program_error("unexpected type")

    end select
  end function

end module
