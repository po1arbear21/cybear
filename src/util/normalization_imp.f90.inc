#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#ifdef NORM
#define OPSIGN               /
#define OPNAME               norm
#endif

#ifdef DENORM
#define OPSIGN               *
#define OPNAME               denorm
#endif

#define OP_SCALAR_TYPE       CONCAT3(OPNAME,scalar,T)
#define OP_ARRAY_TYPE        CONCAT3(OPNAME,array,T)
#define OP_MATRIX_TYPE       CONCAT3(OPNAME,matrix,T)

function OP_SCALAR_TYPE(value, unit, n) result(nvalue)
  !! (de-)normalize scalar

  TT,                            intent(in) :: value
    !! value to normalize
  character(*),                  intent(in) :: unit
    !! physical unit token
  type(normalization), optional, intent(in) :: n
    !! optional normalization object (default: use global normconst)
  TT                                        :: nvalue
    !! return normalized value

  nvalue = value OPSIGN get_norm_value(unit, n)
end function

function OP_ARRAY_TYPE(values, unit, n) result(nvalues)
  !! (de-)normalize array

  TT,                            intent(in) :: values(:)
    !! values to normalize
  character(*),                  intent(in) :: unit
    !! physical unit token
  type(normalization), optional, intent(in) :: n
    !! optional normalization object (default: use global normconst)
  TT                                        :: nvalues(size(values))
    !! return normalized values

  nvalues = values OPSIGN get_norm_value(unit, n)
end function

function OP_MATRIX_TYPE(values, unit, n) result(nvalues)
  !! (de-)normalize matrix

  TT,                            intent(in) :: values(:,:)
    !! values to normalize
  character(*),                  intent(in) :: unit
    !! physical unit token
  type(normalization), optional, intent(in) :: n
    !! optional normalization object (default: use global normconst)
  TT                                        :: nvalues(size(values,1),size(values,2))
    !! return normalized values

  nvalues = values OPSIGN get_norm_value(unit, n)
end function

#undef T
#undef TT
#undef NORM
#undef DENORM

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef OPSIGN
#undef OPNAME
#undef OP_SCALAR_TYPE
#undef OP_ARRAY_TYPE
#undef OP_MATRIX_TYPE
