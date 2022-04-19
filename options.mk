# build configuration (debug, release or profile)
BUILD = debug

# compiler (intel or gnu)
COMPILER = gnu

# default integer size (32 or 64)
INTSIZE = 32

# sparse matrix index integer size (32 or 64; use 64 if you are using huge matrices; FIXME: 64 broken for mkl_ilu)
IDXSIZE = 32

# libraries (enable used libraries, incompatibility between ILUPACK and MUMPS)
# USE_ARPACK   = true
# USE_EXPOKIT  = true
# USE_FEAST    = true
# USE_ILUPACK  = true
# USE_MUMPS    = true
# USE_QUADPACK = true

# python 3
PYTHON = python3
