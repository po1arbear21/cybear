# build configuration (debug, release or profile)
BUILD = debug

# compiler (intel or gnu)
COMPILER = intel

# libraries (disable unused, ILUPACK and MUMPS are incompatible)
USE_ARPACK   = true
USE_EXPOKIT  = true
USE_FEAST    = true
#USE_ILUPACK  = true
USE_MUMPS    = true
USE_QUADPACK = true
