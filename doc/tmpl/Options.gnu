# compiler
F90      = ftn -D__GNU_FORTRAN__

# linker
LD       = $(F90)

# global compilation flags
FFLAGS   = -fno-fast-math # enforce IEEE behaviour
FFLAGS  += -ffree-form # specify free-format source
FFLAGS  += -ffree-line-length-none # specify infinitely long lines
FFLAGS  += -fno-backslash # treat backslash as normal character
FFLAGS  += -fimplicit-none # no implicit typing is allowed
FFLAGS  += -frange-check # range checking of constants at compile time
#FFLAGS  += -std=f95 # conform to Fortran 95 standard
#FFLAGS  += -std=f2003 # conform to Fortran 2003 standard
FFLAGS  += -pedantic # warn about extensions to Fortran 95 standard

# verbose compilation flags
VFLAGS   = -v # displays invocations of compiler, assembler and linker
VFLAGS  += -Wall # enable commonly used warning options
VFLAGS  += -Wconversion # warn about implicit conversion of types

# force preprocessing flags
PFLAGS   = -cpp # perform preprocessing for all source files

# optimization flags (FOPT0=none, FOPT1=strong, FOPT2=strongest)
FOPT0    = -O0 # inhibit any optimizations

FOPT1    = -O3 # aggressive global optimization
FOPT1   += -ftree-vectorize #
#FOPT1   += -ftree-loop-linear #
FOPT1   += -funroll-loops # 

FOPT2    = $(FOPT1)

# debugging flags
FDBG     = -O0 # inhibit any optimizations
FDBG    += -g # generate information for debugger (disable optimizations)
FDBG    += -fbacktrace # provide backtrace of error upon crash
FDBG    += -fdump-core # dump core file on crash
FDBG    += -ffpe-trap=invalid,zero,overflow # floating point exceptions (invalid,zero,underflow)
FDBG    += -fbounds-check # check array bounds
FDBG    += -fcheck-array-temporaries # warn about array temporaries (at runtime)
FDBG    += -finit-real=nan # initialize all reals with NaN
FDBG    += -finit-integer=999999 # initialize ints with a rediculous value
FDBG    += -finit-logical=true # un-intuitive logical initialization
FDBG    += -finit-character=35 # un-intuitive character initialization

# global linking flags
LFLAGS   = 

# global libraries and includes
LIB      = 

INC      = -I.
INC     += -I$(ROOT)/$(SRCDIR)
INC     += -I$(MPICH_DIR)/include

# optimized flags, libraries and includes
OPTL     = 
OPTI     = 

# debugging flags, libraries and includes
DBGL     = 
DBGI     = 

# MPI library
MPIL     = 
MPII     = -D__MPICH2

# Grib1 library
GRIBDWDL = 
GRIBDWDI = 

# Grib-API library
GRIBAPIL =
GRIBAPII =

# NetCDF library
NETCDFL  = 
NETCDFI  = 

# Synsat library
RTTOV7L  = 
RTTOV7I  = 
RTTOV10L = 
RTTOV10I = 

