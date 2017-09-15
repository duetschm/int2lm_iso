# Makefile for compiling COSMO
# fuo, 07.03.2008, oliver.fuhrer@meteoswiss.ch
#
# Description: This makefile automagically searches for source files in the src
#              directory and creates dependencies. Afterwards files are compiled
#              and linked together.
#              Suitable for both, cosmo and int2lm
# History:
#   fuo  07.03.2008   first release


### Macros ###########################

# define shell
SHELL      := /bin/sh

# define the name of the target
TARGET     := int2lm
 
# directory definitions
SRCDIR     := src
OBJDIR     := obj
ifeq (0,${MAKELEVEL})
  ROOT      := $(shell pwd)
  VPATH     := .:$(ROOT)/$(SRCDIR)
endif

# generate list of object files
-include $(ROOT)/Objfiles
ifndef OBJ
  $(error Could not load Objfiles file)
endif

# dynamically generated dependency file
DEPF       := .depend
IGN        := --ignore netcdf --ignore grib_api --ignore oas_cos_vardef

# select machine dependent stuff
-include $(ROOT)/Options
ifndef F90
  $(error Must create platform specific Options file)
endif
FFLAGS0 = $(FFLAGS)
FFLAGS1 = $(FFLAGS)
FFLAGS2 = $(FFLAGS)

# handle options
ifdef VERBOSE
  FFLAGS0  += $(VFLAGS)
  FFLAGS1  += $(VFLAGS)
  FFLAGS2  += $(VFLAGS)
endif
ifdef DEBUG
  FFLAGS0  += $(FDBG)
  FFLAGS1  += $(FDBG)
  FFLAGS2  += $(FDBG)
  LIB      += $(DBGL)
  INC      += $(DBGI)
endif
ifdef OPT
  FFLAGS0  += $(FOPT0)
  FFLAGS1  += $(FOPT1)
  FFLAGS2  += $(FOPT2)
  LIB      += $(OPTL)
  INC      += $(OPTI)
endif
ifdef MPI
  LIB      += $(MPIL)
  INC      += $(MPII)
else
  OBJ      += dummy_mpi.o
  INC      += -I$(ROOT)/$(SRCDIR)/mpi
endif
ifdef DB
  LIB      += $(DBL)
  INC      += $(DBI)
else
  OBJ      += dummy_db.o
endif
ifdef GRIBDWD
  PFLAGS   += -DGRIBDWD
  LIB      += $(GRIBDWDL)
  INC      += $(GRIBDWDI)
endif
ifdef GRIBAPI
  PFLAGS   += -DGRIBAPI
  LIB      += $(GRIBAPIL)
  INC      += $(GRIBAPII)
  IGN      += --ignore grib_api
else
  IGN      += --ignore grib_api
endif
ifdef NETCDF
  PFLAGS   += -DNETCDF
  LIB      += $(NETCDFL)
  INC      += $(NETCDFI)
  IGN      += --ignore netcdf
else
  IGN      += --ignore netcdf
endif
export ROOT VPATH VERBOSE DEBUG OPT MPI DB GRIBDWD GRIBAPI NETCDF

### Phony targets ###########################

.PHONY : default depend paropt pardebug seqopt seqdebug clean info

default : paropt

depend :
	@echo "generating dependencies"
	@$(ROOT)/bin/sfmakedepend --case down --longpath $(INC) $(IGN) --file $(ROOT)/$(DEPF) $(ROOT)/$(SRCDIR)/*.f90

info :
	@echo "generating compile information"
	@-rm -f $(ROOT)/.fconfig
	@echo "Target           : $(TARGET)" > $(ROOT)/.fconfig
	@echo "Compiler command : $(F90)" > $(ROOT)/.fconfig
	@echo "Compiler version : "`$(F90) -V 2>/dev/null | grep pgf` >> $(ROOT)/.fconfig
	@echo "Compiler includes: $(INC)" >> $(ROOT)/.fconfig
	@echo "Compiler flags   : $(PFLAGS) $(FFLAGS1)" >> $(ROOT)/.fconfig
	@echo "Linker command   : $(LD)" >> $(ROOT)/.fconfig
	@echo "Linker version   : "`$(LD) -V 2>/dev/null | grep pgf` >> $(ROOT)/.fconfig
	@echo "Linker flags     : $(LFLAGS) $(FFLAGS1)" >> $(ROOT)/.fconfig
	@echo "Linker libraries : $(LIB)" >> $(ROOT)/.fconfig
	@$(ROOT)/bin/gen_info.sh $(ROOT)/.fconfig $(ROOT)/$(SRCDIR)
	@-rm -f $(ROOT)/.fconfig

paropt :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 MPI=1 OPT=1 depend info
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 MPI=1 OPT=1 $(ROOT)/$(TARGET)

pardebug :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 MPI=1 DEBUG=1 depend info
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 MPI=1 DEBUG=1 $(ROOT)/$(TARGET)

seqopt :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 OPT=1 depend info
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 OPT=1 $(ROOT)/$(TARGET)

seqdebug :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 DEBUG=1 depend info
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile GRIBDWD=1 NETCDF=1 DEBUG=1 $(ROOT)/$(TARGET)

$(TARGET) :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile depend info
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile $(ROOT)/$(TARGET)

clean :
	-rm -f $(DEPF) $(DEPF).old $(OBJDIR)/*

$(ROOT)/$(TARGET) : $(OBJ)
	@echo "linking $@"
	@$(LD) $(LFLAGS) $(FFLAGS1) $(INC) $(OBJ) $(LIB) -o $@

### Suffix Rules ###########################

.SUFFIXES: .o .mod .f90

# standard suffix rules (with optimization)
%.o : %.f90
	@echo "compiling $(patsubst %.o,%.f90,$(notdir $@))"
	@$(F90) -c $(PFLAGS) $(FFLAGS1) $(INC) -o $@ $(ROOT)/$(SRCDIR)/$(patsubst %.o,%.f90,$(notdir $@))

# include dynamically generated dependency file
-include $(ROOT)/$(DEPF)

# goodbye earthling!
