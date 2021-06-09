#=======================================================================
# User definitions
#

NAME = test

F90FILES = control.f90 ic.f90 rnorm.f90
CFILES = fileio_c.c

LIBS  =
FLAGS =

#=======================================================================
# Standard definitions
#
OFILES = $(F90FILES:.f90=.o) $(CFILES:.c=.o)

# SGI compilers
F90 = gfortran -ffree-line-length-none
LD  = gfortran -ffree-line-length-none
CC  = gcc

VPATH = obj
.SUFFIXES:
.SUFFIXES: .out .o .f90 .f .c .inc .h

#=======================================================================
# Targets and dependencies
#
default:: opt
all:: a.out
debug ::
	@make all "FLAGS = -g -CB -CV -traceback -ftrapuv -CS -check all -O0 $(FLAGS)"
prof::
	@make all "LIBS = -p $(LIBS)"
opt::
	@make all "FLAGS = -O2  $(FLAGS)"
clean::
	rm -rf obj a.out *.kmo *.mod $(NAME)

a.out: obj $(OFILES)
	@echo "$(LD) $(OFILES) $(LIBS)"
	@cd obj; $(LD) $(OFILES) $(LIBS); mv a.out ../$(NAME)
obj:
	mkdir obj

#=======================================================================
# Compilation rules
#
.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o obj/$*.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o obj/$*.o
