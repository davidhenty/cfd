# System dependent definitions

# For ARCHER
FC=	ftn
CFLAGS=
LFLAGS=

# For Cirrus or any standard Linux system
#FC=	ifort
#FFLAGS= -qopenmp -O3
#LFLAGS=

# System independent definitions

MF=	Makefile

EXE=	cfd

SRC= \
	boundary.f90 \
	cfd.f90 \
	cfdio.f90 \
	jacobi.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

cfd.o:	boundary.o jacobi.o cfdio.o

tar:
	tar cvf cfd.tar $(MF) $(INC) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) *.mod velocity.dat colourmap.dat cfd.plt core
