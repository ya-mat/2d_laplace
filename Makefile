LINKER	     = g++
CC	     = g++
LDFLAGS     =
COPTS	     = -cpp -O3
#COPTS	     = -cpp -O2 -fmax-errors=3 -g

FLINKER	     = gfortran
FORTRAN	     = gfortran
FLDFLAGS     = -llapack -lrefblas
#FLDFLAGS     = -lopenblas
FOPTS	     = -cpp -O3 -ffree-line-length-none -fmax-errors=3
#FOPTS	     = -cpp -g -O0 -fbounds-check -fmax-errors=3 -ffree-line-length-none -Wall

OBJS          = main.o\

FOBJS          = lp_laplace.o\
		force_raise.o\
		fmain.o\

PROGRAM      = a.out
FPROGRAM      = f.out

all:		$(PROGRAM) $(FPROGRAM)

$(PROGRAM): $(OBJS)
		$(LINKER) $(COPTS) $(OBJS) -o $(PROGRAM) $(DFLAGS)

$(FPROGRAM): $(FOBJS)
		$(FLINKER) $(FOPTS) $(FOBJS) -o $(FPROGRAM) $(FLDFLAGS)

clean:
		rm -f $(PROGRAM) *.o *~ *.mod;\

.SUFFIXES: .o .f90 .cc

.cc.o :
		$(CC) $(COPTS) -c -o $*.o $*.cc

.f90.o :
		$(FORTRAN) $(FOPTS) -c -o $*.o $*.f90
