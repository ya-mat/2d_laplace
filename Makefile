LINKER	     = gfortran
FORTRAN	     = gfortran
LDFLAGS	     = -L/usr/lib/ -lblas -llapack
FOPTS	     = -cpp -O3 -ffree-line-length-none -fmax-errors=3

OBJS          = lp_laplace.o\
		force_raise.o\
		main.o\

PROGRAM	      = a.out

all:		$(PROGRAM)

$(PROGRAM): $(OBJS)
		$(LINKER) $(FOPTS) $(OBJS) -o $(PROGRAM) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) *.o *~ *.mod;\

.SUFFIXES: .o .f90

.f90.o :
		$(FORTRAN) $(FOPTS) -c -o $*.o $*.f90

