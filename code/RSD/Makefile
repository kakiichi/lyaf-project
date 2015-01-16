# For GCC Suite
FC = gfortran
#FC = ifort
FFLAGS = -fdefault-real-8
#FFLAGS = -g -pg -Wall -fbounds-check
#FFLAGS = -pg
#FFLAGS = -O3


SRC = param.f90 \
      RSD.f90 \
      quadpack.f90

OBJS = $(SRC:.f90=.o)

RSD: $(OBJS) 
	  $(FC) $(FFLAGS) -o RSD.exe *.o

%.mod: %.o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	-rm $(OBJS)
	-rm *~                                                                       

