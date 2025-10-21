FC = ifx

CMPLFLAGS = -132 -O3 #-xW -tpp7

LIB = -mkl

LIBRARIES = $(LIB)

SOURCE = MFGDy2.F90 LIB.F system.f90 mainMFG.f90

OBJECT=$(SOURCE:.f=.o)

.SUFFIXES: .f90 .f .o

exe: $(OBJECT) 
	$(FC) $(CMPLFLAGS) -o exe $(OBJECT) $(LIBRARIES)

.f90.o:
	$(FC) $(CMPLFLAGS) -c $<
.f.o:
	$(FC) $(CMPLFLAGS)  -c $<

clean:
	rm -f exe *.mod *.o
