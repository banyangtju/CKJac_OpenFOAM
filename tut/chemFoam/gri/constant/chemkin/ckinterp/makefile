##########################
### makefile in skle
##########################

FF = gfortran

FFLAGS = -c -o $@ -O3

OBJ = ckinterp.o

all : $(OBJ)
	$(FF) -O3 $(OBJ) -o ckinterp

####################################################
# Chemkin files
####################################################

.SUFFIXES: .F90 .f90

.f.o:
	$(FF) $(FFLAGS) $<

.f90.o:
	$(FF) $(FFLAGS) $<

#dvodeoriginal.o: dvodeoriginal.f
#	$(FF) $(FFLAGS) $<

cklib.o: cklib.f ckstrt.h
	$(FF) $(FFLAGS) $<

#main.o: main.f90
#	$(FF) $(FFLAGS) $<


cleanall:
	rm *.o run

clean:
	rm $(OBJ) run
