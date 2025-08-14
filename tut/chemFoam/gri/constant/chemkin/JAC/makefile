##########################
### makefile in skle
##########################

FF = gfortran

FFLAGS = -c -o $@ -O3

OBJ = JAC.o cklib.o

all : $(OBJ)
	$(FF) -O3 $(OBJ) -o run

####################################################
# Chemkin files
####################################################

JAC.o: JAC.F90
	$(FF) $(FFLAGS) $<

cklib.o: cklib.f ckstrt.h
	$(FF) $(FFLAGS) $<


cleanall:
	rm *.o *.mod run

clean:
	rm $(OBJ) run