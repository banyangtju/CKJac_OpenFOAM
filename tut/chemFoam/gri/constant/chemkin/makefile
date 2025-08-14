##########################
### makefile in skle
##########################

FF = gfortran
CC = gcc
CLFLAGS = -shared -fPIC -o ./libchemkin.so -Wl,--no-undefined -lpthread

FFLAGS = -c -o $@ -O3 -fPIC

OBJ = out.o

all : $(OBJ)
	g++ $(CLFLAGS) $(OBJ) -lm -lgfortran
#	$(CC) $(CLFLAGS) $(OBJ) -L . $(dylib) -lm -lgfortran

####################################################
# ZDPlasKin files
####################################################

.SUFFIXES: .F90 .f90

.F90.o:
	$(FF) $(FFLAGS) $<
.f90.o:
	$(FF) $(FFLAGS) $<

#dvode_f90_m.o: dvode_f90_m.F90
#	$(FF) $(FFLAGS) $<
#dvode_f90_m.mod : dvode_f90_m.o
#	@echo Compiling .\dvode_f90_m.F90
#	@$(FF) $(FFLAGS) dvode_f90_m.o 

#zdplaskin_m.o: zdplaskin_m.F90 dvode_f90_m.mod
#	$(FF) $(FFLAGS) $<
#zdplaskin_m.mod : zdplaskin_m.o
#	@echo Compiling .\zdplaskin_m.F90
#	@$(FF) $(FFLAGS) zdplaskin_m.o 

#callZD.o: callZD.f90
#	$(FF) $(FFLAGS) $<


cleanall:
	rm *.o *.mod run

clean:
	rm $(OBJ) *.mod run
