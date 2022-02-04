.PHONY : clean all chemcpd.exe
#-Ofast 
CFLAG = -fbounds-check -fbacktrace -Wall -pedantic -g3 -ffpe-trap=invalid,zero,overflow 
DFLAG = -c -fbounds-check -fbacktrace -Wall -pedantic  -g3 -ffpe-trap=invalid,zero,overflow
OBJ = quadpack.o minpack.o constant_table.o particular_functions.o profiles.o  
DIR = ./src

all: chemcpd.exe

chemcpd.exe : $(DIR)/main.f95  quadpack.o minpack.o constant_table.o env.o particular_functions.o profiles.o solver.o 
	gfortran $(CFLAG) $(DIR)/main.f95 quadpack.o minpack.o constant_table.o profiles.o particular_functions.o env.o solver.o -o chemcpd.exe  
	make clean
	make reset


minpack.o : $(DIR)/minpack.f95 
	gfortran $(DFLAG) $(DIR)/minpack.f95

quadpack.o : $(DIR)/quadpack.f95 
	gfortran $(DFLAG) $(DIR)/quadpack.f95

profiles.o : $(DIR)/profiles.f95 
	gfortran $(DFLAG) $(DIR)/profiles.f95

constant_table.o : $(DIR)/constant_table.f95 
	gfortran $(DFLAG) $(DIR)/constant_table.f95

particular_functions.o : $(DIR)/particular_functions.f95
	gfortran $(DFLAG) $(DIR)/particular_functions.f95

env.o : $(DIR)/env.f95
	gfortran $(DFLAG) $(DIR)/env.f95

solver.o : $(DIR)/solver.f95
	gfortran $(DFLAG) $(DIR)/solver.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ./Data/*