.PHONY : clean all chemcpd.exe

FLAG = -fbounds-check -fbacktrace -g
OBJ = quadpack.o minpack.o constant_table.o particular_functions.o profiles.o  
DIR = ./src

all: chemcpd.exe

chemcpd.exe : $(DIR)/main.f95  quadpack.o minpack.o constant_table.o env.o particular_functions.o profiles.o 
	gfortran $(FLAG) $(DIR)/main.f95 quadpack.o minpack.o constant_table.o profiles.o particular_functions.o env.o -o chemcpd.exe 
	make clean
	make reset


minpack.o : $(DIR)/minpack.f95 
	gfortran -c $(DIR)/minpack.f95

quadpack.o : $(DIR)/quadpack.f95 
	gfortran -c $(DIR)/quadpack.f95

profiles.o : $(DIR)/profiles.f95 
	gfortran -c $(DIR)/profiles.f95

constant_table.o : $(DIR)/constant_table.f95 
	gfortran -c $(DIR)/constant_table.f95

particular_functions.o : $(DIR)/particular_functions.f95
	gfortran -c $(DIR)/particular_functions.f95

env.o : $(DIR)/env.f95
	gfortran -c $(DIR)/env.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ../Data/*