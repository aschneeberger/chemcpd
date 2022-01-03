.PHONY : clean all chemcpd.exe

CFLAG = -fbounds-check -fbacktrace -g3 -ffpe-trap=invalid,zero,overflow -Wall -o3
OBJ = quadpack.o minpack.o constant_table.o particular_functions.o profiles.o  
DIR = ./src

all: chemcpd.exe

chemcpd.exe : $(DIR)/main.f95  quadpack.o minpack.o constant_table.o env.o particular_functions.o profiles.o 
	gfortran $(CFLAG) $(DIR)/main.f95 quadpack.o minpack.o constant_table.o profiles.o particular_functions.o env.o -o chemcpd.exe 
	make clean
	make reset


minpack.o : $(DIR)/minpack.f95 
	gfortran -c -g3 $(DIR)/minpack.f95

quadpack.o : $(DIR)/quadpack.f95 
	gfortran -c -g3 $(DIR)/quadpack.f95

profiles.o : $(DIR)/profiles.f95 
	gfortran -c -g3 $(DIR)/profiles.f95

constant_table.o : $(DIR)/constant_table.f95 
	gfortran -c -g3 $(DIR)/constant_table.f95

particular_functions.o : $(DIR)/particular_functions.f95
	gfortran -c -g3 $(DIR)/particular_functions.f95

env.o : $(DIR)/env.f95
	gfortran -c -g3 $(DIR)/env.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ./Data/*