.PHONY : clean all chemcpd.exe

FLAG = -O2
OBJ = quadpack.o minpack.o constant_table.o particular_functions.o profiles.o  

all: chemcpd.exe

chemcpd.exe : main.f95  quadpack.o minpack.o constant_table.o particular_functions.o profiles.o  
	gfortran $(FLAG) main.f95 quadpack.o minpack.o constant_table.o profiles.o particular_functions.o -o chemcpd.exe 
	make clean
	make reset


minpack.o : minpack.f95 
	gfortran -c minpack.f95

quadpack.o : quadpack.f95 
	gfortran -c quadpack.f95

profiles.o : profiles.f95 
	gfortran -c profiles.f95

constant_table.o : constant_table.f95 
	gfortran -c constant_table.f95

particular_functions.o : particular_functions.f95
	gfortran -c particular_functions.f95

clean : 
	rm *.o *.mod *.dat

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ../Data/*