all: chemcpd.exe

chemcpd.exe : main.f95  quadpack.o minpack.o constant_table.o profiles.o  
	gfortran  main.f95 quadpack.o minpack.o constant_table.o profiles.o  -o chemcpd.exe 
	make clean

.PHONY : clean

minpack.o : minpack.f95 
	gfortran -c minpack.f95

quadpack.o : quadpack.f95 
	gfortran -c quadpack.f95

profiles.o : profiles.f95 
	gfortran -c profiles.f95

constant_table.o : constant_table.f95 
	gfortran -c constant_table.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe