all: chemcpd.exe

chemcpd.exe : main.f95  minipack.o residue.o solver.o 
	gfortran  main.f95 minipack.o solver.o residue.o -o chemcpd.exe 
	make clean

.PHONY : clean

solver.o : solver.f95 
	gfortran -c solver.f95

minipack.o : minipack.f95 
	gfortran -c minipack.f95

residue.o : residue.f95 
	gfortran -c residue.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe