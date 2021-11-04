all: chemcpd.exe

chemcpd.exe : main.f95 residue.o solver.o 
	gfortran  main.f95 solver.o residue.o -o chemcpd.exe 
	make clean

.PHONY : clean

solver.o : solver.f95 
	gfortran -c solver.f95

residue.o : residue.f95 
	gfortran -c residue.f95

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe