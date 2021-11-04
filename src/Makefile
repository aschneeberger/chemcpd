all: chemcpd.exe

chemcpd.exe : main.o residue.o solver.o 
	gfortran -o chemcpd.exe solver.o residue.o 

main.o : main.f95 
	gfortran -c main.f95 