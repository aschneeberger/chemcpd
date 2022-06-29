.PHONY : clean all chemcpd.exe


# Python tool to generate dependencies
MAKEDEPEND = ~/.local/bin/fortdepend
# File where all dependencies are stored 
DEP_FILE = chemcpd.dep

OBJECTS  = src/*.f90 #src/ACE/*.f90

DBFLAG = -fbounds-check -fbacktrace  -g3 -ffpe-trap=invalid,zero,overflow  #-Wall -pedantic
OFLAG =  -O3 -ffast-math -flto -march=native -funroll-loops -fallow-store-data-races #-fopenmp
CFLAG =  $(DBFLAG) $(OFLAG)
DFLAG = -c $(CFLAG)


# MAke sure everything depends on the .dep file 
all: chemcpd.exe $(DEP_FILE)

# Make dependencies
.PHONY. : depend

depend : $(DEP_FILE)

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(DEP_FILE) : $(OBJECTS)
	@echo "Making dependencies"
	$(MAKEDEPEND) -w -o $(DEP_FILE) -f $(OBJECTS)

# Build all dependencies 
.f90.o : 
	gfortran -c $(CFLAG) $(OBJECTS)

# Build the executable 
chemcpd.exe :  $(OBJECTS)
	gfortran $(CFLAG) $^ -o  $@
	make clean 

# remove all mod files 
clean : 
	rm *.mod

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ./Data/*


include $(DEP_FILE)