# .PHONY : clean all chemcpd.exe


# # Python tool to generate dependencies
# MAKEDEPEND = ~/.local/bin/fortdepend
# # File where all dependencies are stored 
# DEP_FILE = chemcpd.dep

# OBJECTS  = src/*.f90 src/ACE/*.f90

# DBFLAG = -fbounds-check -fbacktrace  -g3 -ffpe-trap=invalid,zero,overflow  #-Wall -pedantic
# OFLAG =  -O3 -ffast-math -flto -march=native -funroll-loops -fallow-store-data-races #-fopenmp
# CFLAG =  $(DBFLAG) $(OFLAG)
# DFLAG = -c $(CFLAG)


# # MAke sure everything depends on the .dep file 
# all: chemcpd.exe $(DEP_FILE)

# # Make dependencies
# .PHONY. : depend

# depend : chemcpd.exe $(DEP_FILE)

# # The .dep file depends on the source files, so it automatically gets updated
# # when you change your source
# $(DEP_FILE) : $(OBJECTS)
# 	@echo "Making dependencies"
# 	$(MAKEDEPEND) -w -o $(DEP_FILE) -f $(OBJECTS)

# # Build all dependencies 
# include $(DEP_FILE)


# # Build the executable 
# chemcpd.exe :  $(OBJECTS)
# 	gfortran -N $(CFLAG) $^ -o  $@
# 	make clean 

# # remove all mod files 
# clean : 
# 	-rm *.mod

# mrproper : clean 
# 	-rm chemcpd.exe

# reset : 
# 	-rm -rf ./Data/*


.PHONY : clean all chemcpd.exe
DBFLAG = -fbounds-check -fbacktrace -Wall -pedantic -g3 -ffpe-trap=invalid,zero,overflow 
OFLAG =  -O3 -ffast-math -flto -march=native -funroll-loops -fallow-store-data-races #-fopenmp
CFLAG =  $(DBFLAG)# $(OFLAG)
DFLAG = -c $(CFLAG)

DIR = ./src

all: chemcpd.exe

chemcpd.exe : $(DIR)/main.f90  quadpack.o minpack.o constant_table.o env.o particular_functions.o profiles.o solver.o Md_Types_Numeriques.o Md_Constantes.o Md_numerical_recipes.o Md_Utilitaires.o Md_parametres.o  Md_ACE.o 
	gfortran $(CFLAG) $(DIR)/main.f90 quadpack.o minpack.o constant_table.o env.o particular_functions.o profiles.o solver.o Md_Types_Numeriques.o Md_Constantes.o Md_numerical_recipes.o Md_Utilitaires.o Md_parametres.o  Md_ACE.o -o chemcpd.exe  
	make clean
	make reset

minpack.o : $(DIR)/minpack.f90 
	gfortran $(DFLAG) $(DIR)/minpack.f90

quadpack.o : $(DIR)/quadpack.f90 
	gfortran $(DFLAG) $(DIR)/quadpack.f90

profiles.o : $(DIR)/profiles.f90 
	gfortran $(DFLAG) $(DIR)/profiles.f90

constant_table.o : $(DIR)/constant_table.f90 
	gfortran $(DFLAG) $(DIR)/constant_table.f90

particular_functions.o : $(DIR)/particular_functions.f90
	gfortran $(DFLAG) $(DIR)/particular_functions.f90

env.o : $(DIR)/env.f90
	gfortran $(DFLAG) $(DIR)/env.f90

solver.o : $(DIR)/solver.f90
	gfortran $(DFLAG) $(DIR)/solver.f90

Md_ACE.o : $(DIR)/ACE/Md_ACE.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_ACE.f90

Md_Constantes.o : $(DIR)/ACE/Md_Constantes.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_Constantes.f90

Md_numerical_recipes.o : $(DIR)/ACE/Md_numerical_recipes.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_numerical_recipes.f90

Md_parametres.o : $(DIR)/ACE/Md_parametres.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_parametres.f90

Md_Types_Numeriques.o : $(DIR)/ACE/Md_Types_Numeriques.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_Types_Numeriques.f90

Md_Utilitaires.o : $(DIR)/ACE/Md_Utilitaires.f90
	gfortran $(DFLAG) $(DIR)/ACE/Md_Utilitaires.f90

clean : 
	rm *.o *.mod 

mrproper : clean 
	rm chemcpd.exe

reset : 
	rm -rf ./Data/*


