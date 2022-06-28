# Chemcpd

2D (r,z) Model evolution of Jupiter's CPD. Implementing chemical reaction chains, gas  evolution based on Makalkin & Dorofeeva 2014. 
The model will proceed in three steps:
- [x] Calcul of the midplane temperature as function of the disk mass, accretion rate, planetary luminosity and planetray mass 
- [x] Calcul the entire 2D (r,z) profile with the midplane temperature and surface temperature
- [ ] In this 2D profile we add the evolution of dust and chemical species 

Each step can be precomputed independantly from the others 

The data explorer system is present in the folder `posprocessing`. It is composed of a parser and a plotter to explore the solution and intermediary files. 

The physical explanation is in the wiki of the repo. 
