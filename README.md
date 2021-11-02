# Chemcpd

2D (r,z) Model evolution of Jupiter's CPD. Implementing chemical reaction chains, gas  evolution based on Makalkin & Dorofeeva 2014. 
The model will proceed in three steps:
- [ ] Calcul of the midplane temperature as function of the disk mass, accretion rate, planetary luminosity and planetray mass 
- [ ] Calcul the entire 2D (r,z) profile with the midplane temperature and surface temperature
- [ ] In this 2D profile we add the evolution of dust and chemical species 

Each step can be precomputed independantly from the others 

Currently the model is tested with the python file `poc.py` and its test unit tests `test_poc.py` 
