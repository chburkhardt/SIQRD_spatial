# SIQRD_spatial

The code presented here belongs to the paper "Meso-scale modeling of COVID-19 
spatio-temporal outbreak dynamics in Germany" (DOI: 10.1101/2020.06.10.20126771)
by Andreas Kergassner, Christian Burkhardt, Dorothee Lippold, Sarah Nistler, 
Matthias Kergassner, Paul Steinmann, Dominik Budday, Silvia Budday

Due to the dynamic development process of the implementation, the names of the 
compartments in the paper differ from those in the implementation.

		paper 			--	implementation

		Susceptibles	S 	--	S
		Quarantined 	Q 	--	I
		Recovered 	R 	--	R
		Infected 	I 	--	E
		Dead 		D	--	D

---------------------------------------------------------------------
## Structure

### Daten
* population.txt| population_bayern.txt contain the data from the Statistischen Bundesamt for counties with coordinated, inhabitants,size etc.
* RKI_Corona_Landkreise.txt contains data from RKI with countywise infection numbers
	
### Results
* This is where all results are written in, namely vtks, the parameter file used, *.mat files with a struct for each county and its solution

### Code
is structured in a spatial, a local and a particle swarm optimization (PSO) part

NOTE:
vtkwrite.m is based on Chaoyuan Yeh, 2016  %  Version 2.3 and  
ode23d_fixed.m is based on the corresponding ode23d.m file from the OdePkg package of GNU Octave and modified


#### Local
* SIR_Run.m is the mainfile
* siredMod.m is the file where the ODE is defined and further files needed for the parameter fitting

#### Spatial
* sir_spatial.m is the mainfile
* read_*.m and setup_system.m files input and build all necessary data 
* sir_eqn_spatial.m is the spatial extended SIQRDd Model
* standAlonePosteProcessing.m, generateVTK.m and vtkwrite.m are for postprocessing

#### PSO
* defines the routine for Particel Swarm Optimization
	




