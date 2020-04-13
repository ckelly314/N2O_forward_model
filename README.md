# Forward-running ETNP N2O box model
Requires Python 3.7 and Jupyter 4.4 or above
@author: Colette Kelly,
Stanford University

### DESCRIPTION:
    Simple time-varying model to model the expected d15N2O-alpha, d15N2O-beta, and d18O-N2O from different
    reductive substrates (nitrate and nitrite). Utilizes the known isotope values of the substrate and a
    range of expected isotope effects to model the expected isotopes of N2O.
    
    Model the time-dependent change in N2O isotopes from different substrate pools.
    Created on Mon October 28th, 2019
    
### STATE VARIABLES:
	- n2o_14 = concentration of 14N-14N-16O, umol N/L
	- n2o_15A = concentration of 14N-15N-16O, umol N/L
	- n2o_15B = concentration of 15N-14N-16O, umol N/L
	- n2o_16 = concentration of 14N-14N-16O, umol N/L (i.e., 1/2 O per every N2O)
	- n2o_18 = concentration of 14N-14N-18O, umol N/L (i.e., 1/2 O per every N2O)

### (POTENTIALLY) TIME-VARYING PARAMETERS:
	- Rate constants for production and consumption of N2O
	- Delta values of substrates (nitrate and nitrite)

### OUTPUT:
	- Value of each state variable at each time step.

@author: Colette Kelly, Stanford University (clkelly@stanford.edu)

### Running the model
#### Open the Jupyter Notebook "200413 forward model v4.ipynb"
    Run the model by running the cells in the Jupyter notebook.
       
 #### convert_delta.py
      Convert isotope ratio in delta notation to 15R or 18R

 #### model_functions.py
      Functions to initialize and run model, as well as function for time-varying
      parameters.
       
 #### N2O_forward_model_v4_/Figures
 	  Where model output figures are stored.