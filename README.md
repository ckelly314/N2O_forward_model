# Forward-running ETNP N2O box model
Requires Python 3.7 and Jupyter 4.4 or above

### DESCRIPTION:
    Model the time-dependent change in N2O isotopes from different substrate pools.
    Created on Mon October 28th, 2019
    
    Simple time-varying model to model the expected d15N2O-alpha, d15N2O-beta, and d18O-N2O from different
    reductive substrates (nitrate and nitrite). Utilizes the known isotope values of the substrate and a
    range of expected isotope effects to model the expected isotopes of N2O.
    
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
 	  
@author: Colette L. Kelly, Stanford University (clkelly@stanford.edu)

MIT License

Copyright (c) 2020 Colette LaMonica Kelly

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.