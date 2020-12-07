import pandas as pd # dataframe management
import numpy as np # numerical functions
import matplotlib.pyplot as plt # plotting
import statsmodels.api as sm # regression analysis

from convert_delta import convert_delta # custom - convert deltas to isotope ratios

def exponential(time_vector=None,initial=None,decayrate=None,minimum=None):
    '''
    INPUTS:
    time_vector = vector of times (positive values)
    initial = value at time zero, minus the minimum
    decayrate = how quickly it goes to the minimum value.
    minimum = minimum value
    ___
    NOTE: smaller decay rate = closer to linear function over short timescales.
    
    OUTPUT:
    output = vector w/ same shape as time_vector
    '''
    
    exponential = initial*2.**(-decayrate*time_vector)+minimum
    
    return exponential

def logfxn(time_vector=None, scaling=None, maximum=None, minimum=None):
    '''
    INPUTS:
    time_vector = vector of times (positive values)
    scaling = sharpness of the bend in the plot (smaller=sharper)
    maximum = asymptote value (can be positive or negative)
    minimum = minimum value
    ___
    NOTE: smaller decay rate = closer to linear function over short timescales.
    
    OUTPUT:
    logfxn = vector w/ same shape as time_vector
    '''
    
    logfxn = (maximum-minimum)*(time_vector/(time_vector+scaling))+minimum
    
    return logfxn

class initialize(object):
    '''
    INPUTS:
    alpha_values = isotope effects in per mil (list):
                        [NOx -> N2Oa, N2Oa -> N2,NOx -> N2Ob, N2Ob -> N2,
                        NO3- -> NO2- branching isotope effect,
                        NO2- -> N2O, branching isotope effect,
                        NO2- -> N2O, kinetic isotope effect,
                        N2O -> N2, kinetic isotope effect,
                        equilibration with seawater, exhange]    
    
    n2o_cons_constants = [scaling, maximum] for logfxn on N2O consumption
    n2o_prod_constants = 
    
    substrate = "NO3-" or "NO2-"
    
    OUTPUT:
    output = object with properties alpha_values,
    n2o_cons_constants, n2o_prod_constants, substrate,
    init_from, initials
    
    '''
    
    def __init__(self, alpha_values=None, n2o_cons_constants=None,
                n2o_prod_constants=None, substrate=None,
                 init_from=None, initials=None):

        [self.N2O_init, self.d15N2Oa_init,
         self.d15N2Ob_init, self.d18N2O_init,
         self.NO2_init, self.d15NO2_init,
         self.d18NO2_init, self.NO3_init,
         self.d15NO3_init, self.d18NO3_init] = self.get_initial_isotopes(init_from, initials)
        
        [self.alpha15noxA, self.alpha15N2OCONSA,
        self.alpha15noxB, self.alpha15N2OCONSB,
        self.alpha18NO3TON2O_b, self.alpha18NO2TON2O_b,
        self.alpha18NO2TON2O, self.alpha18N2OCONS,
        self.alpha18h2o, self.alphaexch] = self.define_isotope_effects(alpha_values)
        
        [self.n2o_cons_scaling,
         self.n2o_cons_max,
         self.n2o_cons_min] = self.get_kN2OCONS(n2o_cons_constants)
        
        [self.NO2TON2O_scaling, self.NO2TON2O_max,
         self.NO2TON2O_min, self.NO3TON2O_scaling,
         self.NO3TON2O_max, self.NO3TON2O_min] = self.get_kN2OPROD(n2o_prod_constants,
                                           substrate)
    
    def get_initial_isotopes(self, init_from, initials):
        
        if init_from=='ODZ' or init_from==None:

            N2O_init = 92.4 # nmol/L; mean at stn. PS1 in sigma-theta range 26.0-27.0
            d15N2Oa_init = 12.3 # per mil
            d15N2Ob_init = -1.6 # per mil
            d18N2O_init = 51.3 # per mil

            NO2_init = 2.0 # umol/L; mean at stn. PS3 in sigma-theta range 26.0-27.0
            d15NO2_init = -22.0 # per mil
            d18NO2_init = 16.6 # per mil

            NO3_init = 26.5 # umol/L; mean at stn. PS3 in sigma-theta range 26.0-27.0
            d15NO3_init = 16.3 # per mil
            d18NO3_init = 15.0 # per mil
        
        elif init_from=='ATM':
            
            N2O_init = 6.3 # nmol/L; saturation value based on mean surface T&S
            d15N2Oa_init = 15.7 # per mil; atmospheric N2O (Mohn et al., 2014)
            d15N2Ob_init = -3.3 # per mil; atmospheric N2O (Mohn et al., 2014)
            d18N2O_init = 44.3 # per mil; atmospheric N2O (Mohn et al., 2014)

            NO2_init = 0.5 # umol/L; mean at [NO2-]>=0.1uM & sigma-theta<=25.0
            d15NO2_init = -1.0 # per mil; mean at [NO2-]>=0.1uM & sigma-theta<=25.0
            d18NO2_init = 12. # per mil; mean at [NO2-]>=0.1uM & sigma-theta<=25.0

            NO3_init = 10.7 # umol/L; mean at [NO2-]<0.1uM & sigma-theta<=25.0
            d15NO3_init = 7.0 
            d18NO3_init = 19.0 
            
        elif init_from=='custom':
            
            N2O_init = initials['[N2O]_nM']
            d15N2Oa_init = initials.d15N2Oa
            d15N2Ob_init = initials.d15N2Ob
            d18N2O_init = initials.d18N2O

            NO2_init = initials['[NO2-]_uM']
            d15NO2_init = initials['d15NO2-']
            d18NO2_init = initials['d18NO2-']

            NO3_init = initials['[NO3-]_uM']
            d15NO3_init = initials['d15NO3-']
            d18NO3_init = initials['d18NO3-']
        
        return [N2O_init,d15N2Oa_init,d15N2Ob_init,d18N2O_init,
                NO2_init,d15NO2_init,d18NO2_init,
                NO3_init,d15NO3_init,d18NO3_init]
    
    def define_isotope_effects(self, alpha_values):
        
        if alpha_values==None:
            # Isotope effects, trans into fractionation factors
            alpha15noxA = 0./1000.+1 # NOx -> N2Oa - need 40 per mil SP for denit to get observed profiles
            alpha15N2OCONSA = 11.8/1000.+1  # N2Oa -> N2
            alpha15noxB = 0./1000.+1 # NOx -> N2Ob
            alpha15N2OCONSB = 0./1000.+1 #1.64/1000.+1  # N2Ob -> N2
            alpha18NO3TON2O_b = 24./1000.+1 # NO3- -> NO2-, branching isotope effect
            alpha18NO2TON2O_b = 12./1000.+1 # NO2- -> N2O, branching isotope effect
            alpha18NO2TON2O = 0./1000.+1 # NO2- -> N2O, kinetic isotope effect
            alpha18N2OCONS = 20.2/1000.+1 # N2O -> N2, kinetic isotope effect
            alpha18h2o = 10/1000.+1  # equilibration with seawater
            alphaexch = 1.0153
        else:
            alpha15noxA = alpha_values[0]/1000.+1
            alpha15N2OCONSA = alpha_values[1]/1000.+1
            alpha15noxB = alpha_values[2]/1000.+1
            alpha15N2OCONSB = alpha_values[3]/1000.+1
            alpha18NO3TON2O_b = alpha_values[4]/1000.+1
            alpha18NO2TON2O_b = alpha_values[5]/1000.+1
            alpha18NO2TON2O = alpha_values[6]/1000.+1
            alpha18N2OCONS = alpha_values[7]/1000.+1
            alpha18h2o = alpha_values[8]/1000.+1
            alphaexch = alpha_values[9]/1000.+1
            
        return [alpha15noxA, alpha15N2OCONSA,
        alpha15noxB, alpha15N2OCONSB,
        alpha18NO3TON2O_b, alpha18NO2TON2O_b,
        alpha18NO2TON2O, alpha18N2OCONS,
        alpha18h2o, alphaexch]
    
    def get_kN2OCONS(self, n2o_cons_constants):
    
        if n2o_cons_constants == None:
            n2o_cons_scaling = 1000.
            n2o_cons_max = 0.21
            n2o_cons_min = 0.0
        else:
            n2o_cons_scaling = n2o_cons_constants[0]
            n2o_cons_max = n2o_cons_constants[1]
            n2o_cons_min = n2o_cons_constants[2]
                
        return [n2o_cons_scaling, n2o_cons_max, n2o_cons_min]
    
    def get_kN2OPROD(self, n2o_prod_constants, substrate):
        
        if n2o_prod_constants == None:
            if substrate=='NO2-':
                # NO2 -> N2O ("NO3TON2O")
                NO2TON2O_scaling = 0.
                NO2TON2O_max = .0000525 
                NO2TON2O_min = 0.0
                
                # NO3 -> N2O ("NO3TON2O")
                NO3TON2O_scaling = 0.0
                NO3TON2O_max = 0.0
                NO3TON2O_min = 0.0
                
            elif substrate=='NO3-' or substrate==None:
                # NO2 -> N2O ("NO3TON2O")
                NO2TON2O_scaling = 0.0
                NO2TON2O_max = 0.0
                NO2TON2O_min = 0.0
                
                # NO3 -> N2O ("NO3TON2O")
                NO3TON2O_scaling = 0.
                NO3TON2O_max = .000003962
                NO3TON2O_min = 0.0
     
            else:
                print('Please input valid choice of substrate')
        else:
            NO2TON2O_scaling = n2o_prod_constants[0]
            NO2TON2O_max = n2o_prod_constants[1]
            NO2TON2O_min = n2o_prod_constants[2]
            
            NO3TON2O_scaling = n2o_prod_constants[3]
            NO3TON2O_max = n2o_prod_constants[4]
            NO3TON2O_min = n2o_prod_constants[5]
                
        return [NO2TON2O_scaling, NO2TON2O_max, NO2TON2O_min,
               NO3TON2O_scaling, NO3TON2O_max, NO3TON2O_min]
    
    def __repr__(self):
        
        return "Model parameters initialized"

def run_model(params):
    
    '''
    OUTPUT:
    output = pandas dataframe with columns:
                '[NO3-]_uM', 'd15NO3-', 'd18NO3-',
                '[NO2-]_uM', 'd15NO2-', 'd18NO2-',
                '[N2O]_nM', 'd15N2Oa', 'd15N2Ob', 'd18N2O'
     '''           
     
    ### INITIALIZE MODEL PARAMETERS ###

    # time step (d)
    dt = 0.2

    # number of timesteps (y)
    # increasing n(timesteps) by a factor of 10 decreases rate constants by the same factor
    T = 1000
    times = np.array(list(range(1,T+1))) # vector of timesteps
    
    ### TIME-VARYING N2O CONSUMPTION ###
    kN2OCONS = logfxn(time_vector=np.array(list(range(1,T+1))),scaling=params.n2o_cons_scaling,
                      maximum=params.n2o_cons_max, minimum=params.n2o_cons_min) # N2O -> N2 ("N2OCONS")

    ### TIME-VARYING N2O PRODUCTION ###
    kNO2TON2O = logfxn(time_vector=np.array(list(range(1,T+1))),
                       scaling=params.NO2TON2O_scaling,
                       maximum=params.NO2TON2O_max,
                       minimum=params.NO2TON2O_min) # N2O -> N2 ("N2OCONS")
    kNO3TON2O = logfxn(time_vector=np.array(list(range(1,T+1))),
                       scaling=params.NO3TON2O_scaling,
                       maximum=params.NO3TON2O_max,
                       minimum=params.NO3TON2O_min)

    ### SUBSTRATE ISOTOPES ###
    d15NO2  = np.ones(shape=(T,))
    d15NO2[:,] = params.d15NO2_init
    R15NO2 = convert_delta(d=d15NO2,i='d15N')
    no2_14  = np.ones(shape=(T,))
    no2_14[:,] = params.NO2_init
    no2_15 = no2_14*R15NO2

    d18NO2  = np.ones(shape=(T,))
    d18NO2[:,] = params.d18NO2_init
    R18NO2 = convert_delta(d=d18NO2,i='d18O')
    no2_16 = 2*no2_14 # 2 O atoms per every N
    no2_18 = no2_16*R18NO2

    d15NO3 = np.ones(shape=(T,))
    d15NO3[:,] = params.d15NO3_init
    R15NO3 = convert_delta(d=d15NO3,i='d15N')
    no3_14  = np.ones(shape=(T,))
    no3_14[:,] = params.NO3_init
    no3_15 = no3_14*R15NO3

    d18NO3 = np.ones(shape=(T,))
    d18NO3[:,] = params.d18NO3_init
    R18NO3 = convert_delta(d=d18NO3,i='d18O')
    no3_16 = 3*no3_14 # 3 O atoms per every N
    no3_18 = no3_16*R18NO3

    ### ISOTOPE CONSTANTS ###
    R15std = 0.00367647 # air N2
    R18std = 0.00200517 # VSMOW
    # exchange factor for 18O in no2
    E = 0.09 # exchange
    kexch = E*np.ones(shape=(T,T)) # get TxT array of exchange
    d18h2o = 0.
    R18h2o = convert_delta(d=d18h2o, i="d18O")


    ### INITIALIZE STATE VARIABLES ###

    # initialize with PS1 boundary conditions
    [N2O_i, d15n2oA, d15n2oB, d18n2o] = [params.N2O_init, params.d15N2Oa_init,
                                         params.d15N2Ob_init, params.d18N2O_init]

    R15n2oA_i = convert_delta(d=d15n2oA,i="d15N")
    R15n2oB_i = convert_delta(d=d15n2oB, i="d15N")
    R18n2o_i = convert_delta(d=d18n2o, i="d18O")

    # 5 state variables
    n2o_14_i = N2O_i/1000. # umol/L
    n2o_15A_i = R15n2oA_i*n2o_14_i
    n2o_15B_i = R15n2oB_i*n2o_14_i
    n2o_16_i = 0.5*n2o_14_i # 2 nitrogen atoms per oxygen!
    n2o_18_i = R18n2o_i*n2o_16_i

    # initialize arrays of state variables
    n2o_14 = np.zeros(shape = (T,1))
    n2o_15A = np.zeros(shape = (T,1))
    n2o_15B = np.zeros(shape = (T,1))
    n2o_16 = np.zeros(shape = (T,1))
    n2o_18 = np.zeros(shape = (T,1))

    # initial values of state variables
    n2o_14[0,:] = n2o_14_i
    n2o_15A[0,:] = n2o_15A_i
    n2o_15B[0,:] = n2o_15B_i
    n2o_16[0,:] = n2o_16_i
    n2o_18[0,:] = n2o_18_i

    ### TIME STEPPING ###
    for iT in range(T-1):

        # update values
        n2o_14[iT+1,:] = n2o_14[iT,:] + dt*(kNO3TON2O[iT]*no3_14[iT] # NO3- -> N2O
                                            +kNO2TON2O[iT]*no2_14[iT] # NO2- -> N2O
                                            -kN2OCONS[iT]*n2o_14[iT,:]) # N2O- -> N2

        n2o_15A[iT+1,:] = n2o_15A[iT,:] + dt*(kNO3TON2O[iT]/params.alpha15noxA*no3_15[iT]
                                            +kNO2TON2O[iT]/params.alpha15noxA*no2_15[iT]
                                            -kN2OCONS[iT]/params.alpha15N2OCONSA*n2o_15A[iT,:])

        n2o_15B[iT+1,:] = n2o_15B[iT,:] + dt*(kNO3TON2O[iT]/params.alpha15noxB*no3_15[iT]
                                            +kNO2TON2O[iT]/params.alpha15noxB*no2_15[iT]
                                            -kN2OCONS[iT]/params.alpha15N2OCONSB*n2o_15B[iT,:])

        n2o_16[iT+1,:] = n2o_16[iT,:] + dt*(kNO3TON2O[iT]*0.5*no3_14[iT]
                                            +kNO2TON2O[iT]*0.5*no2_14[iT]
                                            -kN2OCONS[iT]*n2o_16[iT,:])

        n2o_18[iT+1,:] = n2o_18[iT,:] + dt*(kNO3TON2O[iT]*params.alpha18NO3TON2O_b*params.alpha18NO2TON2O_b/params.alpha18NO2TON2O*0.5/3*no3_18[iT]
                                            +kNO2TON2O[iT]*params.alpha18NO2TON2O_b/params.alpha18NO2TON2O*0.5/2*no2_18[iT]
                                            -kN2OCONS[iT]/params.alpha18N2OCONS*n2o_18[iT,:])

    '''
    -kexch[iT,iT]*no3_18[iT]
    +kexch[iT,iT]*no3_16[iT]*alphaexch*R18h2o
    '''

    ### CALCULATE OUTPUT ###
    no3_concentration = pd.DataFrame(no3_14, columns={'[NO3-]_uM'})
    d15no3 = pd.DataFrame(d15NO3, columns={'d15NO3-'})
    d18no3 = pd.DataFrame(d18NO3, columns={'d18NO3-'})

    no2_concentration = pd.DataFrame(no2_14, columns={'[NO2-]_uM'})
    d15no2 = pd.DataFrame(d15NO2, columns={'d15NO2-'})
    d18no2 = pd.DataFrame(d18NO2, columns={'d18NO2-'})

    n2o_concentration = pd.DataFrame(n2o_14*1000, columns={'[N2O]_nM'})
    d15N2O_A = pd.DataFrame((((n2o_15A/(n2o_14+n2o_15B))/R15std)-1)*1000, columns={'d15N2Oa'})
    d15N2O_B = pd.DataFrame((((n2o_15B/(n2o_14+n2o_15A))/R15std)-1)*1000, columns={'d15N2Ob'})
    d18N2O = pd.DataFrame((((n2o_18/n2o_16)/R18std)-1)*1000, columns={'d18N2O'})

    output_array = no3_concentration.join([d15no3,d18no3,no2_concentration,d15no2,d18no2,
                                     n2o_concentration,d15N2O_A,d15N2O_B,d18N2O])
    
    return output_array