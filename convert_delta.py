def convert_delta(d=None, i=None):
    
    '''
    INPUTS:
    d = delta value
    i = isotope ("d15N" or "d18O") - to determine standard
    
    OUTPUT:
    R = 15R or 18R corresponding to input delta value
    '''
    
    R15std = 0.00367647 # air N2
    R18std = 0.00200517 # VSMOW
    
    if i == "d15N":
        R = ((d/1000)+1)*R15std
    elif i == "d18O":
        R = ((d/1000)+1)*R18std
    else:
        print('Please enter i as d15N or d18O')
        R = -999
    
    return R