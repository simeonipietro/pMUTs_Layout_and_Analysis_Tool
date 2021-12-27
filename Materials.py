import numpy as np
from numpy.linalg import inv

# Material properties taken from COMSOL unless otherwise specified

# Piezoelectric Materials

AlN = {
    'Name' :'AlN',
    'Young' : 345e9, #[Pa]
    'Density' : 3300, #[kg/m^3]
    'e31_eff' : 1.08, # [C/m^2]
    'epsilon_r' : 9
}

def ScAlN(x=30):
 
    '''
    Sources: 
    - Caro et al., "Piezoelectric coefficients and spontaneous polarization of ScAlN"
    - Piazza and Bhugra, "Piezoelectric MEMS resonators", pp. 23 and pp. 156
    '''

    assert 0<=x<=50, 'The Sc percentage must be given as a number between 0 and 50!'

    x=x/100

    c11 = (410.2*(1-x)+295.3*x-210.3*x*(1-x))*1e9
    c12 = (142.4*(1-x)+198.6*x-61.9*x*(1-x))*1e9
    c13 = (110.1*(1-x)+135.5*x+78.9*x*(1-x))*1e9
    c33 = (385.0*(1-x)-23.8*x-101.4*x*(1-x))*1e9
    c44 = (122.9*(1-x)+169.5*x-137.3*x*(1-x))*1e9
    c66 = (c11-c12)/2

    #stiffness matrix
    c = np.array(
        [
            [c11, c12, c13, 0, 0, 0],
            [c12, c11, c13, 0, 0, 0],
            [c13, c13, c33, 0, 0, 0],
            [0, 0, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, 0], 
            [0, 0, 0, 0, 0, (c11-c12)/2]
        ]
    )

    #compliance matrix
    s = inv(c)

    e15 = -0.39*(1-x)-0.458*x+0.417*x*(1-x)
    e31 = -0.63*(1-x)-0.492*x-0.615*x*(1-x)
    e33 = 1.46*(1-x)+8.193*x-5.912*x*(1-x)

    young = 1/s[0,0]
    density = (4.6703e-05*x**2+6.3187e-04*x+3.512)*1e3
    e31_eff = e31 - e33*c13/c33
    epsilon_r = 40.293*x**2+3.7161*x+9.37

    ScAlN = {
        'Name' : str(x*100)+'% Scandium AlN',
        'Young' : young, #[Pa]
        'Density' : density, #[kg/m^3]
        'e31_eff' : e31_eff, # [C/m^2]
        'epsilon_r' : epsilon_r
        }

    return ScAlN


# Dielectric Materials

Silicon = {
    'Name' : 'Silicon',
    'Young' : 170e9, #[Pa]
    'Density' : 2329, #[kg/m^3]
    'epsilon_r' : 11.7
}

Silicon_Oxide = {
    'Name' : 'Silicon Oxide',
    'Young' : 70e9, #[Pa]
    'Density' : 2200, #[kg/m^3]
    'epsilon_r' : 4.2
}


# Metals

Platinum = {
    'Name' : 'Platinum',
    'Young' : 168e9, #[Pa]
    'Density' : 21450 #[kg/m^3]
}

Gold = {
    'Name' : 'Gold',
    'Young' : 70e9, #[Pa]
    'Density' : 19300 #[kg/m^3]
}

Tungsten = {
    'Name' : 'Tungsten',
    'Young' : 411e9, #[Pa]
    'Density' : 19350 #[kg/m^3]
}


Aluminum = {
    'Name' : 'Aluminum',
    'Young' : 70e9, #[Pa]
    'Density' : 2700 #[kg/m^3]
}

Molybdenum = {
    'Name' : 'Molybdenum',
    'Young' : 312e9, #[Pa]
    'Density' : 10200 #[kg/m^3]
}


# Fluids

Air = {
    'Name' : 'Air',
    'Speed' : 343, #[m/s]
    'Density' : 1.3 #[kg/m^3]
}

Water = {
    'Name' : 'Water',
    'Speed' : 1500, #[m/s]
    'Density' : 1000 #[kg/m^3]
}

