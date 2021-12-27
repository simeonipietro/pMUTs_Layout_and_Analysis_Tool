#Imports
from re import T
import numpy as np
import os
import math
import cmath
import warnings
from phidl import Device
from phidl import quickplot as qp 
import phidl.geometry as pg
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams.update({'font.size': 16})
from sympy import symbols, diff, integrate
from scipy.special import jv, iv, struve
from scipy.optimize import fsolve
from pandas import Series, DataFrame
import schemdraw
import schemdraw.elements as elm

#Circular pMUT Geometry

def draw_center_electrode(
    pMUT_Radius,
    layer,
    center_electrode_radius_fraction,
    electrodes_gap
    ):

    Center_Electrode = Device()
    
    electrode_shape = pg.circle(
        radius=pMUT_Radius*center_electrode_radius_fraction-electrodes_gap/2,
        layer=layer
    )
    
    Center_Electrode.add_ref(electrode_shape)

    return Center_Electrode

def draw_outer_electrode_single_opening(
    pMUT_Radius, 
    layer,
    outer_electrode_radius_fraction,
    electrodes_gap,
    trace_width,
    orientation
    ):

    Outer_Electrode = Device()
    
    Opening_Width = trace_width+2*electrodes_gap
    
    electrode_ring = pg.ring(
                        radius=pMUT_Radius*(1-outer_electrode_radius_fraction)
                                +electrodes_gap/2
                                +(pMUT_Radius*outer_electrode_radius_fraction 
                                -electrodes_gap/2)/2,
                        width=pMUT_Radius*outer_electrode_radius_fraction-electrodes_gap/2,
                        layer=layer
                        )
    
    electrode_opening = pg.rectangle(size=(Opening_Width,pMUT_Radius))

    electrode_shape = pg.boolean(
                        A=electrode_ring, 
                        B=electrode_opening.movex(-Opening_Width/2),
                        operation='not', 
                        precision=1e-6,
                        num_divisions=[1, 1],
                        layer=layer
                        )

    # Adds geometry to the GDSII with phidl << operator
    Outer_Electrode << electrode_shape

    # rotate the pMUT as needed
    Outer_Electrode.rotate(orientation)

    return Outer_Electrode

def draw_outer_electrode_double_opening(
    pMUT_Radius, 
    layer,
    outer_electrode_radius_fraction,
    electrodes_gap,
    trace_width,
    orientation
    ):

    Outer_Electrode = Device()
    
    Opening_Width = trace_width+2*electrodes_gap
    
    electrode_ring = pg.ring(
                        radius=pMUT_Radius*(1-outer_electrode_radius_fraction)
                                +electrodes_gap/2
                                +(pMUT_Radius*outer_electrode_radius_fraction 
                                -electrodes_gap/2)/2,
                        width=pMUT_Radius*outer_electrode_radius_fraction-electrodes_gap/2,
                        layer=layer
                        )
    
    electrode_opening = pg.rectangle(size=(Opening_Width,2*pMUT_Radius))
    electrode_opening = electrode_opening.move([-Opening_Width/2,-pMUT_Radius])
    
    electrode_opening = electrode_opening.rotate(orientation)

    electrode_shape = pg.boolean(
                        A=electrode_ring, 
                        B=electrode_opening,
                        operation='not', 
                        precision=1e-6,
                        num_divisions=[1, 1],
                        layer=layer
                        )

    # Adds geometry to the GDSII with phidl << operator
    Outer_Electrode << electrode_shape

    return Outer_Electrode

def draw_etch_pit(
    pMUT_Radius, 
    layer
    ):

    Etch_Pit = Device()
    pit_shape = pg.circle(radius=pMUT_Radius, 
                          layer=layer
                          )
    Etch_Pit << pit_shape

    return Etch_Pit

def draw_center_electrode_single_trace(
    pMUT_Radius,
    layer,    
    trace_width, 
    orientation
    ):

    Connection = Device()

    connection_shape = pg.rectangle(
        size=(trace_width, pMUT_Radius), 
        layer=layer
        )

    connection_shape.ports = {}  # Removes default ports
    Connection << connection_shape  # Adds geometry to the GDSII 
    Connection.movex(-trace_width/2)

    Connection.rotate(orientation)

    return Connection

def draw_center_electrode_double_trace(
    pMUT_Radius,
    layer,    
    trace_width, 
    orientation
    ):

    Connection = Device()

    connection_shape = pg.rectangle(
        size=(trace_width, 2*pMUT_Radius), 
        layer=layer
        )

    connection_shape.ports = {}  # Removes default ports
    Connection << connection_shape  # Adds geometry to the GDSII 
    Connection.move([-trace_width/2,-pMUT_Radius])

    Connection.rotate(orientation)

    return Connection

def draw_rectangular_attachment(
    pMUT_Radius,
    layer,    
    trace_width, 
    orientation):

    Connection = Device()

    connection_shape1 = pg.rectangle(
        size=(trace_width, 0.1*pMUT_Radius), 
        layer=layer
        )

    connection_shape2 = pg.rectangle(
        size=(trace_width, 0.1*pMUT_Radius), 
        layer=layer
        )

    connection_shape1.ports = {}  # Removes default ports
    connection_shape2.ports = {}
    connection_shape1.move([-trace_width/2,-pMUT_Radius])
    connection_shape2.move([-trace_width/2,0.9*pMUT_Radius])
    Connection << connection_shape1  # Adds geometry to the GDSII 
    Connection << connection_shape2

    Connection.rotate(orientation)

    return Connection


# Circular pMUT Device Properties

def neutral_axis(thicknesses, young_moduli):
    '''
    Implementation taken from Piezoelectric MEMS book, pp. 154.
    thicknesses and young moduli are lists with the values for each layer
    '''
    zn = [0.5*thicknesses[0]]
    
    for i in range(1,len(thicknesses)):
        zn.append(sum(thicknesses[0:i]) + 0.5*thicknesses[i])
    
    z_neutral = np.sum(np.array(zn)*np.array(thicknesses)*np.array(young_moduli))/\
                np.sum(np.array(thicknesses)*np.array(young_moduli))

    return z_neutral

def flexural_rigidity(thicknesses, 
                      young_moduli,
                      poisson = 0.32
                      ):
    '''
    Implementation taken from Piezoelectric MEMS book, pp. 154.
    thicknesses and young moduli are lists with the values for each layer
    '''
    na = neutral_axis(thicknesses, young_moduli)
    hn = [thicknesses[0]]
    plate_moduli = np.array(young_moduli)/(1-poisson**2)

    for i in range(1,len(thicknesses)):
        hn.append(sum(thicknesses[0:i]) + thicknesses[i])
    
    hn_dash = np.array(hn)-na
    hn_diff = [hn_dash[0]**3+na**3]

    for i in range(1,len(thicknesses)):
        hn_diff.append(hn_dash[i]**3-hn_dash[i-1]**3)

    D = 1/3*np.sum(plate_moduli*np.array(hn_diff))

    return D

def pre_stressed_eigenfrequency(radius,
                                stresses,
                                thicknesses,
                                young_moduli,
                                densities,
                                poisson=0.32
                                ):

    '''
    Function implemented referencing the model presented in:

    'An Analytical Analysis of the Sensitivity of Circular
    Piezoelectric Micromachined Ultrasonic Transducers
    to Residual Stress'

    '''

    # Define plate tension, flexural rigidity, surface density, and neutral axis
    T = np.sum(stresses*thicknesses)
    D = flexural_rigidity(thicknesses, young_moduli, poisson)
    rho = np.sum(densities*thicknesses)

    # Solve numerically for alpha_k and beta_k
    buckled = False
    phi = T/(14.68*D/radius**2)

    if phi<-1:
        buckled = True

    def equations(p):
        alphak, betak = p
        
        equation1 = (
                    alphak*(jv(1, alphak)/jv(0, alphak)) 
                    +betak*(iv(1, betak)/iv(0, betak))
                    )
        equation2 = betak**2 - alphak**2 - 14.68*phi
        
        return (equation1, equation2)

    # I define the starting point for the solver based on the mode
    # we are targeting. Depending on the stress levels  and geometry
    # you are using you might get a warning that the solver had trouble 
    # converging and you might need to play with different starting points
    # See this link for reference values I used to select the starting points:
    # https://asa.scitation.org/doi/pdf/10.1121/1.1928110

    a_init = 3.2
    b_init = 3.2
    alpha_k, beta_k =  fsolve(equations, (a_init, b_init))    

    # I re-express the frequency-dependent term at the denominator of equation (25) to get
    # a more compact expression for the resonance frequency
    A = alpha_k**2/radius**2
    B = T/2/D

    return np.sqrt(D/rho*A*(2*B+A))/2/np.pi, buckled

def equivalent_mass(radius, thicknesses, densities):
    mu = np.sum(np.array(thicknesses)*np.array(densities))
    return np.pi*radius**2*mu/5

def equivalent_stiffness(radius,
                         stresses,
                         thicknesses,
                         young_moduli,
                         densities,
                         poisson=0.32
                         ):

    f0, buckled = pre_stressed_eigenfrequency(radius,
                                     stresses,
                                     thicknesses,
                                     young_moduli,
                                     densities,
                                     poisson
                                     )
    Meq = equivalent_mass(radius, thicknesses, densities)
    
    return (2*np.pi*f0)**2*Meq

def coupling_coefficient(e31_effective, 
                         thicknesses, 
                         young_moduli,
                         central_electrode_radius_fraction, 
                         piezo_position, 
                         electrode_type 
                         ):

    '''
    Taken from Piezoelectric MEMS Resonators, pp. 168.
    I don't understand where the 0.5 factor comes from in their expression, and
    if included it does not match their own test shown in Fig. 6.2. So I removed it.
    '''
    
    na = neutral_axis(thicknesses, young_moduli)
    zn = [0.5*thicknesses[0]]
    for i in range(1,len(thicknesses)):
        zn.append(sum(thicknesses[0:i]) + 0.5*thicknesses[i])

    zp = abs(zn[piezo_position]-na)

    x = symbols('x')
    modeshape = (1-x**2)**2
    expression = x*diff(modeshape, x, 2) + diff(modeshape, x)
    I_piezo = 2*np.pi*(\
                electrode_type[0]*integrate(expression, (x,0, central_electrode_radius_fraction))\
               -electrode_type[1]*integrate(expression, (x, central_electrode_radius_fraction, 1))
               )

    return float(e31_effective*zp*abs(I_piezo)) 

def electrical_capacitance(radius, 
                           piezo_thickness,
                           central_electrode_radius_fraction, 
                           electrode_type, 
                           electrodes_gap,
                           relative_dielectric_constant
                           ):

    Cel = 8.854e-12*relative_dielectric_constant/piezo_thickness*np.pi*\
        (\
        electrode_type[0]*(radius*central_electrode_radius_fraction-electrodes_gap/2)**2
        + electrode_type[1]*(radius**2-(radius*central_electrode_radius_fraction+electrodes_gap/2)**2)
        )

    return Cel

def effective_area(radius):
    return np.pi*radius**2/3

def radiation_impedance(medium_density, medium_speed, pMUTradius, frequency):
    
    'pMUT radiation impedance in the mechanical domain'

    a = np.sqrt(effective_area(pMUTradius)/np.pi)
    ka = a*2*np.pi*frequency/medium_speed

    return medium_density*medium_speed*np.pi*a**2 * (1 - jv(1,2*ka)/ka + 1j*struve(1, 2*ka)/ka)#medium_density*medium_speed*np.pi*a**2 * (2 - jv(1,4*ka)/ka + 1j*struve(1, 4*ka)/ka)

def pMUT_DF(theta, medium_speed, frequency, pMUTradius):

    ''' 
    This function calculates the directivity function of
    the circular pMUT based on its dimensions, frequency,
    and properties of the medium.

    Theta is the angle between the xy plane and the z-axis
    as the convention used in physics. For any point on the 
    positive z-axis, theta = 0.
    '''
    a = np.sqrt(effective_area(pMUTradius)/np.pi)
    ka = a*2*np.pi*frequency/medium_speed

    return jv(1, ka*np.sin(theta))/(ka*np.sin(theta))

def Array_DF(thetas, phi, coord, phases, medium_speed, pMUTradius, frequency, FF):

    '''
    Array_DF outputs an array of points between 0 and 1, representing the directivity
    strength for different angles theta and at a fixed azimuthal angle phi. 
    The values of thetas and phi are given as an input.
    The desired azimuthal angle phi is also given as an input to select the slice
    of radiation pattern that is observed.
    '''
    
    '''
    coord represents the positions of the pMUTs on the xy plane.
    coord must be an array of tuples.
    Each tuple has two slots with coordinates (x, y).
    coord should have the origin at the centerof the array.
    '''

    '''
    FF Defines the far-field radius from the coordinates (0,0)
    '''

    # Defines wave number in the medium
    k  = 2*np.pi*frequency/medium_speed

    # Defines distance ri between each pMUT and the far-field evaluation points.
    # Evaluating ri is necessary to account of the phase differences between
    # the contributions of the individual pMUTs.

    # ri is a 2D array with dimensions (# evaluation points, # pMUTs), which
    # correspond to len(thetas) and len(coord), respectively.
    
    ri = np.zeros((len(thetas), len(coord)))
    iso_sum = np.zeros(len(thetas))

    evalX = FF*np.cos(phi)*np.sin(thetas)
    evalY = FF*np.sin(phi)*np.sin(thetas)
    evalZ = FF*np.cos(thetas)
    
    for i in range(len(thetas)):

        ri[i] = np.sqrt((evalX[i]-coord[:,0])**2 + (evalY[i]-coord[:,1])**2 + evalZ[i]**2)

        # Sum all the isotropic sources contributions

        iso = 0

        for j in range(len(coord)):
            r = ri[i,j]
            p = phases[j]
            iso += cmath.exp(-1j*(k*r-p))/r

        iso_sum[i] = abs(iso)

    DF = abs(pMUT_DF(thetas, medium_speed, frequency, pMUTradius) * iso_sum)
    # Use acoustic product theorem to combine the radiation pattern
    # of the single pMUT with the isotropic sources       
  
    return DF/np.amax(DF)


# Square Microphone - Triangular Cantilevers Geometry

def triangular_electrodes(side, offset, layer):

    electrode = pg.rectangle(size = (side,side),
                             layer = layer)
    electrode.move([-side/2,-side/2])

    cut = pg.cross(length = 2*side, 
                   width = 2*offset, 
                   layer = layer)
    cut.rotate(45)

    electrode = pg.boolean(A = electrode, 
                           B = cut, operation = 'not', 
                           precision = 1e-6,
                           num_divisions = [1,1], 
                           layer = layer)

    return electrode

def triangular_mic_slots(side, slot_width, layer):

    slots = pg.cross(length = np.sqrt(2)*side, 
                   width = slot_width, 
                   layer = layer)

    slots.rotate(45)

    return slots

def triangular_mic_back_etch(side, layer):
    
    etch = pg.rectangle(size=(side,side),layer=layer)
    
    return etch.move([-side/2, -side/2])