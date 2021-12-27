from Materials import *
from pMUT_Functions import *

# Basic pMUT Element Class 
class Basic_pMUT:

    '''
    Variables description:
    pMUT_Radius: Radius of the pMUT
    *Layer: Layer number in the layout
    Center_Electrode_Fraction: Fraction of the pMUT radius occupied by the center electrode
    Electrodes_Gap: Gap between outer and center electrode
    Angular_Coverage: Fraction of the outer electrode inner circumference degrees 
    Orientation: Angle that controls the orientation of the pMUT
    '''

    def __init__(
        self, 
        pMUT_Radius=40e-6, 
        BottomMetalLayer=0,
        ViasLayer=1,
        TopMetalLayer=2,
        BackEtchLayer=3, 
        Center_Electrode_Fraction=0.55,
        Electrodes_Gap=2.5e-6,
        Connection_Width=5e-6,
        Orientation=0,
        Elastic_Material = Silicon_Oxide,
        Elastic_Thickness = 1e-6,
        Elastic_Stress = 0e6,
        Bottom_Electrode_Material = Molybdenum,
        Bottom_Electrode_Thickness = 0e-9,
        Bottom_Electrode_Stress = 0e6,
        Piezo_Material = ScAlN(0),
        Piezo_Thickness = 0.5e-6,
        Piezo_Stress = 0e6,
        Top_Electrode_Material = Aluminum,
        Top_Electrode_Thickness = 0e-9,
        Top_Electrode_Stress = 0e6,
        Medium_Material = Water
        ):
        
        assert isinstance(BottomMetalLayer,int), 'Bottom metal layer must be int!'
        assert isinstance(TopMetalLayer,int), 'Top metal layer must be int!'
        assert isinstance(BackEtchLayer,int), 'Backetch layer must be int!'
        assert pMUT_Radius>0, 'The pMUT radius must be positive!'
        assert 0.1<Center_Electrode_Fraction<0.9, 'Center_Electrode_Fraction must be between 0.1 and 0.9'
        assert 0<=Electrodes_Gap<pMUT_Radius, 'The Electrodes_Gap value is invalid'
        assert Connection_Width<2*pMUT_Radius , 'Connection_Width must be smaller than the center electrode diameter'
        assert isinstance(Orientation, int) or\
               isinstance(Orientation, float), 'Orientation must be a real number'
        assert isinstance(Elastic_Material, dict), 'Elastic_Material must be a dictionary with material properties'
        assert isinstance(Piezo_Material, dict), 'Piezo_Material must be a dictionary with material properties'
        assert isinstance(Bottom_Electrode_Material, dict), 'Bottom_Electrode_Material must be a dictionary with material properties'
        assert isinstance(Top_Electrode_Material, dict), 'Top_Electrode_Material must be a dictionary with material properties'
        assert isinstance(Medium_Material, dict), 'Medium_Material must be a dictionary with material properties'      
        
        self.pMUT_Radius = pMUT_Radius
        self.bottom_metal_layer = BottomMetalLayer
        self.vias_layer = ViasLayer
        self.top_metal_layer = TopMetalLayer
        self.back_etch_layer = BackEtchLayer
        self.electrodes_gap = Electrodes_Gap
        self.center_electrode_radius_fraction = Center_Electrode_Fraction
        self.connection_width = Connection_Width
        self.orientation = Orientation
        self.center_outer_electrode_presence = None
        self.e31_eff = Piezo_Material['e31_eff']
        self.relative_dielectric_constant = Piezo_Material['epsilon_r']
        self.layout = Device()

        self.thicknesses = np.array([Elastic_Thickness, 
                            Bottom_Electrode_Thickness, 
                            Piezo_Thickness, 
                            Top_Electrode_Thickness
                            ])

        self.young_moduli = np.array([Elastic_Material['Young'], 
                             Bottom_Electrode_Material['Young'], 
                             Piezo_Material['Young'], 
                             Top_Electrode_Material['Young']
                             ])

        self.densities = np.array([Elastic_Material['Density'], 
                          Bottom_Electrode_Material['Density'], 
                          Piezo_Material['Density'], 
                          Top_Electrode_Material['Density']
                          ])

        self.stresses = np.array([Elastic_Stress,
                         Bottom_Electrode_Stress,
                         Piezo_Stress,
                         Top_Electrode_Stress
                         ])

        self.mat_names = np.array([Elastic_Material['Name'], 
                          Bottom_Electrode_Material['Name'], 
                          Piezo_Material['Name'], 
                          Top_Electrode_Material['Name']
                          ])
        
        self.medium_density = Medium_Material['Density']
        self.medium_speed = Medium_Material['Speed']

    def __repr__(self):

        # Geometry and layers properties
        properties = [
                      'pMUT radius',
                      'Elastic layer material',
                      'Elastic layer thickness',
                      'Elastic layer stress',
                      'Bottom electrode material',
                      'Bottom electrode thickness',
                      'Bottom electrode stress',
                      'Piezo layer material',
                      'Piezo layer thickness',
                      'Piezo layer stress',
                      'Top electrode material',
                      'Top electrode thickness',
                      'Top electrode stress'
                      ]

        prop_values = [
                      '{:.0f}'.format(self.pMUT_Radius*1e6),
                      self.mat_names[0],
                      '{:.1f}'.format(self.thicknesses[0]*1e6),
                      '{:.0f}'.format(self.stresses[0]/1e6),
                      self.mat_names[1],
                      '{:.0f}'.format(self.thicknesses[1]*1e9),
                      '{:.0f}'.format(self.stresses[1]/1e6),
                      self.mat_names[2],
                      '{:.1f}'.format(self.thicknesses[2]*1e6),
                      '{:.0f}'.format(self.stresses[2]/1e6),
                      self.mat_names[3],
                      '{:.0f}'.format(self.thicknesses[3]*1e9),
                      '{:.0f}'.format(self.stresses[3]/1e6),
                      ]

        prop_unit = [
                     'µm',
                     '-',
                     'µm',
                     'MPa',
                     '-',
                     'nm',
                     'MPa', 
                     '-',
                     'µm',
                     'MPa',
                     '-',
                     'nm',
                     'MPa'                     
                     ]

        df_prop = DataFrame({
                         'VALUE': prop_values, 
                         'UNIT': prop_unit
                         },
                         index=properties
                         )

        # Equivalent parameters
        parameters = [
                      'f_eigen',
                      'f_resonance',
                      'Keq',
                      'Meq',
                      'Cel',
                      'η',
                      'Aeff',
                      'R_rad',
                      'M_rad'
                      ]

        values = [
                '{:.3f}'.format(self.eigenfrequency/1e6),
                '{:.3f}'.format(self.resonance_frequency/1e6),
                '{:.3f}'.format(self.equivalent_stiffness),
                '{:.3f}e-12'.format(self.equivalent_mass*1e12),
                '{:.3f}'.format(self.electrical_capacitance*1e12),
                '{:.3f}e-7'.format(self.eta*1e7),
                '{:.3f}e-9'.format(self.effective_area*1e9),
                '{:.3f}e-6'.format(self.radiation_impedance.real*1e6),
                '{:.3f}e-12'.format(self.radiation_impedance.imag*1e12/
                                   2/np.pi/self.resonance_frequency)
                ]

        units = [
                '[MHz]',
                '[MHz]',
                '[N/m]',
                '[kg]',
                '[pF]',
                '[N/V]',
                '[m^2]',
                '[kg s]',
                '[kg]'
                ]

        description = [
                      'Resonance frequency in vacuum',
                      'Resonance frequency in the medium',
                      'Equivalent stiffness of the pMUT',
                      'Equivalent mass of the pMUT',
                      'Electrodes capacitance',
                      'Coupling coefficient',
                      'pMUT effective area',
                      'Radiation resistance',
                      'Radiation mass at resonance'
                      ]

        df_eq = DataFrame({
                         'VALUE': values, 
                         'UNIT': units,
                         'DESCRIPTION': description
                         },
                         index=parameters
                         )

        return df_prop.to_string() + '\n\n' + df_eq.to_string()

    @property
    def neutral_axis(self):
        return neutral_axis(self.thicknesses, self.young_moduli)

    @property    
    def flexural_rigidity(self):
        return flexural_rigidity(self.thicknesses, self.young_moduli)

    @property
    def eigenfrequency(self):
        
        f,buckled = pre_stressed_eigenfrequency(self.pMUT_Radius,
                                             self.stresses,
                                             self.thicknesses,
                                             self.young_moduli,
                                             self.densities
                                             )
        return f

    @property    
    def equivalent_stiffness(self):
        
        k = equivalent_stiffness(self.pMUT_Radius,
                                 self.stresses,
                                 self.thicknesses,
                                 self.young_moduli,
                                 self.densities
                                 )
        
        return k

    @property    
    def equivalent_mass(self):
        
        m = equivalent_mass(self.pMUT_Radius,
                            self.thicknesses,
                            self.densities
                            )
        return m

    @property
    def eta(self):
        
        eta = coupling_coefficient(self.e31_eff,
                                   self.thicknesses,
                                   self.young_moduli,
                                   self.center_electrode_radius_fraction,
                                   2,
                                   self.center_outer_electrode_presence
                                   )
        return abs(eta)

    @property 
    def electrical_capacitance(self):
        
        Cel = electrical_capacitance(self.pMUT_Radius, 
                                     self.thicknesses[2],
                                     self.center_electrode_radius_fraction, 
                                     self.center_outer_electrode_presence, 
                                     self.electrodes_gap,
                                     self.relative_dielectric_constant
                                     )
        return Cel

    @property
    def effective_area(self):

        return np.pi*self.pMUT_Radius**2/3

    @property
    def resonance_frequency(self):

        # Resonance frequency of the device accounting for the medium loading

        f0 = 0
        f1 = self.eigenfrequency
        err = 0.1

        keq = self.equivalent_stiffness
        meq = self.equivalent_mass

        while err>0.001:

                m_medium = radiation_impedance(self.medium_density, 
                                               self.medium_speed, 
                                               self.pMUT_Radius, 
                                               f1
                                               ).imag/2/np.pi/f1

                f0 = f1
                f1 = 1/2/np.pi*np.sqrt(keq/(meq+m_medium))  
                err = abs(f1-f0)/f0 

        correction_coefficient = 1.12**2
        m_medium = m_medium*correction_coefficient
        f1 = 1/2/np.pi*np.sqrt(keq/(meq+m_medium)) 

        return f1
    
    @property
    def radiation_impedance(self):

        return radiation_impedance(self.medium_density, 
                                self.medium_speed, 
                                self.pMUT_Radius, 
                                self.resonance_frequency)

    def eq_circuit(self):

        d = schemdraw.Drawing()

        # Left side of first transformer
        HD = d.add(elm.Dot(open=True).label('V', 'left'))
        d += elm.Line().right().length(1)
        d += elm.Capacitor().down().label('Cel', loc='top')
        d += elm.Ground()
        d.push()
        d += elm.Line().left().length(1)
        LD = d.add(elm.Dot(open=True).left())
        d += elm.Line().right().length(2)
        d += elm.Line().up().length(0.6)
        d += elm.Line().right().length(0.3)
        T = d.add(elm.Transformer().right().label('1:η', 'top'))
        d += elm.Line().left().length(0.3).at(T.p1)
        d += elm.Line().up().toy(HD.start)
        d += elm.Line().left().length(1)

        # Right side of first transformer
        #Upper branch
        d += elm.Line().right().length(0.2).at(T.s1)
        d += elm.Line().up().toy(HD.start)
        d += elm.Line().right().length(0.5)
        d += elm.Resistor().right().length(1.5).label('ζeq')
        d += elm.Inductor().right().length(1.5).label('Meq')
        d += elm.Capacitor().right().length(1.5).label('1/Keq')
        Zrad = d.add(elm.ResistorIEC().right().length(2).label('Zrad', 'bottom'))
        #Lower branch
        d += elm.Line().right().length(0.3).at(T.s2)
        d += elm.Line().down().toy(LD.start)
        d += elm.Line().right().tox(Zrad.end)

        #Second transformer
        d += elm.Line().down().toy(T.p1).at(Zrad.end)
        d += elm.Line().right().length(0.3)
        T2 = d.add(elm.Transformer().anchor('p1').label('Aeff:1', 'top'))
        d += elm.Line().left().length(0.3).at(T2.p2)
        d += elm.Line().down().toy(LD.start)
        d += elm.Line().right().length(0.3).at(T2.s2)
        d += elm.Line().down().toy(LD.start)
        d += elm.Line().right().length(1)
        d += elm.Dot(open=True)
        d += elm.Line().right().length(0.3).at(T2.s1)
        d += elm.Line().up().toy(HD.start)
        d += elm.Line().right().length(1)
        d += elm.Dot(open=True).label('P','right')

        return d.draw()

# Individual pMUTs Design
class Individual_Bipolar_pMUT(Basic_pMUT):

    def __init__(self, **kwargs):
            super(Individual_Bipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 1]

    def draw_pMUT(self):
        
        #Draw bottom metal
        self.layout << draw_outer_electrode_single_opening(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6,
                self.orientation
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                self.orientation
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                self.connection_width*1e6, 
                self.orientation+180
                )
        
        #Draw top metal
        self.layout << draw_outer_electrode_single_opening(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6, 
                self.orientation+180
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,  
                self.orientation+180
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.connection_width*1e6, 
                self.orientation
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                        self.pMUT_Radius*1e6,
                        self.back_etch_layer
                        )

        self.layout = pg.union(self.layout, by_layer = True)

        return None

class Individual_Unipolar_pMUT(Basic_pMUT):

    def __init__(self, **kwargs):
            super(Individual_Unipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 0]            

    def draw_pMUT(self):

        #Draw bottom metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                1, 
                electrodes_gap=0
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                self.orientation
                )
        
        #Draw top metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                electrodes_gap=0
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,   
                self.orientation+180
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                self.pMUT_Radius*1e6,
                self.back_etch_layer
                )
        
        self.layout = pg.union(self.layout, by_layer = True)

        return None

# Complete Individual pMUT
class Individual_pMUT():

        def __init__(
            self,
            pMUT_Type = Individual_Unipolar_pMUT(),
            Pad_Size = 150e-6,
            Pad_Distance = 20e-6,
            Routing_Width = 10e-6,
            Routing_Distance = 10e-6,
            Via_Size = 10e-6,
            Label_Size = 30e-6
            ):

           self.layout = Device()
           self.pmut_type = pMUT_Type
           self.pad_size = Pad_Size
           self.pad_distance = Pad_Distance
           self.routing_width = Routing_Width
           self.routing_distance = Routing_Distance
           self.via_size = Via_Size
           self.label_size = Label_Size

        def __repr__(self):
            return self.pmut_type.__repr__()

        @property
        def label(self):
            return 'pMUT radius: '\
                   + str(int(self.pmut_type.pMUT_Radius*1e6)) + ' µm'

        def add_pMUT(self):

            self.pmut_type.draw_pMUT()    
            self.layout<<self.pmut_type.layout

            return self

        def add_pads(self):
            
            #Bottom metal
            BL_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.bottom_metal_layer
                                  ).move([-1e6*self.pad_size*3/2
                                          -1e6*self.routing_distance,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])

            BR_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.bottom_metal_layer
                                  ).move([1e6*self.pad_size/2
                                          +1e6*self.routing_distance,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])
            
            BC_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.bottom_metal_layer
                                  ).move([-1e6*self.pad_size/2,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])

            self.layout << BL_pad
            self.layout << BR_pad
            if type(self.pmut_type) == Individual_Bipolar_pMUT:
                    self.layout << BC_pad

            #Top metal
            BL_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.top_metal_layer
                                  ).move([-1e6*self.pad_size*3/2
                                          -1e6*self.routing_distance,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])

            BR_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.top_metal_layer
                                  ).move([1e6*self.pad_size/2
                                          +1e6*self.routing_distance,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])

            BC_pad = pg.rectangle(size = (1e6*self.pad_size,
                                          1e6*self.pad_size),
                                  layer = self.pmut_type.top_metal_layer
                                  ).move([-1e6*self.pad_size/2,
                                          -1e6*self.pad_size
                                          -1e6*self.pmut_type.pMUT_Radius
                                          -1e6*self.pad_distance
                                          ])

            self.layout << BL_pad
            self.layout << BR_pad
            self.layout << BC_pad

            return self

        def add_routing(self):
            
            #Top metal
            BC_connection = pg.polygon_ports(xpts=[-1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.routing_width/2,
                                                   -1e6*self.routing_width/2], 
                                             ypts=[-1e6*self.pmut_type.pMUT_Radius,
                                                   -1e6*self.pmut_type.pMUT_Radius,
                                                   -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance,
                                                   -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance],
                                             layer = self.pmut_type.top_metal_layer)
            
            TC_connection = pg.polygon_ports(xpts=[-1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.routing_width/2,
                                                   -1e6*self.routing_width/2], 
                                             ypts=[1e6*self.pmut_type.pMUT_Radius,
                                                   1e6*self.pmut_type.pMUT_Radius,
                                                   1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance,
                                                   1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance],
                                             layer = self.pmut_type.top_metal_layer)

            Top_bar_connection = pg.rectangle(size = (1e6*2*self.pad_size
                                                      +1e6*2*self.routing_width
                                                      +1e6*2*self.routing_distance,
                                                      1e6*self.routing_width),
                                              layer = self.pmut_type.top_metal_layer
                                             ).move([-1e6*self.pad_size
                                                     -1e6*self.routing_distance
                                                     -1e6*self.routing_width,
                                                     1e6*self.pmut_type.pMUT_Radius
                                                     +1e6*self.pad_distance
                                                    ])

            Left_bar_connection = pg.rectangle(size = (1e6*self.routing_width,
                                                      1e6*2*self.pmut_type.pMUT_Radius
                                                      +1e6*2*self.pad_distance
                                                      +1e6*self.routing_width),
                                              layer = self.pmut_type.top_metal_layer
                                             ).move([-1e6*self.pad_size
                                                     -1e6*self.routing_distance
                                                     -1e6*self.routing_width,
                                                     -1e6*self.pmut_type.pMUT_Radius
                                                     -1e6*self.pad_distance
                                                    ])

            Right_bar_connection = pg.rectangle(size = (1e6*self.routing_width,
                                                      1e6*2*self.pmut_type.pMUT_Radius
                                                      +1e6*2*self.pad_distance
                                                      +1e6*self.routing_width),
                                              layer = self.pmut_type.top_metal_layer
                                             ).move([1e6*self.pad_size
                                                     +1e6*self.routing_distance,
                                                     -1e6*self.pmut_type.pMUT_Radius
                                                     -1e6*self.pad_distance
                                                    ])

            self.layout << BC_connection
            if type(self.pmut_type) == Individual_Bipolar_pMUT:
                self.layout << TC_connection
                self.layout << Top_bar_connection
                self.layout << Left_bar_connection
                self.layout << Right_bar_connection
        
            #Bottom metal

            BC_connection = pg.polygon_ports(xpts=[-1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.routing_width/2,
                                                   -1e6*self.routing_width/2], 
                                             ypts=[-1e6*self.pmut_type.pMUT_Radius,
                                                   -1e6*self.pmut_type.pMUT_Radius,
                                                   -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance,
                                                   -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance],
                                             layer = self.pmut_type.bottom_metal_layer)

            TC_connection = pg.polygon_ports(xpts=[-1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.pmut_type.connection_width/2,
                                                   +1e6*self.routing_width/2,
                                                   -1e6*self.routing_width/2], 
                                             ypts=[1e6*self.pmut_type.pMUT_Radius,
                                                   1e6*self.pmut_type.pMUT_Radius,
                                                   1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance,
                                                   1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance],
                                             layer = self.pmut_type.bottom_metal_layer)
            
            Top_bar_connection = pg.rectangle(size = (1e6*2*self.pad_size
                                                      +1e6*2*self.routing_width
                                                      +1e6*2*self.routing_distance,
                                                      1e6*self.routing_width),
                                  layer = self.pmut_type.bottom_metal_layer
                                  ).move([-1e6*self.pad_size
                                          -1e6*self.routing_distance
                                          -1e6*self.routing_width,
                                          1e6*self.pmut_type.pMUT_Radius
                                          +1e6*self.pad_distance
                                          ])

            Left_bar_connection = pg.rectangle(size = (1e6*self.routing_width,
                                                      1e6*2*self.pmut_type.pMUT_Radius
                                                      +1e6*2*self.pad_distance
                                                      +1e6*self.routing_width),
                                              layer = self.pmut_type.bottom_metal_layer
                                             ).move([-1e6*self.pad_size
                                                     -1e6*self.routing_distance
                                                     -1e6*self.routing_width,
                                                     -1e6*self.pmut_type.pMUT_Radius
                                                     -1e6*self.pad_distance
                                                    ])

            Right_bar_connection = pg.rectangle(size = (1e6*self.routing_width,
                                                      1e6*2*self.pmut_type.pMUT_Radius
                                                      +1e6*2*self.pad_distance
                                                      +1e6*self.routing_width),
                                              layer = self.pmut_type.bottom_metal_layer
                                             ).move([1e6*self.pad_size
                                                     +1e6*self.routing_distance,
                                                     -1e6*self.pmut_type.pMUT_Radius
                                                     -1e6*self.pad_distance
                                                    ])
            
            if type(self.pmut_type) == Individual_Bipolar_pMUT:
                        self.layout << BC_connection

            self.layout << TC_connection
            self.layout << Top_bar_connection
            self.layout << Left_bar_connection
            self.layout << Right_bar_connection

            return self

        def add_vias(self):

            number_of_vias = 5
            via_height = 10

            via_size = (1e6*self.pad_size
                        -1e6*(1+number_of_vias)*self.routing_distance)\
                        /number_of_vias   

            for i in range(number_of_vias):

                    vias_BL = pg.rectangle(size = (via_size,via_height),
                                           layer = self.pmut_type.vias_layer
                                           ).move([-1e6*1.5*self.pad_size
                                                   +i*(1e6*self.routing_distance
                                                   +via_size),
                                                   -1e6*self.pmut_type.pMUT_Radius
                                                   -1e6*self.pad_distance
                                                   -1e6*self.routing_distance
                                                   -via_height
                                                   ]) 

                    vias_BC = pg.rectangle(size = (via_size,via_height),
                                           layer = self.pmut_type.vias_layer
                                           ).move([-1e6*self.pad_size/2
                                                   +1e6*self.routing_distance
                                                   +i*(1e6*self.routing_distance
                                                   +via_size),
                                                   -1e6*self.pmut_type.pMUT_Radius
                                                   -1e6*self.pad_distance
                                                   -1e6*self.routing_distance
                                                   -via_height
                                                   ])

                    vias_BR = pg.rectangle(size = (via_size,via_height),
                                           layer = self.pmut_type.vias_layer
                                           ).move([1e6*self.pad_size/2
                                                   +2*1e6*self.routing_distance
                                                   +i*(1e6*self.routing_distance
                                                   +via_size),
                                                   -1e6*self.pmut_type.pMUT_Radius
                                                   -1e6*self.pad_distance
                                                   -1e6*self.routing_distance
                                                   -via_height
                                                   ]) 

                    self.layout << vias_BL
                    if type(self.pmut_type) == Individual_Bipolar_pMUT:
                            self.layout << vias_BC 
                    self.layout << vias_BR    

            return self

        def add_label(self):
            
            text = pg.text(self.label, 
                       size = self.label_size*1e6,
                       layer = self.pmut_type.top_metal_layer
                      )

            text.movex(-self.label_size*1e6*6.5)

            text.movey(self.pmut_type.pMUT_Radius*1e6
                       +self.pad_distance*1e6
                       +self.routing_width*1e6
                       +self.label_size/2*1e6)

            self.layout<<text

            return self

        def draw_layout(self):
            
            self.add_pMUT()
            self.add_pads()
            self.add_routing()
            self.add_vias()
            self.add_label()
            self.layout = pg.union(self.layout, by_layer = True)

            return None 

        def draw_save_layout(self, name='test_layout'):
                
            self.draw_layout()
            self.layout.write_gds(name)
                
            return None


# pMUTs Arrays Design
class Parallel_Column_Array_Element_Bipolar_pMUT(Basic_pMUT):

    def __init__(self, **kwargs):
            super(Parallel_Column_Array_Element_Bipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 1]

    def draw_pMUT(self):
        
        #Draw bottom metal
        self.layout << draw_outer_electrode_double_opening(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6,
                90
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_double_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                90
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                0
                )
        
        #Draw top metal
        self.layout << draw_outer_electrode_double_opening(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6, 
                0
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_double_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,  
                0
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,  
                90
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                self.pMUT_Radius*1e6,
                self.back_etch_layer
                )

        return None

    @property
    def lcc(self): #left connection coordinates
        left_down = (-self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        left_up = (-self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [left_down, left_up]

    @property
    def rcc(self): #right connection coordinates
        right_down = (self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        right_up = (self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [right_down, right_up]

    @property
    def ucc(self): #up connection coordinates
        up_left = (-self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        up_right = (self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        return [up_left, up_right]

    @property
    def dcc(self): #down connection coordinates
        down_left = (-self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        down_right = (self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        return [down_left, down_right]

class Parallel_Column_Array_Element_Unipolar_pMUT(Basic_pMUT):

    def __init__(self,**kwargs):
            super(Parallel_Column_Array_Element_Unipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 0]

    def draw_pMUT(self):
        
        #Draw bottom metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                1, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_double_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                90
                )
        
        #Draw top metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_double_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,  
                0
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                self.pMUT_Radius*1e6,
                self.back_etch_layer
                )

        return None

    @property
    def lcc(self): #left connection coordinates
        left_down = (-self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        left_up = (-self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [left_down, left_up]

    @property
    def rcc(self): #right connection coordinates
        right_down = (self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        right_up = (self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [right_down, right_up]

    @property
    def ucc(self): #up connection coordinates
        up_left = (-self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        up_right = (self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        return [up_left, up_right]

    @property
    def dcc(self): #down connection coordinates
        down_left = (-self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        down_right = (self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        return [down_left, down_right]

class Series_Column_Array_Element_Bipolar_pMUT(Basic_pMUT):

    def __init__(self, **kwargs):
            super(Series_Column_Array_Element_Bipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 1]

    def draw_pMUT(self):
        
        #Draw bottom metal
        self.layout << draw_outer_electrode_single_opening(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6,
                self.orientation
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                self.orientation
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                self.connection_width*1e6, 
                self.orientation+180
                )
        
        #Draw top metal
        self.layout << draw_outer_electrode_single_opening(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                1-self.center_electrode_radius_fraction,
                self.electrodes_gap*1e6, 
                self.connection_width*1e6, 
                self.orientation+180
                )

        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                self.electrodes_gap*1e6
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,  
                self.orientation+180
                )

        self.layout << draw_rectangular_attachment(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.connection_width*1e6, 
                self.orientation
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                        self.pMUT_Radius*1e6,
                        self.back_etch_layer
                        )

        self.layout = pg.union(self.layout, by_layer = True)

        return None

    @property
    def lcc(self): #left connection coordinates
        left_down = (-self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        left_up = (-self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [left_down, left_up]

    @property
    def rcc(self): #right connection coordinates
        right_down = (self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        right_up = (self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [right_down, right_up]

    @property
    def ucc(self): #up connection coordinates
        up_left = (-self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        up_right = (self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        return [up_left, up_right]

    @property
    def dcc(self): #down connection coordinates
        down_left = (-self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        down_right = (self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        return [down_left, down_right]

class Series_Column_Array_Element_Unipolar_pMUT(Basic_pMUT):

    def __init__(self, **kwargs):
            super(Series_Column_Array_Element_Unipolar_pMUT, self).__init__(**kwargs)
            self.center_outer_electrode_presence = [1, 0]            

    def draw_pMUT(self):

        #Draw bottom metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.bottom_metal_layer,
                1, 
                electrodes_gap=0
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.bottom_metal_layer, 
                self.connection_width*1e6,  
                self.orientation
                )
        
        #Draw top metal
        self.layout << draw_center_electrode(
                self.pMUT_Radius*1e6,
                self.top_metal_layer,
                self.center_electrode_radius_fraction, 
                electrodes_gap=0
                )

        self.layout << draw_center_electrode_single_trace(
                self.pMUT_Radius*1e6, 
                self.top_metal_layer, 
                self.connection_width*1e6,   
                self.orientation+180
                )
        
        #Draw back etch
        self.layout << draw_etch_pit(
                self.pMUT_Radius*1e6,
                self.back_etch_layer
                )
        
        self.layout = pg.union(self.layout, by_layer = True)

        return None

    @property
    def lcc(self): #left connection coordinates
        left_down = (-self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        left_up = (-self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [left_down, left_up]

    @property
    def rcc(self): #right connection coordinates
        right_down = (self.pMUT_Radius*1e6, -self.connection_width*1e6/2) 
        right_up = (self.pMUT_Radius*1e6, self.connection_width*1e6/2) 
        return [right_down, right_up]

    @property
    def ucc(self): #up connection coordinates
        up_left = (-self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        up_right = (self.connection_width*1e6/2, self.pMUT_Radius*1e6) 
        return [up_left, up_right]

    @property
    def dcc(self): #down connection coordinates
        down_left = (-self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        down_right = (self.connection_width*1e6/2, -self.pMUT_Radius*1e6) 
        return [down_left, down_right]

# Complete Arrays
class Columns_Array():

        def __init__(
            self,
            pMUT_Type = Series_Column_Array_Element_Bipolar_pMUT(),
            Rows = 10,
            Columns = 15,
            Pitch = 200e-6,
            Pad_Size = 100e-6,
            Pad_Distance = 30e-6,
            Routing_Width = 20e-6,
            Routing_Distance = 10e-6,
            Via_Size = 10e-6,
            Label_Size = 50e-6,
            Connected_Columns=False # only for pMUT columns with electrodes connected in parallel
            ):
            
            assert isinstance(Rows, int), 'Rows must be an integer'
            assert isinstance(Columns, int), 'Columns must be an integer'
            assert Pitch>2*pMUT_Type.pMUT_Radius, 'The pitch is must be larger than the pMUTs diameter'
            assert isinstance(pMUT_Type, Basic_pMUT), 'pMUT_Type must be an instance of the class Basic_pMUT'
            assert Pitch>Pad_Size and not(Connected_Columns), 'The pitch is too small to have separate pads'
            assert Via_Size<Routing_Width, 'The vias should be smaller than the connections width'

            self.rows_number = Rows
            self.columns_number = Columns
            self.pitch = Pitch
            self.pmut_type = pMUT_Type
            self.layout = Device()
            self.pad_size = Pad_Size
            self.pad_distance = Pad_Distance
            self.route_width = Routing_Width
            self.route_distance = Routing_Distance
            self.via_size = Via_Size
            self.label_size = Label_Size
            self.connected_columns = Connected_Columns

        
        def __repr__(self):
            
            if (type(self.pmut_type)==Series_Column_Array_Element_Unipolar_pMUT or
                type(self.pmut_type)==Series_Column_Array_Element_Bipolar_pMUT):
                connection_type = 'Series'
                capacitance_info = ('Array capacitance: ' +
                                     '{:.2f}'.format(self.pmut_type.electrical_capacitance*1e12/
                                         self.rows_number*self.columns_number)+
                                     ' pF\n'
                                     )

            if (type(self.pmut_type)==Parallel_Column_Array_Element_Unipolar_pMUT or
                type(self.pmut_type)==Parallel_Column_Array_Element_Bipolar_pMUT):
                connection_type = 'Parallel'

                if self.connected_columns:
                        capacitance_info = ('Array capacitance: ' +
                                            '{:.2f}'.format(self.pmut_type.electrical_capacitance*1e12*
                                            self.rows_number*self.columns_number)+
                                            ' pF\n'
                                           )
                else:
                        capacitance_info = ('Single column capacitance: ' +
                                            '{:.2f}'.format(self.pmut_type.electrical_capacitance*1e12*
                                            self.rows_number)+
                                            ' pF\n'
                                           )

            return (
                    'ARRAY PROPERTIES\n' + 
                    'Array size: ' + str(self.rows_number) + ' by ' + str(self.columns_number) + '\n' +
                    'Column pMUTs connection: ' + connection_type + '\n' +
                    capacitance_info + '\n' + 
                    'pMUT ELEMENT PROPERTIES\n' +
                    self.pmut_type.__repr__() 
                    )

        @property
        def coordinates(self):
            
            coord = np.empty((self.rows_number, self.columns_number), dtype=object)
            
            for i in range(self.rows_number):
                    for j in range(self.columns_number):
                        coord[i,j]=(j*self.pitch*1e6,i*self.pitch*1e6)

            return coord

        @property
        def label(self):

            radius_label = 'pMUTs Radius: '+ str(int(1e6*self.pmut_type.pMUT_Radius))+' µm'
            array_label  = 'Array Size  : '+ str(self.rows_number) + ' x ' +\
                            str(self.columns_number)

            return radius_label + '     ' + array_label

        def add_pMUTs(self):       

            self.pmut_type.draw_pMUT()

            for i in range(self.rows_number):
                    for j in range(self.columns_number):
                        
                        cell = Device() 
                        cell << self.pmut_type.layout

                        if (type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or  
                           type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT):

                                cell.rotate(i*180)

                        cell.move((self.coordinates[i,j][0],\
                                   self.coordinates[i,j][1])
                                 )
                        
                        self.layout<<cell

            return self

        def add_pads(self):
            
            #Bottom pads
            for i in range(self.columns_number):
                    pad = pg.rectangle(size=(self.pad_size*1e6, 
                                             self.pad_size*1e6),
                                             layer=self.pmut_type.top_metal_layer)
                    pad.movey(-self.pmut_type.pMUT_Radius*1e6
                              -self.pad_size*1e6
                              -self.pad_distance*1e6
                              )
                    pad.movex(-self.pad_size*1e6/2
                              +i*self.pitch*1e6
                              )
                    self.layout<<pad

            #Top pads
            for i in range(self.columns_number):
                    pad = pg.rectangle(size=(self.pad_size*1e6, 
                                             self.pad_size*1e6),
                                             layer=self.pmut_type.top_metal_layer)
                    pad.movey(self.pmut_type.pMUT_Radius*1e6
                              +(self.rows_number-1)*self.pitch*1e6
                              +self.pad_distance*1e6
                              )
                    pad.movex(-self.pad_size*1e6/2
                              +i*self.pitch*1e6
                              )
                    self.layout<<pad

            #Fully Connected Columns
            
            if (self.connected_columns or 
                type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or
                type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT):

                    pad = pg.rectangle(size=(self.pad_size/2*1e6
                                             +(self.columns_number-1)*self.pitch*1e6, 
                                             self.pad_size*1e6),
                                       layer=self.pmut_type.top_metal_layer)
                    
                    pad.movey((self.rows_number-1)*self.pitch*1e6
                             +self.pmut_type.pMUT_Radius*1e6
                             +self.pad_distance*1e6)

                    pad.movex(-self.pad_size/2)
                    
                    self.layout<<pad

                    pad = pg.rectangle(size=(self.pad_size/2*1e6
                                             +(self.columns_number-1)*self.pitch*1e6, 
                                             self.pad_size*1e6),
                                       layer=self.pmut_type.top_metal_layer)
                    
                    pad.movey(-self.pmut_type.pMUT_Radius*1e6
                              -self.pad_distance*1e6
                              -self.pad_size*1e6)

                    pad.movex(-self.pad_size/2)
                    
                    self.layout<<pad
                
            # Ground pads    
            if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                type(self.pmut_type) == Parallel_Column_Array_Element_Unipolar_pMUT):

                #Ground pad left
                pad = pg.rectangle(size=(self.pad_size*1e6, 
                                                self.pad_size*1e6
                                                +(self.rows_number-1)*self.pitch*1e6),
                                                layer=self.pmut_type.top_metal_layer)
                pad.movey(-self.pad_size*1e6/2)
                pad.movex(-self.pmut_type.pMUT_Radius*1e6
                        -self.pad_size*1e6
                        -self.pad_distance*1e6
                        )
                self.layout<<pad

                #Ground pad right
                pad = pg.rectangle(size=(self.pad_size*1e6, 
                                        self.pad_size*1e6
                                        +(self.rows_number-1)*self.pitch*1e6),
                                layer=self.pmut_type.top_metal_layer)
                pad.movey(-self.pad_size*1e6/2)
                
                pad.movex(self.pmut_type.pMUT_Radius*1e6
                                +(self.columns_number-1)*self.pitch*1e6
                                +self.pad_distance*1e6
                                )
                self.layout<<pad

            return self

        def add_routing(self):
            
            if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                type(self.pmut_type) == Parallel_Column_Array_Element_Unipolar_pMUT):

                # Routing between columns
                for i in range(self.columns_number-1):
                        for j in range(self.rows_number):
                                trace = pg.rectangle(size=(self.pitch*1e6
                                                        -2*self.pmut_type.pMUT_Radius*1e6
                                                        -2*self.route_distance*1e6,
                                                        self.route_width*1e6),
                                                        layer = self.pmut_type.top_metal_layer)
                                trace.movey(-self.route_width*1e6/2
                                        +self.pitch*1e6*j
                                        )
                                trace.movex(self.pmut_type.pMUT_Radius*1e6
                                        +self.route_distance*1e6
                                        +i*self.pitch*1e6
                                        )
                                self.layout<<trace

                for i in range(self.columns_number-1):
                        for j in range(self.rows_number):
                                trace = pg.rectangle(size=(self.pitch*1e6
                                                        -2*self.pmut_type.pMUT_Radius*1e6
                                                        -2*self.route_distance*1e6,
                                                        self.route_width*1e6),
                                                        layer = self.pmut_type.bottom_metal_layer)
                                trace.movey(-self.route_width*1e6/2
                                        +self.pitch*1e6*j
                                        )
                                trace.movex(self.pmut_type.pMUT_Radius*1e6
                                        +self.route_distance*1e6
                                        +i*self.pitch*1e6
                                        )
                                self.layout<<trace

            # Routing between rows
            for i in range(self.rows_number-1):
                    for j in range(self.columns_number):
                        trace = pg.rectangle(size=(self.route_width*1e6,
                                                   self.pitch*1e6
                                                -2*self.pmut_type.pMUT_Radius*1e6
                                                -2*self.route_distance*1e6),
                                                layer = self.pmut_type.top_metal_layer)
                        trace.movex(-self.route_width*1e6/2
                                    +self.pitch*1e6*j
                                   )
                        trace.movey(self.pmut_type.pMUT_Radius*1e6
                                   +self.route_distance*1e6
                                   +i*self.pitch*1e6
                                   )
                        self.layout<<trace

            for i in range(self.rows_number-1):
                    for j in range(self.columns_number):
                        trace = pg.rectangle(size=(self.route_width*1e6,
                                                   self.pitch*1e6
                                                -2*self.pmut_type.pMUT_Radius*1e6
                                                -2*self.route_distance*1e6),
                                                layer = self.pmut_type.bottom_metal_layer)
                        trace.movex(-self.route_width*1e6/2
                                    +self.pitch*1e6*j
                                   )
                        trace.movey(self.pmut_type.pMUT_Radius*1e6
                                   +self.route_distance*1e6
                                   +i*self.pitch*1e6
                                   )

                        if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                            type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or
                            type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT):

                                self.layout<<trace

            # Trapezoidal connections to routings and pads
            for i in range(self.rows_number):
                    for j in range(self.columns_number):
                            
                            #define [x,y] coordinates of routing connections 
                            bl = [-1e6*self.route_width/2,
                                  -1e6*self.pmut_type.pMUT_Radius-1e6*self.route_distance] #bottom left
                            br = [1e6*self.route_width/2,
                                  -1e6*self.pmut_type.pMUT_Radius-1e6*self.route_distance]#bottom rigth
                            rb = [1e6*self.pmut_type.pMUT_Radius+1e6*self.route_distance,
                                  -1e6*self.route_width/2]#right bottom
                            rt = [1e6*self.pmut_type.pMUT_Radius+1e6*self.route_distance,
                                  1e6*self.route_width/2]#right top
                            tl = [-1e6*self.route_width/2,
                                  1e6*self.pmut_type.pMUT_Radius+1e6*self.route_distance]#top left
                            tr = [1e6*self.route_width/2,
                                  1e6*self.pmut_type.pMUT_Radius+1e6*self.route_distance]#top right
                            lb = [-1e6*self.pmut_type.pMUT_Radius-1e6*self.route_distance,
                                  -1e6*self.route_width/2]#left bottom
                            lt = [-1e6*self.pmut_type.pMUT_Radius-1e6*self.route_distance,
                                  1e6*self.route_width/2]#left top
                            
                            if j==0:
                                    lb = [-1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance,
                                          -1e6*self.pad_size/2]
                                    lt = [-1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance,
                                          1e6*self.pad_size/2]

                            if i==0:
                                    bl = [-1e6*self.pad_size/2,
                                          -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance] 
                                    br = [1e6*self.pad_size/2,
                                          -1e6*self.pmut_type.pMUT_Radius-1e6*self.pad_distance]

                            if j==(self.columns_number-1):
                                    rb = [1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance,
                                        -1e6*self.pad_size/2]
                                    rt = [1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance,
                                        1e6*self.pad_size/2]

                            if i==(self.rows_number-1):
                                    tl = [-1e6*self.pad_size/2,
                                          1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance]#top left
                                    tr = [1e6*self.pad_size/2,
                                          1e6*self.pmut_type.pMUT_Radius+1e6*self.pad_distance]

                            bottom_connection=pg.polygon_ports(xpts=[bl[0],
                                                                     br[0],
                                                                     self.pmut_type.dcc[1][0],
                                                                     self.pmut_type.dcc[0][0]], 
                                                               ypts=[bl[1],
                                                                     br[1],
                                                                     self.pmut_type.dcc[1][1],
                                                                     self.pmut_type.dcc[0][1]],
                                                                     layer = self.pmut_type.top_metal_layer)

                            right_connection=pg.polygon_ports(xpts=[rb[0],
                                                                    rt[0],
                                                                    self.pmut_type.rcc[1][0],
                                                                    self.pmut_type.rcc[0][0]], 
                                                               ypts=[rb[1],
                                                                     rt[1],
                                                                     self.pmut_type.rcc[1][1],
                                                                     self.pmut_type.rcc[0][1]],
                                                                     layer = self.pmut_type.bottom_metal_layer)

                            top_connection=pg.polygon_ports(xpts=[tl[0],
                                                                  tr[0],
                                                                  self.pmut_type.ucc[1][0],
                                                                  self.pmut_type.ucc[0][0]], 
                                                           ypts=[tl[1],
                                                                 tr[1],
                                                                 self.pmut_type.ucc[1][1],
                                                                 self.pmut_type.ucc[0][1]],
                                                                 layer = self.pmut_type.top_metal_layer)

                            left_connection=pg.polygon_ports(xpts=[lb[0],
                                                                   lt[0],
                                                                   self.pmut_type.lcc[1][0],
                                                                   self.pmut_type.lcc[0][0]], 
                                                            ypts=[lb[1],
                                                                  lt[1],
                                                                  self.pmut_type.lcc[1][1],
                                                                  self.pmut_type.lcc[0][1]],
                                                                  layer = self.pmut_type.bottom_metal_layer)

                            # Removes default ports
                            bottom_connection.ports = {}
                            right_connection.ports = {}
                            top_connection.ports = {}
                            left_connection.ports = {}
                            
                            #move connections to pMUT location
                            bottom_connection.move((self.coordinates[i,j][0],\
                                   self.coordinates[i,j][1]))
                            right_connection.move((self.coordinates[i,j][0],\
                                   self.coordinates[i,j][1]))
                            top_connection.move((self.coordinates[i,j][0],\
                                   self.coordinates[i,j][1]))
                            left_connection.move((self.coordinates[i,j][0],\
                                   self.coordinates[i,j][1]))

                            self.layout<<bottom_connection
                            self.layout<<top_connection

                            if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                                type(self.pmut_type) == Parallel_Column_Array_Element_Unipolar_pMUT):

                                self.layout<<left_connection
                                self.layout<<right_connection

                            #A superimposing trapezoidal connection on the other metal is added
                            #left/right only for the unipolar case and on all sides for the 
                            # bipolar case
                             
                            right_connection=pg.polygon_ports(xpts=[rb[0],
                                                                    rt[0],
                                                                    self.pmut_type.rcc[1][0],
                                                                    self.pmut_type.rcc[0][0]], 
                                                               ypts=[rb[1],
                                                                     rt[1],
                                                                     self.pmut_type.rcc[1][1],
                                                                     self.pmut_type.rcc[0][1]],
                                                                     layer = self.pmut_type.top_metal_layer)

                            left_connection=pg.polygon_ports(xpts=[lb[0],
                                                                   lt[0],
                                                                   self.pmut_type.lcc[1][0],
                                                                   self.pmut_type.lcc[0][0]], 
                                                            ypts=[lb[1],
                                                                  lt[1],
                                                                  self.pmut_type.lcc[1][1],
                                                                  self.pmut_type.lcc[0][1]],
                                                                  layer = self.pmut_type.top_metal_layer)
                            
                            right_connection.ports = {}
                            left_connection.ports = {}

                            right_connection.move((self.coordinates[i,j][0],\
                                                   self.coordinates[i,j][1]))
                            left_connection.move((self.coordinates[i,j][0],\
                                                  self.coordinates[i,j][1]))

                            if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                                type(self.pmut_type) == Parallel_Column_Array_Element_Unipolar_pMUT):

                                self.layout<<left_connection
                                self.layout<<right_connection

                            if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                                type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or
                                type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT):
                                    
                                    bottom_connection=pg.polygon_ports(xpts=[bl[0],
                                                                     br[0],
                                                                     self.pmut_type.dcc[1][0],
                                                                     self.pmut_type.dcc[0][0]], 
                                                               ypts=[bl[1],
                                                                     br[1],
                                                                     self.pmut_type.dcc[1][1],
                                                                     self.pmut_type.dcc[0][1]],
                                                                     layer = self.pmut_type.bottom_metal_layer)

                                    top_connection=pg.polygon_ports(xpts=[tl[0],
                                                                  tr[0],
                                                                  self.pmut_type.ucc[1][0],
                                                                  self.pmut_type.ucc[0][0]], 
                                                           ypts=[tl[1],
                                                                 tr[1],
                                                                 self.pmut_type.ucc[1][1],
                                                                 self.pmut_type.ucc[0][1]],
                                                                 layer = self.pmut_type.bottom_metal_layer)

                                    # Removes default ports
                                    bottom_connection.ports = {}
                                    top_connection.ports = {}
                                        
                                    #move connections to pMUT location
                                    bottom_connection.move((self.coordinates[i,j][0],\
                                                            self.coordinates[i,j][1]))
                                    top_connection.move((self.coordinates[i,j][0],\
                                                         self.coordinates[i,j][1]))

                                    self.layout<<bottom_connection
                                    self.layout<<top_connection

            return self

        def add_vias(self):
                
                for i in range(self.rows_number):
                    for j in range(self.columns_number):
                                bottom_via = pg.rectangle(size=(1e6*self.via_size,
                                                                1e6*self.via_size),
                                                                layer = self.pmut_type.vias_layer)

                                bottom_via.move([-1e6*self.via_size/2
                                                 +j*1e6*self.pitch,
                                                -1e6*self.pmut_type.pMUT_Radius
                                                -1e6*self.via_size
                                                -1e6*self.route_distance
                                                +i*1e6*self.pitch])

                                top_via = pg.rectangle(size=(1e6*self.via_size,
                                                        1e6*self.via_size),
                                                        layer = self.pmut_type.vias_layer)

                                top_via.move([-1e6*self.via_size/2
                                              +j*1e6*self.pitch,
                                              1e6*self.pmut_type.pMUT_Radius
                                              +1e6*self.route_distance
                                              +i*1e6*self.pitch])

                                if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                                    type(self.pmut_type) == Parallel_Column_Array_Element_Unipolar_pMUT):

                                        right_via = pg.rectangle(size=(1e6*self.via_size,
                                                                1e6*self.via_size),
                                                                layer = self.pmut_type.vias_layer)

                                        right_via.move([1e6*self.pmut_type.pMUT_Radius
                                                        +1e6*self.route_distance
                                                        +j*1e6*self.pitch,
                                                        -1e6*self.via_size/2
                                                        +i*1e6*self.pitch])              

                                        left_via = pg.rectangle(size=(1e6*self.via_size,
                                                                1e6*self.via_size),
                                                                layer = self.pmut_type.vias_layer)

                                        left_via.move([-1e6*self.pmut_type.pMUT_Radius
                                                -1e6*self.via_size
                                                -1e6*self.route_distance
                                                +j*1e6*self.pitch,
                                                -1e6*self.via_size/2
                                                +i*1e6*self.pitch])

                                        self.layout<<right_via
                                        self.layout<<left_via

                                if (type(self.pmut_type) == Parallel_Column_Array_Element_Bipolar_pMUT or
                                    type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or
                                    type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT):
                                        
                                        self.layout<<top_via
                                        self.layout<<bottom_via
                return self

        def add_label(self):

                text = pg.text(self.label, 
                               size = self.label_size*1e6,
                               layer = self.pmut_type.top_metal_layer
                               )
                text.movey(-self.pad_distance*1e6
                           -self.pmut_type.pMUT_Radius*1e6
                           -self.pad_size*1e6
                           -1.5*self.label_size*1e6)

                self.layout<<text

                return self                        

        def draw_layout(self):

                self.add_pMUTs()
                self.add_pads()
                self.add_routing()
                self.add_vias()
                self.add_label()
                self.layout = pg.union(self.layout, by_layer = True)

                return None 

        def draw_save_layout(self, name='test_layout'):
                
                self.draw_layout()
                self.layout.write_gds(name)
                
                return None

        def plot_directivity_function(self,
                                      thetas = np.linspace(-np.pi/2,np.pi/2,300),
                                      phi = 0,
                                      phases = None,
                                      frequency = None,
                                      farfield = 0.5
                                      ):

            #Passing self arguments as default arguments
            if (type(self.pmut_type) == Series_Column_Array_Element_Unipolar_pMUT or
                type(self.pmut_type) == Series_Column_Array_Element_Bipolar_pMUT or
                self.connected_columns):
                phases = None

            if phases is None:
                phases = np.zeros(self.columns_number)
            else:
                phases = np.array(phases)

            if frequency is None:
                frequency = self.pmut_type.resonance_frequency

            # phases array dimensions check depending on the type of array
            if self.connected_columns:
                if len(phases)>1:
                   warnings.warn('pMUTs are all connected in parallel, the displayed\
                                  radiation pattern assumes they are all in phase')
            else:
                if len(phases)!=self.columns_number:
                   raise Exception('the phases array should have the same length as\
                                    the number of columns')
            
            #Pre-processing to change the pMUT coordinates and place 
            #the radiation pattern center at the center of the array 
            originX = self.pitch*(self.columns_number-1)/2
            originY = self.pitch*(self.rows_number-1)/2
            coordX   = np.zeros((self.rows_number,self.columns_number)) 
            coordY   = np.zeros((self.rows_number,self.columns_number)) 
            ph       = np.array([])
            
            #assign coordinates and phases to flattened arrays
            for i in range(self.rows_number):
                    for j in range(self.columns_number):

                        coordX[i,j] = 1e-6*self.coordinates[i,j][0]-originX
                        coordY[i,j] = 1e-6*self.coordinates[i,j][1]-originY
             
            for j in range(self.columns_number):
                    if not(self.connected_columns):
                           ph = np.append(ph,phases[j]*np.ones(self.rows_number))

            coordX = np.transpose(coordX).flatten()
            coordY = np.transpose(coordY).flatten()
            coord = np.array([coordX, coordY]).transpose()

            #Plotting
            b = Array_DF(thetas, 
                         phi, 
                         coord, 
                         ph, 
                         self.pmut_type.medium_speed, 
                         self.pmut_type.pMUT_Radius,
                         frequency,
                         farfield)

            ax = plt.subplot(111, polar=True)
            ax.scatter(thetas, b)
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            ax.set_thetamin(-90)
            ax.set_thetamax(90)
            ax.set_title('Array Directivity Function')
            
            return None


# Triangular Microphone
class Triangular_Cantilever_Mic:

    def __init__(
        self, 
        Side=500e-6, 
        BottomMetalLayer=0,
        ViasLayer=1,
        TopMetalLayer=2,
        BackEtchLayer=3, 
        EtchThroughLayer=4,
        CutWidth = 2e-6,
        BottomElectrodeOffset=2e-6,
        TopElectrodeOffset=3e-6,
        RoutingWidth=10e-6,
        RoutingDistance=10e-6,
        PadSize=150e-6,
        PiezoMaterial = ScAlN(36),
        PiezoThickness = 500e-6,
        LabelSize = 30e-6
        ):    
        
        self.layout = Device()
        self.side = Side
        self.bottom_metal_layer = BottomMetalLayer
        self.vias_layer = ViasLayer
        self.top_metal_layer = TopMetalLayer
        self.back_etch_layer = BackEtchLayer
        self.etch_through_layer = EtchThroughLayer
        self.cut_width = CutWidth      
        self.bottom_electrode_offset = BottomElectrodeOffset
        self.top_electrode_offset = TopElectrodeOffset
        self.routing_width = RoutingWidth
        self.routing_distance = RoutingDistance
        self.pad_size = PadSize
        self.pad_distance = 3*self.routing_distance+2*self.routing_width
        self.piezo_material = PiezoMaterial
        self.piezo_thickness = PiezoThickness
        self.label_size = LabelSize
    
    @property
    def electrical_capacitance(self):
            
        return 8.854e-12*self.piezo_material['epsilon_r']*\
               self.side**2/self.piezo_thickness

    @property
    def label(self):
        return 'Side length: '+ str(int(self.side*1e6)) + ' µm'

    def add_electrodes(self):

        #Bottom Electrode
        self.layout<<triangular_electrodes(self.side*1e6,
                                           self.cut_width/2*1e6
                                           +self.bottom_electrode_offset*1e6, 
                                           self.bottom_metal_layer)
        
        #Top Electrode
        self.layout<<triangular_electrodes(self.side*1e6,
                                           self.cut_width/2*1e6
                                           +self.top_electrode_offset*1e6, 
                                           self.top_metal_layer)

        return self

    def add_cuts(self):

        self.layout<<triangular_mic_slots(self.side*1e6, 
                                          self.cut_width*1e6,
                                          self.etch_through_layer)

        return self

    def add_backetch(self):
        
        self.layout<<triangular_mic_back_etch(1e6*self.side, 
                                              self.back_etch_layer)
        
        return self

    def add_routing(self):
        
        #Bottom routing - TopRight corner
        connection1 = pg.rectangle(size=(1e6*self.routing_width,1e6*self.routing_distance),
                                    layer = self.bottom_metal_layer)

        connection1<< pg.rectangle(size=(1e6*2*self.routing_distance+1e6*2*self.routing_width,
                                         1e6*self.routing_width),
                                   layer = self.bottom_metal_layer
                                   ).movey(1e6*self.routing_distance)

        connection1<< pg.rectangle(size=(1e6*self.routing_width,
                                         1e6*2*self.routing_distance+1e6*2*self.routing_width),
                                   layer = self.bottom_metal_layer
                                   ).move([1e6*2*self.routing_distance+1e6*self.routing_width,
                                           -1e6*self.routing_distance-1e6*self.routing_width])

        connection1<< pg.rectangle(size=(1e6*self.routing_distance,1e6*self.routing_width),
                                    layer = self.bottom_metal_layer
                                   ).move([1e6*self.routing_distance+1e6*self.routing_width,
                                           -1e6*self.routing_width-1e6*self.routing_distance])
        
        #Move it to the TopRight corner
        connection1.movex(-1e6*self.routing_width-1e6*self.routing_distance+1e6*self.side/2)
        connection1.movey(1e6*self.side/2)

        #Generate the other 3 connections
        connection2 = Device()
        connection2 << connection1
        connection2.rotate(90)

        connection3 = Device()
        connection3 << connection2
        connection3.rotate(90)

        connection4 = Device()
        connection4 << connection3
        connection4.rotate(90)

        self.layout<<connection1
        self.layout<<connection2
        self.layout<<connection3
        self.layout<<connection4

        #Top routing - TopRight corner
        connection5 = pg.rectangle(size=(1e6*self.routing_width,1e6*2*self.routing_distance+1e6*self.routing_width),
                                    layer = self.top_metal_layer)

        connection5<< pg.rectangle(size=(1e6*4*self.routing_distance+1e6*4*self.routing_width,
                                         1e6*self.routing_width),
                                   layer = self.top_metal_layer
                                   ).movey(1e6*2*self.routing_distance+1e6*self.routing_width)

        connection5<< pg.rectangle(size=(1e6*self.routing_width,
                                         1e6*4*self.routing_distance+1e6*4*self.routing_width),
                                   layer = self.top_metal_layer
                                   ).move([1e6*4*self.routing_distance+1e6*3*self.routing_width,
                                           -1e6*2*self.routing_distance-1e6*2*self.routing_width])

        connection5<< pg.rectangle(size=(1e6*self.routing_width+1e6*2*self.routing_distance,1e6*self.routing_width),
                                    layer = self.top_metal_layer
                                   ).move([1e6*2*self.routing_distance+1e6*2*self.routing_width,
                                           -1e6*2*self.routing_width-1e6*2*self.routing_distance])

        #Move it to the TopRight corner
        connection5.movex(-1e6*2*self.routing_width-1e6*2*self.routing_distance+1e6*self.side/2)
        connection5.movey(1e6*self.side/2)

        #Generate the other 3 connections
        connection6 = Device()
        connection6 << connection5
        connection6.rotate(90)

        connection7 = Device()
        connection7 << connection6
        connection7.rotate(90)

        connection8 = Device()
        connection8 << connection7
        connection8.rotate(90)

        self.layout<<connection5
        self.layout<<connection6
        self.layout<<connection7
        self.layout<<connection8

        return self 

    def add_pads(self):   

        #Bottom metal
        connection_UL = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.bottom_metal_layer
                                      ).move([-1e6*(self.side/2-
                                                    2*self.routing_width-
                                                    3*self.routing_distance),
                                              1e6*self.side/2])

        connection_UR = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(self.side/2-
                                                   3*self.routing_width-
                                                   3*self.routing_distance),
                                              1e6*self.side/2])

        connection_BL = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.bottom_metal_layer
                                      ).move([-1e6*(self.side/2-
                                                    2*self.routing_width-
                                                    3*self.routing_distance),
                                              -1e6*self.side/2-1e6*self.pad_distance])

        connection_BR = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(self.side/2-
                                                   3*self.routing_width-
                                                   3*self.routing_distance),
                                              -1e6*self.side/2-1e6*self.pad_distance])

        bottom_pad_UL = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(-self.side/2-
                                                   2*self.routing_width-
                                                   2*self.routing_distance),
                                              1e6*self.side/2+1e6*self.pad_distance])

        bottom_pad_UR = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(self.pad_size/2+
                                                   self.routing_distance),
                                              1e6*self.side/2+1e6*self.pad_distance])

        bottom_pad_BL = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(-self.side/2-
                                                   2*self.routing_width-
                                                   2*self.routing_distance),
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance-
                                              1e6*self.pad_size])

        bottom_pad_BR = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.bottom_metal_layer
                                      ).move([1e6*(self.pad_size/2+
                                                   self.routing_distance),
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance-
                                              1e6*self.pad_size])
        
        self.layout<<connection_UL
        self.layout<<connection_UR
        self.layout<<connection_BL
        self.layout<<connection_BR

        self.layout<<bottom_pad_UL
        self.layout<<bottom_pad_UR
        self.layout<<bottom_pad_BL
        self.layout<<bottom_pad_BR

        #Top metal
        connection_UC = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.top_metal_layer
                                      ).move([-1e6*self.routing_width/2,
                                              1e6*self.side/2])

        connection_BC = pg.rectangle(size = (1e6*self.routing_width,
                                             1e6*self.pad_distance),
                                      layer = self.top_metal_layer
                                      ).move([-1e6*self.routing_width/2,
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance])  

        top_pad_UL = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([1e6*(-self.side/2-
                                                   2*self.routing_width-
                                                   2*self.routing_distance),
                                              1e6*self.side/2+1e6*self.pad_distance])

        top_pad_UR = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([1e6*(self.pad_size/2+
                                                   self.routing_distance),
                                              1e6*self.side/2+1e6*self.pad_distance])

        top_pad_BL = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([1e6*(-self.side/2-
                                                   2*self.routing_width-
                                                   2*self.routing_distance),
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance-
                                              1e6*self.pad_size])

        top_pad_BR = pg.rectangle(size = (1e6*self.side/2+
                                             1e6*self.routing_distance+
                                             1e6*2*self.routing_width-
                                             1e6*self.pad_size/2,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([1e6*(self.pad_size/2+
                                                   self.routing_distance),
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance-
                                              1e6*self.pad_size])

        top_pad_BC = pg.rectangle(size = (1e6*self.pad_size,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([-1e6*self.pad_size/2,
                                              -1e6*self.side/2-
                                              1e6*self.pad_distance-
                                              1e6*self.pad_size])

        top_pad_UC = pg.rectangle(size = (1e6*self.pad_size,
                                             1e6*self.pad_size),
                                      layer = self.top_metal_layer
                                      ).move([-1e6*self.pad_size/2,
                                              1e6*self.side/2+
                                              1e6*self.pad_distance])

        self.layout<<connection_UC
        self.layout<<connection_BC

        self.layout<<top_pad_UL
        self.layout<<top_pad_UR
        self.layout<<top_pad_BL
        self.layout<<top_pad_BR
        self.layout<<top_pad_BC
        self.layout<<top_pad_UC

        return self

    def add_vias(self):

        number_of_vias = 5
        via_height = 10

        via_size = (1e6*self.side/2-
                    1e6*number_of_vias*self.routing_distance+
                    1e6*2*self.routing_width-
                    1e6*self.pad_size/2)\
                    /number_of_vias   

        for i in range(number_of_vias):

                vias_UL = pg.rectangle(size = (via_size,via_height),
                                        layer = self.vias_layer
                                        ).move([1e6*(-self.side/2-
                                                self.routing_width-
                                                2*self.routing_distance+
                                                i*(self.routing_distance+
                                                1e-6*via_size)
                                                ),
                                                1e6*self.side/2+
                                                1e6*self.routing_distance+
                                                1e6*self.pad_distance])

                vias_UR = pg.rectangle(size = (via_size,via_height),
                                        layer = self.vias_layer
                                        ).move([1e6*(self.pad_size/2+
                                                2*self.routing_distance+
                                                i*(self.routing_distance+
                                                1e-6*via_size)
                                                ),
                                                1e6*self.side/2+
                                                1e6*self.routing_distance+
                                                1e6*self.pad_distance])

                vias_BL = pg.rectangle(size = (via_size,via_height),
                                        layer = self.vias_layer
                                        ).move([1e6*(-self.side/2-
                                                self.routing_width-
                                                2*self.routing_distance+
                                                i*(self.routing_distance+
                                                1e-6*via_size)
                                                ),
                                                -1e6*self.side/2-
                                                1e6*self.routing_distance-
                                                1e6*self.pad_distance-
                                                via_height])

                vias_BR = pg.rectangle(size = (via_size,via_height),
                                        layer = self.vias_layer
                                        ).move([1e6*(self.pad_size/2+
                                                2*self.routing_distance+
                                                i*(self.routing_distance+
                                                1e-6*via_size)
                                                ),
                                                -1e6*self.side/2-
                                                1e6*self.routing_distance-
                                                1e6*self.pad_distance-
                                                via_height])

                self.layout<<vias_UL
                self.layout<<vias_UR
                self.layout<<vias_BL
                self.layout<<vias_BR

        return self

    def add_label(self):
        
        text = pg.text(self.label, 
                       size = self.label_size*1e6,
                       layer = self.top_metal_layer
                      )

        text.movex(-self.label_size*1e6*6
                  )

        text.movey(-self.side/2*1e6
                   -self.pad_distance*1e6
                   -self.pad_size*1e6
                   -2*self.label_size*1e6)

        self.layout<<text

        return self

    def draw_layout(self):

        self.add_electrodes()
        self.add_cuts()
        self.add_backetch()
        self.add_routing()
        self.add_pads()
        self.add_vias()
        self.add_label()
        self.layout = pg.union(self.layout, by_layer = True)

        return self.layout

    def draw_save_layout(self, name='test_layout'):
        
        self.draw_layout()
        self.layout.write_gds(name)
                
        return qp(self.layout)