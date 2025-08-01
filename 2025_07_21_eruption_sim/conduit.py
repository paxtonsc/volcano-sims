''' Sample input file for conduit in Vulcanian eruption. See below for instructions. '''

import run_globals
import numpy as np

# SolutionOrder = 2 means that the solution is piecewise quadratic.
Numerics = {'AVParameter': 30,
 'ApplyLimiters': 'PositivityPreservingMultiphasevpT',
 'ArtificialViscosity': True,
 'ElementQuadrature': 'GaussLegendre',
 'FaceQuadrature': 'GaussLegendre',
 'L2InitialCondition': True,
 'SolutionBasis': 'LagrangeSeg',
 'SolutionOrder': 0,
 'Solver': 'DG'}


# Specify the mesh here. The domain goes from xmin to xmax in NumElemsX
# elements.
Mesh = {'ElementShape': 'Segment',
 'File': None,
 'NumElemsX': 1000,
 'xmax': 0,   # Top of the conduit is placed at -150 m for 2D coupling, but you can choose whatever if not coupling to 2D.
 'xmin': -1000.0
}


Physics = {'ConvFluxNumerical': 'LaxFriedrichs',
 'Gas1': {'R': 287.0, 'gamma': 1.4},
 'Gas2': {'R': 461.3762486126526, 'c_p': 2288.0},
 'Liquid': {'E_m0': 0,                # Ground state energy (0 is fine)
            'K': 10000000000.0,       # Bulk modulus of magma phase (Pa)
            'c_m': 1000.0,            # Heat capacity of magma phase (J / kg K)
            'p0': 35999999.99999999,  # Linearized magma EOS constant
            'rho0': 2600.0},          # Linearized magma EOS constant
 'Solubility': {'k': 2.8e-06, 'n': 0.5}, # Henry's law: max soluble mass fraction is = (k * p ** n) * yL
 'Type': 'MultiphasevpT',
 'Viscosity': {'mu0': 300000.0},      # Unused option (viscosity is specified in the source code for FrictionVolFracVariableMu)
 'tau_d': Exception("Deprecated")     # Unused. For exsolution timescale (s), see SourceTerms below.
}


# Set common parameters
p_chamber =  20627361.74
T_chamber = 950 + 273.15 # 1223.15
yC = 0.4   # Crystal mass fraction
yWt = 0.008 # Total water mass fraction


chi_water = (1.0 - yC) * yWt / (1 - yWt)
radio = 5
f_plug = 1.9e8
len_plug = 50
t_plug = f_plug / (2*np.pi*radio*len_plug)
trac_par = 2*t_plug/radio


SourceTerms = {'source1': {'Function': 'GravitySource', # Gravity
             'gravity': 9.8,
             'source_treatment': 'Explicit'},
 'source2': {'Function': 'FrictionVolFracVariableMu',   # Friction source term for given conduit radius
             'use_default_viscosity': True,
             'default_viscosity': 1e6,
             'conduit_radius': 5.,
             'viscosity_factor': 1/5,#1/20,
             'source_treatment': 'Explicit',
             'model_plug': True,
             'plug_boundary_0': -50,
             },
 'source3': {'Function': 'ExsolutionSource',            # Exsolution source
             'source_treatment': 'Explicit',
             'tau_d': 0.01,                             # Exsolution timescale
             },
 "source4": {
        "Function": "FragmentationStrainRateSource",
        "tau_f": 0.01, # This is the fragmentation timescale
        "G": 1e10,
        "k": 0.0001,
        "fragsmooth_scale": 0.010, # This is the smoothing scale that was also in the volume fraction source
        "which_criterion": "both", # This is either "shear" for only using shear strain to determine whether the magma fragments, "tensile" for only using the longitudinal (du/dz) strain rate, or "both" for using the larger of the two strain rates
        "conduit_radius": 5, # Conduit radius (make sure this is the same as what you provide to the Friction source term!
        'source_treatment': 'Explicit',
    },
"slip_source": {
    "Function": "SlipSource",
    "source_treatment": "Explicit",
    
    },
"source5": {
  "Function": "FrictionVolSlip",
  "source_treatment": "Explicit", 
  "conduit_radius": 5, 
  "tau_peak": t_plug,
  "tau_r": 0  ,
  "D_c": 0.01,
  "plug_boundary_0": -50,
  "use_constant_tau": True,
  "exponential_tau": False,
},
}


Output = {'AutoPostProcess': False,
 'Prefix': 'tungurahua_atm_1m',              # Output filename
 'WriteInitialSolution': True,
 'WriteInterval': run_globals.write_interval,                                   # Output frequency (this many timesteps pass before file is written)
}




# Define the cosine taper function
def cosine_taper(x, x1, x2, y1, y2):
    return np.where(x < x1, y1,
                    np.where(x > x2, y2,
                             y1 + (y2 - y1) * 0.5 * (1 - np.cos(np.pi * (x - x1) / (x2 - x1)))))
 # Define the transition region
x1 = -len_plug - 10  # Start of transition
x2 = -len_plug + 10  # End of transition


# Initial condition parameters
InitialCondition = {'Function': 'StaticPlug',          # Specify to call physics/multiphasevpT/functions > SteadyState as initial condition
 'p_chamber': p_chamber,
  # Magma properties
 'K_magma': 1e10,
 'c_v_magma': 1000.0,
 'neglect_edfm': True,
 'p0_magma': 35999999.99999999,
 'rho0_magma': 2600.0,
 # Solubility properties
 'solubility_k': Physics["Solubility"]["k"],
 'solubility_n': Physics["Solubility"]["n"],
 'x_global': np.linspace(Mesh["xmin"], Mesh["xmax"], Mesh["NumElemsX"]), # All x meshes linked together


 # Define the functions using cosine taper
 'traction_fn': lambda x: cosine_taper(x, x1, x2, 0, -trac_par),
 'yWt_fn': lambda x: cosine_taper(x, x1, x2, yWt, 1e-5),
 'yC_fn': lambda x: cosine_taper(x, x1, x2, yC, 0.95),
 # 'yF_fn': lambda x: cosine_taper(x, x1, x2, 0, 1),
 'T_fn': lambda x: cosine_taper(x, x1, x2, T_chamber, T_chamber - 20),

 # p_vent is set slightly above 1e5 to make sure flow is outflow
 'enforce_p_vent': None,                    # If not None, Scales traction_fn iteratively so that vent pressure has this value
 # 'is_solve_direction_downward': False,
 }


# This is needed by Quail. This is not the exact solution, just something callable.
ExactSolution = {'Function': 'RiemannProblem'}


# Set boundary conditions here.
#BoundaryConditions = {
#  'x1': {'BCType': 'SlipWall',   # Inlet boundary condition
                          # Henry's law exponent
#},
    
BoundaryConditions = {
  'x1': {'BCType': 'PressureStableLinearizedInlet1D',   # Inlet boundary condition
        'T_chamber': T_chamber,                         # Chamber temperature
        'cVFamp': 0.0,                                  # Amplitude of crystal fraction variation
        'cVFav': yC,                                    # Crystal fraction input
        'chi_water': chi_water,                         # Mass concentration of water
        'cos_freq': 0.25,                               # Frequency of crys    tal fraction variation
        'is_gaussian': False,                           # If true, use Gaussian variation instead of cosine
        'p_chamber': p_chamber,                         # Chamber pressure
        'trace_arho': 2.6000000000000002e-05,           # Trace density (for numerical stability)
        # 'approx_mass_fracs': False,                   # Compute exactly the mass fraction of inlet fluid
        # 'solubility_k': 5e-06,                        # Henry's law coefficient
        # 'solubility_n': 0.5,                          # Henry's law exponent
  },
  # For running in parallel, the following boundary condition type is needed:
  'x2': {
      'BCType': 'MultiphasevpT2D1D', 
      'bkey': 'vent'
  }
  # For running serial, using p boundary condition:
  #'x2': {'BCType': 'PressureOutlet1D',                   # Pressure outlet boundary condition (automatically chokes if needed)
  #       'p': 1e5,                  # Boundary pressure (if flow not choked) -- with scale height factor
  #       },
}


# Linked parallel solvers. If running in serial, leave as empty list.
# LinkedSolvers = []
# A parallel solver would have multiple input files like this one, and they would
# be linked as follows:
LinkedSolvers = [
    {
        'BoundaryName': 'vent',
        'DeckName': 'vent_region.py'
    },
]

TimeStepping = {'FinalTime': 10, # Final 
 'InitialTime': 0.0,
 'NumTimeSteps': 60000,# Number of timesteps to run for
 'TimeStepper': 'RK3SR', # 'FE', # 'RK3SR',  # 4-step RK3 scheme that maximizes CFL stability region per function eval
}