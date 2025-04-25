''' Sample input file for conduit in Vulcanian eruption. See below for instructions. '''

import numpy as np
import run_globals

# Numerical solver options. These options worked well in the past.
# Artificial viscosity keeps spurious oscillations in check, and the positivity-
# preserving limiter keeps variables like density, pressure positive.
# SolutionOrder specifies the order of polynomial to use as approximation basis.
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
 'NumElemsX': 500,
 'xmax': 0,   # Top of the conduit is placed at -150 m for 2D coupling, but you can choose whatever if not coupling to 2D.
 'xmin': -5000.0
}

# Specify physical properties here.
# LaxFriedrichs is the numerical flux. This is a standard option, and does not
# typically need to be changed.

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

SourceTerms = {'source1': {'Function': 'GravitySource', # Gravity
             'gravity': 9.8,
             'source_treatment': 'Explicit'},
 'source2': {'Function': 'FrictionVolFracVariableMu',   # Friction source term for given conduit radius
             'conduit_radius': 5.,
             'viscosity_factor': 1/5,#1/20,
             'source_treatment': 'Implicit'
             },
 'source3': {'Function': 'ExsolutionSource',            # Exsolution source
             'source_treatment': 'Explicit',
             'tau_d': 0.001,                             # Exsolution timescale
             },
 "source4": {
        "Function": "FragmentationStrainRateSource",
        "tau_f": 0.003, # This is the fragmentation timescale
        "G": 1e9,
        "k": 0.001,
        "fragsmooth_scale": 0.010, # This is the smoothing scale that was also in the volume fraction source
        "which_criterion": "both", # This is either "shear" for only using shear strain to determine whether the magma fragments, "tensile" for only using the longitudinal (du/dz) strain rate, or "both" for using the larger of the two strain rates
        "conduit_radius": 5, # Conduit radius (make sure this is the same as what you provide to the Friction source term!
        'source_treatment': 'Explicit',
    },
}

Output = {'AutoPostProcess': False,
  'Prefix' : f"{run_globals.output_file_prefix}_cond",              # Output filename
  'WriteInitialSolution': True,
  'WriteInterval': 100,                                   # Output frequency (this many timesteps pass before file is written)
}

# Define a traction function
def gaussian_traction(x:np.array, total_pressure_change=10e6, x0=-350, sigma=50) -> np.array:
  ''' Traction function added to the hydrostatic equation. Negative
   sign indicates downward traction on the fluid. Units are Pa / m. 
   Total pressure change due to traction is amp * sigma.
   Inputs:
     x: array of points at which traction is evaluated
     total_pressure_change: total pressure change across the traction. The
       Gaussian amplitude is calculated from this.
     x0: Gaussian center (m)
     sigma: standard deviation parameter (length scale of traction function)
   '''
  # Compute amplitude of Gaussian TODO:
  amplitude = total_pressure_change / (np.sqrt(np.pi) * sigma)
  _t = (x-x0)/sigma
  return -amplitude * np.exp(-_t * _t)

radio = 5
f_plug = 1.79e8
len_plug = 60
t_plug = f_plug / (2*np.pi*radio*len_plug)
trac_par = 2*t_plug/radio

# Set common parameters
p_chamber = 64374916.25# 52164806.11
T_chamber = 950 + 273.15 # 1223.15
yC = 0.4    # Crystal mass fraction
yWt = 0.02  # Total water mass fraction

chi_water = (1.0 - yC) * yWt / (1 - yWt)

a = 1/7

# Initial condition parameters
# Note that some parameters here are repeated, and must be consistent with the
# parameters specified in SourceTerms. It's hardcoded here because this
# file was generated by a script that writes parallelized input files.
InitialCondition = {'Function': 'StaticPlug',          # Specify to call physics/multiphasevpT/functions > SteadyState as initial condition
 'p_chamber': p_chamber,
  # Magma properties
 'K_magma': 10000000000.0,
 'c_v_magma': 1000.0,
 'neglect_edfm': True,
 'p0_magma': 35999999.99999999,
 'rho0_magma': 2600.0,
 # Solubility properties
 'solubility_k': Physics["Solubility"]["k"],
 'solubility_n': Physics["Solubility"]["n"],
 'x_global': np.linspace(Mesh["xmin"], Mesh["xmax"], Mesh["NumElemsX"]), # All x meshes linked together
 # Traction due to plug; gaussian_traction is defined above
 'traction_fn': lambda x:  -trac_par * (1 / (1 + np.exp(-a * (x + len_plug)))), # -3183.09 # -7283.09 used for producing fragmentation
 # Distribution of total water, crystal, and temperature in conduit
 'yWt_fn': lambda x: yWt  - (0.01 / (1 + np.exp(-a * (x + len_plug)))) , # Value times constant 1 (with array shape like x)
 'yC_fn': lambda x: yC  + (0.5 / (1 + np.exp(-a * (x + len_plug)))),   # Value times constant 1 (with array shape like x)
 'T_fn': lambda x: T_chamber - (20 / (1 + np.exp(-a * (x + len_plug)))),   # Value times constant 1 (with array shape like x)
 'yF_fn': lambda x: 1/(1 + np.exp(-a * (x + len_plug))) , 
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
  # For running serial, using p boundary condition:
#  'x2': {'BCType': 'PressureOutlet1D',                   # Pressure outlet boundary condition (automatically chokes if needed)
        # 'p': 100000.0,                  # Boundary pressure (if flow not choked) -- with scale height factor
        # },
   'x2': {'BCType': 'MultiphasevpT2D1D',
          'bkey': 'vent'}
}

# Linked parallel solvers. If running in serial, leave as empty list.
LinkedSolvers = [
  {
    "DeckName": "vent_region.py",
    "BoundaryName": "vent", # set equal to the "bkey" param of a MultiphasevpT2D1D BoundaryCondition
  },
]

# Set timestepping options. The timestep size (dt) is calculated based on final
# time and NumTimeSteps. If a NonPhysicalError is returned, check here first to
# see if the CFL condition is met.6
TimeStepping = {'FinalTime': 0.001, # Final 
 'InitialTime': 0.0,
 'NumTimeSteps': 1000, # 2500000,# Number of timesteps to run for
 'TimeStepper': 'Strang', # 'FE', # 'RK3SR',  # 4-step RK3 scheme that maximizes CFL stability region per function eval
}

# dx = 5
# max wave speed w <~ 2000
# 1/(2p+1) = 1/5
# With CFL efficiency of 1:
# dt ~ dx/w * (1/5) ~ 1/2000 ~ 120/240000
# (RK3SR efficiency is > 1)