# 2D axisymmetric part near the vent

# TimeStepping is set in input_file.py

import run_globals

BASE_PATH='/Users/paxton/git/quail_volcano'

Numerics = {
	"SolutionOrder" : 0,
	"SolutionBasis" : "LagrangeSeg",
	"Solver" : "DG",
	"ApplyLimiters" : "PositivityPreservingMultiphasevpT",
	"ArtificialViscosity" : True,
	"AVParameter" : 30,
	'L2InitialCondition': False, # False == Use interpolation instead of L2 projection of Riemann data
}

Mesh = {
    "File" : None,
    "ElementShape" : "Segment",
    "NumElemsX" : 1e4,
    "xmin" : 0,
    "xmax" : 1000,
}


Output = {
	"Prefix" : f"{run_globals.output_file_prefix}_atm1",
	"WriteInterval" : run_globals.write_interval_1D,
	"WriteInitialSolution" : True,
	"AutoPostProcess": False,
}

Physics = {
    "Type" : "MultiphasevpT",
    "ConvFluxNumerical" : "LaxFriedrichs",
}

# Source terms: exsolution, fragmentation are not turned on
SourceTerms = {
	#"source1": {
	#	"Function" : "GravitySource",
	#	"gravity": 9.8,
	#	'source_treatment': 'Explicit',
	#},
}

# Set initial condition specified in run_globals.py
InitialCondition = run_globals.InitialCondition

# Dummy function passed to quail
ExactSolution = InitialCondition.copy()

# Set boundary conditions (and provisionally add impedance BC at outer boundary)
BoundaryConditions = {
	"x1" : {
		"BCType" : "MultiphasevpT1D1D",
		"bkey": "vent",
	},
	'x2': {
        'BCType': 'PressureOutlet1D',   # Inlet boundary condition
        'p': 1e5,                      # Pressure
    },
}


LinkedSolvers = []