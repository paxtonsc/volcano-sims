# 2D axisymmetric part near the vent

# TimeStepping is set in input_file.py

import run_globals

BASE_PATH='/Users/paxton/git/quail_volcano'

Numerics = {
	"SolutionOrder" : run_globals.ElementOrder2D,
	"SolutionBasis" : "LagrangeTri",
	"Solver" : "DG",
	"ApplyLimiters" : "PositivityPreservingMultiphasevpT",
	"ArtificialViscosity" : True,
	"AVParameter" : 30,
	'L2InitialCondition': False, # False == Use interpolation instead of L2 projection of Riemann data
}

Mesh = {
	"File" : f"{BASE_PATH}/scenarios/meshes/{run_globals.mesh_prefix}3.msh",
}

Output = {
	"Prefix" : f"{run_globals.output_file_prefix}_atm3",
	"WriteInterval" : run_globals.write_interval_2D,
	"WriteInitialSolution" : True,
	"AutoPostProcess": False,
}

Physics = {
    "Type" : "MultiphasevpT",
    "ConvFluxNumerical" : "LaxFriedrichs",
}

# Source terms: exsolution, fragmentation are not turned on
SourceTerms = {
	"source1": {
		"Function" : "GravitySource",
		"gravity": 9.8,
		'source_treatment': 'Explicit',
	},
	"source2": {
		"Function" : "CylindricalGeometricSource",
		'source_treatment': 'Explicit',
	}
}

# Set initial condition specified in run_globals.py
InitialCondition = run_globals.InitialCondition

# Dummy function passed to quail
ExactSolution = InitialCondition.copy()

# Set boundary conditions (and provisionally add impedance BC at outer boundary)
BoundaryConditions = {
	"ground3" : {
		"BCType" : "SlipWall",
	},
	"symmetry3" : {
		"BCType" : "SlipWall",
	},
	"r3" : {
		"BCType" : "MultiphasevpT2D2D",
		"bkey": "r3",
	},
	"r2" : {
		"BCType" : "MultiphasevpT2D2D",
		"bkey": "r2",
	},
}

LinkedSolvers = [
	{
		"DeckName": "r3r4.py",
		"BoundaryName": "r3",
	},
]

extend_atm = False

if not extend_atm:
	BoundaryConditions["r3"] = {
		"BCType" : "LinearizedImpedance2D",
	}
	LinkedSolvers = []