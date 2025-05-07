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
	"File" : f"{BASE_PATH}/scenarios/meshes/{run_globals.mesh_prefix}2.msh",
}

Output = {
	"Prefix" : f"{run_globals.output_file_prefix}_atm2",
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
	"ground2" : {
		"BCType" : "SlipWall",
	},
	"symmetry2" : {
		"BCType" : "SlipWall",
	},
	"r2" : {
		"BCType" : "MultiphasevpT2D2D",
		"bkey": "r2",
	},
	"r1" : {
		"BCType" : "MultiphasevpT2D2D",
		"bkey": "r1",
	},
}

LinkedSolvers = [
	{
		"DeckName": "r2r3.py",
		"BoundaryName": "r2",
	},
]

extend_atm = True

if not extend_atm:
	BoundaryConditions["r2"] = {
		"BCType" : "LinearizedImpedance2D",
	}
	LinkedSolvers = []