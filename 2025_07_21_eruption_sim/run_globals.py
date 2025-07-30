# Global parameters that are relevant to multiple input files 

ElementOrder1D = 1
ElementOrder2D = 1
InitialCondition = {
   "Function" : "HomogeneousAtmosphere",
  "h0": -5000,
  "p_atm": 1e5,
}

# Output file prefix (no trailing underscore)
output_file_prefix = "test_infrasound_v1"
write_interval = 200

# Mesh file prefix (no trailing underscore)
# There are some meshes in the scenarios/meshes/ folder
# tungurahuaA1 to tungurahuaA8 are different parts of the same mesh
# A1 is the section closest to the vent, A8 is the section farthest
mesh_prefix = "tungurahuaA"