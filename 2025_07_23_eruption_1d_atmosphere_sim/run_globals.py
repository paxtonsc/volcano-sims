# Global parameters that are relevant to multiple input files 

ElementOrder1D = 1
ElementOrder2D = 1
InitialCondition = {
	"Function" : "HomogeneousAtmosphere1D",
  "h0": 5000,
  "p_atm": 1e5,
}

# Output file prefix (no trailing underscore)
output_file_prefix = "vertical_atmosphere_001"
write_interval_1D = int(5e3)

