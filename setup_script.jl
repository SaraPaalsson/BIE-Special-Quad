include("initialization.jl")

# Setup starfish test domain
	Npanels = 10
	geomShape = "starfish"
	filename = "starfish"
	initialization.setup(Npanels,geomShape,filename)
