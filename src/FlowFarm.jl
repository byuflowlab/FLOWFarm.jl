module FlowFarm

# using CCBlade
# using PyPlot
using FLOWMath: linear
# using Statistics

include("turbines.jl")
include("windfarms.jl")
include("wind_resources.jl")
include("wind_shear_models.jl")
include("wake_combination_models.jl")
include("wake_deficit_models.jl")
include("wake_deflection_models.jl")
include("thrust_coefficient_models.jl")
include("local_turbulence_intensity_models.jl")
include("general_models.jl")
include("power_models.jl")
include("user_functions.jl")
include("optimization_functions.jl")
# include("fatigue_functions.jl")

end # module
