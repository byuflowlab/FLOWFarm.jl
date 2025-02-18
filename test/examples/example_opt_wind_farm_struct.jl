using SNOW
using DelimitedFiles
import FLOWFarm; const ff = FLOWFarm
using SparseArrays

function opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)

    # calculate spacing constraint value and jacobian
    ff.calculate_spacing_jacobian!(spacing_struct,x)

    # calculate boundary constraint and jacobian
    ff.calculate_boundary_jacobian!(boundary_struct,x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_struct.spacing_vec; boundary_struct.boundary_vec]
    g[:] .= ((c[:]))
    dg[:] = (([spacing_struct.jacobian; boundary_struct.jacobian]))

    # calculate AEP and gradient
    ff.calculate_aep_gradient!(farm,x)

    AEP = -farm.AEP[1]
    df[:] .= ((-farm.AEP_gradient))

    return AEP
end

function update_function(farm,x,r)
    n = length(farm.turbine_x)
    x .= ((x .* 2.0) .- 1.0) .* r
    farm.turbine_x .= x[1:n]
    farm.turbine_y .= x[n+1:end]
    x .= ((x ./ r) .+ 1.0) ./ 2.0
    return nothing
end

# import model set with wind farm and related details
include("../model_sets/model_set_9_38turb_round_farm.jl")

# set wind farm boundary parameters
boundary_center = [0.0,0.0]
boundary_radius = 1225.8227848101264

x = ((turbine_x ./ boundary_radius) .+ 1.0) ./ 2
y = ((turbine_y ./ boundary_radius) .+ 1.0) ./ 2

# initialize design variable array
x0 = [deepcopy(x);deepcopy(y)]

num_turbines = n = length(turbine_x)

# set up boundary function
f(farm,x) = update_function(farm,x,boundary_radius)

# build wind farm struct
    # opt_x and opt_y are set to true to optimize turbine x and y positions
farm = ff.build_wind_farm_struct(x0,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
            rotor_diameter,ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,
            rated_power,windresource,power_models,model_set,f;input_type="ForwardDiff",opt_x=true,opt_y=true)

spacing = 2.0 * rotor_diameter[1] # 2 rotor diameters spacing between turbines
scaling = 1.0/rotor_diameter[1]^2 # scaling factor for spacing constraint
spacing_struct = ff.build_spacing_struct(x0,num_turbines,spacing,scaling,f)

boundary_function(a,x,y) = ff.circle_boundary!(a,boundary_center,boundary_radius,x,y) # boundary function
n_constraints = num_turbines # number of boundary constraints
scaling = 1.0/boundary_radius^3.5 # scaling factor for boundary constraint
boundary_struct = ff.build_boundary_struct(x0,num_turbines,n_constraints,scaling,boundary_function,f)


ip_options = Dict(
    "max_iter" => 30,
    "tol" => 1e-6
)
solver = IPOPT(ip_options)
options = Options(derivatives=UserDeriv();solver)

lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]
lx = zeros(length(x0))
ux = ones(length(x0))
ng = length(lg)

x_start = deepcopy(x0)

opt!(g,df,dg,x) = opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)
xopt, fopt, info = minimize(opt!, x0, ng, lx, ux, lg, ug, options)