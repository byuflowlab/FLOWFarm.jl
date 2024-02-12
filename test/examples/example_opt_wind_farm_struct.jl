using SNOW
using DelimitedFiles
import FLOWFarm; const ff = FLOWFarm

function opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)

    # calculate spacing constraint value and jacobian
    ff.calculate_spacing!(spacing_struct.spacing_vec,x,farm,spacing_struct)
    ds(a,x) = ff.calculate_spacing!(a,x,farm,spacing_struct)
    ForwardDiff.jacobian!(spacing_struct.spacing_jacobian,ds,spacing_vec,x,spacing_struct.forward_cfg)

    # calculate boundary constraint and jacobian
    ff.calculate_boundary!(boundary_struct.boundary_vec,x,farm,boundary_struct)
    db(a,x) = ff.calculate_boundary!(a,x,farm,boundary_struct)
    ForwardDiff.jacobian!(boundary_struct.boundary_jacobian,db,boundary_vec,x,boundary_struct.forward_cfg)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_struct.spacing_vec; boundary_struct.boundary_vec]
    g[:] .= c[:]

    dcdx = [dropzeros(sparse(spacing_struct.spacing_jacobian)); dropzeros(sparse(boundary_struct.boundary_jacobian))]
    dg[:] .= dcdx[:]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -ff.calculate_aep(x,farm)

    aep(x) = ff.calculate_aep(x,farm)
    ForwardDiff.gradient!(farm.AEP_gradient,df,x,farm.duals.forward_cfg)

    df[:] .= -farm.AEP_gradient

    return AEP
end

# import model set with wind farm and related details
include("../model_sets/model_set_9_38turb_round_farm.jl")

# set wind farm boundary parameters
boundary_center = [0.0,0.0]
boundary_radius = 1225.8227848101264

# initialize design variable array
x0 = [copy(turbine_x);copy(turbine_y)]

function update_function(farm,x,n)
    ff.update_turbine_x!(farm,x[1:n])
    ff.update_turbine_y!(farm,x[n+1:end])
    return nothing
end
n = length(turbine_x)
f(farm,x) = update_function(farm,x,n)


farm = ff.build_wind_farm_struct(x0,turbine_x,turbine_y,turbine_z,hub_height,turbine_yaw,
            rotor_diameter,ct_models,generator_efficiency,cut_in_speed,cut_out_speed,rated_speed,
            rated_power,windresource,power_models,model_set;update_function=f,opt_x=true,opt_y=true)

spacing_struct = ff.build_spacing_struct(x,n,1,1)
boundary_function(a,x,y) = ff.circle_boundary!(a,boundary_centerboundary_radius,x,y)
boundary_struct = ff.build_boundary_struct(x,n,1,boundary_function)

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

opt!(g,df,dg,x) = opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)
xopt, fopt, info = minimize(opt!, x0, ng, lx, ux, lg, ug, options)

println("Done")
