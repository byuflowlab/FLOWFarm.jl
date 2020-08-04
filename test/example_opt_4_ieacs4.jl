using FlowFarm; const ff = FlowFarm
using Snopt
using DelimitedFiles
using Distributed
#using PyPlot
import ForwardDiff
import YAML
using CSV
include("iea37_specific_functions.jl")

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.splined_boundary_discreet_regions(turbine_x, turbine_y, params.bndry_x_clsd, params.bndry_y_clsd, params.num_bndry_verts, params.bndry_corner_indcies, params.turbs_per_region)
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.0*params.rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
@everywhere function aep_wrapper(x, params)
    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, params.turbine_z, params.rotor_diameter,
                params.hub_height, params.turbine_yaw, params.ct_models, params.generator_efficiency, params.cut_in_speed,
                params.cut_out_speed, params.rated_speed, params.rated_power, params.windresource, params.power_models, params.model_set,
                rotor_sample_points_y=params.rotor_points_y,rotor_sample_points_z=params.rotor_points_z, hours_per_year=365.0*24.0)
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt(x)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)
    
    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_con; boundary_con]
    dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# import model set with wind farm and related details
#include("./model_sets/model_set_7_ieacs4.jl")
include("./model_sets/model_set_7_ieacs4_reduced_wind_rose.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-9

# set globals for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    bndry_x_clsd
    bndry_y_clsd
    num_bndry_verts
    bndry_corner_indcies
    turbs_per_region
    obj_scale
    hub_height
    turbine_yaw
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

#--- Read in windfarm boundary data ---#
# Which case study we're doing. 'cs3' or 'cs4'
str_case = "4"
#- Rip the boundary coordinates from the .yaml file -#
file_dir = "./inputfiles/"
bnry_file_name_orig = "iea37-boundary-cs" * str_case * ".yaml"
bnry_file_name = string(file_dir,bnry_file_name_orig)
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name) # Make a matrix of all regions and boundary points for each region
bndry_x_clsd, bndry_y_clsd = ff.closeBndryLists(bndry_x, bndry_y)

#--- Read in random turbine locations ---#
# Make an array of the number of turbines in each region
nNumRegions = 5     # Number of reigons we're using (cs4 = 5, cs3 = 1)
turbs_per_region = zeros(Int8, nNumRegions)  # Preallocated turbines in each region
num_bndry_verts = zeros(Int8, nNumRegions)
for cntr in 1:nNumRegions
    num_bndry_verts[cntr] = length(getCs34VertList(getCs34Name(cntr)))
    turbs_per_region[cntr] = floor(getCs34NumTurbs(getCs34Name(cntr)))
end
bndry_corner_indcies = getCs34VertList("All")
num_tot_turbs = sum(turbs_per_region)

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, bndry_x_clsd, bndry_y_clsd, num_bndry_verts, bndry_corner_indcies, turbs_per_region, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# # Pull the pre-made random starts
# which_prestart = 1
# x0l = [ Float64[] for i in 1:nNumRegions ]
# for i in 1:nNumRegions    # Loop through our regions
#     PreStartsXY = convert(Matrix{Float64}, CSV.read(string(file_dir, "iea37-randostarts-", getCs34Name(i), "-200.csv"), delim=','))
#     x0l[i] = append!(x0l[i], PreStartsXY[which_prestart,:])
# end

# # Convert to arrays of x and y coordinates
# turbine_x = zeros(num_tot_turbs)
# turbine_y = zeros(num_tot_turbs)
# cntr = 0
# prev_index = 1

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# report initial objective value
println("Number of turbines: ", num_tot_turbs)
println("Rotor diameter: ", rotor_diameter[1])
println("Starting AEP value (GWh): ", aep_wrapper(x, params)[1]*1e-9/obj_scale)

# println("Directional AEP at start: ", dir_aep.*1E-6)

# continue
# # add initial turbine location to plot
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true,color="black"))#"C0"))
# end

# set general lower and upper bounds for design variables
lb = zeros(length(x))
ub = zeros(length(x)) .+ Inf

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 1 #3
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary.out"
options["Print file"] = "print.out"

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)

# run and time optimization
println
t1 = time()
xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)
t2 = time()
clkt = t2-t1

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end objective value: ", aep_wrapper(xopt))

# extract final turbine locations
turbine_x = copy(xopt[1:num_tot_turbs])
turbine_y = copy(xopt[num_tot_turbs+1:end])

#-- Save our optimized locations --#
#- Make sure the file doesn't exit -#
directory = "./results/"
file_name = "turblocs-bpm"
file_type = "yaml"
save_filename = ff.getNextFileName(directory, file_name, file_type)

# Necessary variables for writing turb locations
t = "IEA Wind Task 37 case study 4, BYU's BPM/SNOPT optimized layout"
td = "baseline layout for the 25 turbine wind plant model for IEA Task 37 case study 4"
tf ="iea37-10mw.yaml"
lu ="m"
wmu ="iea37-aepcalc.jl"
wrf ="iea37-windrose-cs4.yaml"
#aepd = aep_wrapper(x, params)
aepd = aep_wrapper(xopt, params)
aept = sum(aepd)
aepu ="MWh"
by="./inputfiles/default.yaml"
# Actually write the file
ff.write_turb_loc_YAML(save_filename, turbine_x, turbine_y; title=t, titledescription=td, 
    turbinefile=tf, locunits=lu, wakemodelused=wmu, windresourcefile=wrf, aeptotal=t, 
    aepdirs=aepd, aepunits=aepu, baseyaml=by)


# # add final turbine locations to plot
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
# end

# add wind farm boundary to plot
# plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))

# # set up and show plot
# axis("square")
# xlim(-boundary_radius-200,boundary_radius+200)
# ylim(-boundary_radius-200,boundary_radius+200)
# plt.show()


# #---- Debug plotting ----#
# # Plot the Boundary
# for cntr in 1:nNumRegions
#     plot(bndry_x_clsd[cntr], bndry_y_clsd[cntr])
# end

# # Plot the turbines
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=true,color="black"))
#     plt.text(turbine_x[i]+rotor_diameter,turbine_y[i]+rotor_diameter, string(i))
# end

# axis("square")
# axis("off")
# plt.show()