using Snopt
using DelimitedFiles 
using PyPlot

# import model set with wind farm and related details
include("./model_sets/model_set_10_amalia_wind_park.jl")

AEP = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)

# add turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

println("AEP: ", AEP)

# set up and show plot
axis("square")
plt.xlim([0, 5000])
plt.ylim([0, 5000])
plt.show()
