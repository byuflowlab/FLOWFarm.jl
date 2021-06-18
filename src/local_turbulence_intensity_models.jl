abstract type AbstractLocalTurbulenceIntensityModel end

"""
    LocalTIModelNoLocalTI()

Don't calculate local turbulence intensity. Ambient TI will be used instead for all points

"""
struct LocalTIModelNoLocalTI{} <: AbstractLocalTurbulenceIntensityModel

end


"""
    LocalTIModelMaxTI(astar, bstar, k1, k2)

Calculate local turbulence intensity using the model presented in Niayifar and 
Porte Agel (2015, 2016)

# Arguments
- `astar::Float`: wake spreading parameter from Bastankhah and Porte-Agel Gaussian wake model
- `bstar::Float`: wake spreading parameter from Bastankhah and Porte-Agel Gaussian wake model
- `k1::Float`: slope of k vs TI curve
- `k2::Float`: vertical offset of k vs TI curve
"""
struct LocalTIModelMaxTI{TF} <: AbstractLocalTurbulenceIntensityModel
    astar::TF
    bstar::TF
    k1::TF
    k2::TF
end
LocalTIModelMaxTI(x, y) = LocalTIModelMaxTI(x, y, 0.3837, 0.003678)
LocalTIModelMaxTI() = LocalTIModelMaxTI(2.32, 0.154, 0.3837, 0.003678)

"""
    calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
    turbine_inflow_velcities, turbine_ct, ti_model::LocalTIModelNoLocalTI; turbine_id=1, tol=1E-6)

Returns ambient turbulence intesity value whenever local turbulence intensity is requestesd

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine wind direction locations in the wind direction 
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine cross wind locations in the wind direction 
    reference frame
- `ambient_ti::Float`: ambient turbulence intensity
- `rotor_diameter::Array{Float,nTurbines}`: rotor diameters of all turbines
- `hub_height::Array{Float,nTurbines}`: hub heights of all turbines relative to the ground
- `turbine_yaw::Array{Float,nTurbines}`: yaw of all turbines for the current wind state in radians
- `turbine_local_ti::Array{Float,nTurbines}`: local turbulence intensity of all turbines for the current wind state`
- `sorted_turbine_index::Array{Float,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `turbine_inflow_velcities::Array{Float,nTurbines}`: effective inflow wind speed at each turbine for given state
- `turbine_ct::Array{Float,nTurbines}`: thrust coefficient of each turbine for the given state
- `ti_model::LocalTIModelNoLocalTI`: contains a struct defining the desired turbulence intensity model, no local TI in this case
- `turbine_id::Int`: index of wind turbine of interest. Provide 1 as default.
- `tol::Float`: How far upstream a turbine should be before being included in TI calculations
"""
function calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_inflow_velcities, turbine_ct, ti_model::LocalTIModelNoLocalTI; turbine_id=1, tol=1E-6)
    return ambient_ti
end

"""
    _k_star_func(ti_ust,k1, k2)

Calculate local wake spreading rate based on turbulence intensity using the model presented 
in Niayifar and Porte Agel (2015, 2016)

# Arguments
- `ti_ust::Float`: upstream turbine local turbulence intensity
- `k1::Float`: slope of k vs TI curve
- `k2::Float`: vertical offset of k vs TI curve
"""
# compute wake spread parameter based on local turbulence intensity
function _k_star_func(ti_ust,k1,k2)
    # calculate wake spread parameter from Niayifar and Porte Agel (2015, 2016)
 
    k_star_ust = k1*ti_ust + k2

    return k_star_ust

end
_k_star_func(x) = _k_star_func(x, 0.3837, 0.003678)

"""
    _niayifar_added_ti_function(x, d_dst, d_ust, h_ust, h_dst, ct_ust, kstar_ust, delta_y, 
        ti_amb, ti_ust, ti_dst, ti_area_ratio_in; s=700.0)

Main code for calculating the local turbulence intensity at a turbine using the method of
    Niayifar and Porte Agel (2015, 2016).

# Arguments
- `x::Float`: downstream distance from turbine to point of interest
- `d_dst::Float`: downstream turbine rotor diameter
- `d_ust::Float`: upstream turbine rotor diameter
- `h_ust::Float`: upstream turbine hub height
- `h_dst::Float`: downstream turbine hub height
- `ct_ust::Float`: upstream turbine thrust coefficient
- `kstar_ust::Float`: upstream turbine wake expansion rate
- `delta_y::Float`: cross wind separation from turbine to point of interest
- `ti_amb::Float`: ambient turbulence intensity
- `ti_ust::Float`: upstream turbine local turbulence intensity
- `ti_dst::Float`: downstream turbine local turbulence intensity
- `ti_area_ratio_in::Float`: current value of TI-area ratio for use in calculatin local TI
- `s::Float`: smooth max smootheness parameter
"""
function _niayifar_added_ti_function(x, d_dst, d_ust, h_ust, h_dst, ct_ust, kstar_ust, delta_y, ti_amb, ti_ust, ti_dst, ti_area_ratio_in; s=700.0)
    # Niayifar and Porte Agel 2015, 2016 using smooth max on area TI ratio

    # calculate axial induction based on the Ct value
    axial_induction_ust = _ct_to_axial_ind_func(ct_ust)

    # calculate BPA spread parameters Bastankhah and Porte Agel 2014
    beta = 0.5*((1.0 + sqrt(1.0 - ct_ust))/sqrt(1.0 - ct_ust))
    epsilon = 0.2*sqrt(beta)

    # calculate wake spread for TI calcs
    sigma = kstar_ust*x + d_ust*epsilon
    d_w = 4.0*sigma

    # calculate wake overlap ratio
    # println("wol check: ", h_dst, " ", d_dst, " ", delta_y, " ", h_ust, " ", d_w, " ", sigma)

    wake_overlap = overlap_area_func(0.0, h_dst, d_dst, delta_y, h_ust, d_w)

    # initialize the wake/rotor area overlap ratio
    ti_area_ratio = 0.0

    # only include turbines with area overlap in the softmax
    if wake_overlap > 0.0
        # Calculate the turbulence added to the inflow of the downstream turbine by the
        # wake of the upstream turbine
        ti_added = 0.73*(axial_induction_ust^0.8325)*(ti_ust^0.0325)*((x/d_ust)^(-0.32))

        rotor_area_dst = 0.25*pi*d_dst^2
        ti_area_ratio_tmp = ti_added*(wake_overlap/rotor_area_dst)

        # println("ti check: ", axial_induction_ust, " ", ti_ust, " ", d_ust, " ", rotor_area_dst, " ", ti_added, " ", wake_overlap)

        # Run through the smooth max to get an approximation of the true max TI area ratio
        ti_area_ratio = smooth_max(ti_area_ratio_in, ti_area_ratio_tmp, s=s)

        # Calculate the total turbulence intensity at the downstream turbine based on
        # the result of the smooth max function
        # println("ti check: ", ti_amb, " ", ti_area_ratio)
        ti_dst = norm([ti_amb, ti_area_ratio])

    end

    return ti_dst, ti_area_ratio

end

"""
    calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
    turbine_inflow_velcities, turbine_ct, ti_model::LocalTIModelMaxTI; turbine_id=1, tol=1E-6)

Returns local turbulence intensity calculated using Niayifar and Porte Agel 2015, 2016 using smooth max on area TI ratio

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine wind direction locations in the wind direction 
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine cross wind locations in the wind direction 
    reference frame
- `ambient_ti::Float`: ambient turbulence intensity
- `rotor_diameter::Array{Float,nTurbines}`: rotor diameters of all turbines
- `hub_height::Array{Float,nTurbines}`: hub heights of all turbines relative to the ground
- `turbine_yaw::Array{Float,nTurbines}`: yaw of all turbines for the current wind state in radians
- `turbine_local_ti::Array{Float,nTurbines}`: local turbulence intensity of all turbines for the current wind state`
- `sorted_turbine_index::Array{Float,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `turbine_inflow_velcities::Array{Float,nTurbines}`: effective inflow wind speed at each turbine for given state
- `turbine_ct::Array{Float,nTurbines}`: thrust coefficient of each turbine for the given state
- `ti_model::LocalTIModelMaxTI`: contains a struct defining the desired turbulence intensity model, no local TI in this case
- `turbine_id::Int`: index of wind turbine of interest. Provide 1 as default.
- `tol::Float`: How far upstream a turbine should be before being included in TI calculations
"""
function calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_inflow_velcities, turbine_ct, ti_model::LocalTIModelMaxTI; turbine_id=1, tol=1E-6)

    # calculate local turbulence intensity at turbI

    # initialize the ri_dst and ti_area_ratio to 0.0 for current turbine
    ti_area_ratio = 0.0
    ti_dst = copy(ambient_ti)

    # extract downstream turbine information
    d_dst = rotor_diameter[turbine_id]
    h_dst = hub_height[turbine_id]

    # extract the number of turbines
    nturbines = length(rotor_diameter)

    # loop over upstream turbines
    for u=1:nturbines

        # get index of upstream turbine
        turb = sorted_turbine_index[u]

        # skip turbine's influence on itself
        if turb == turbine_id; continue; end

        # calculate downstream distance between wind turbines
        x = turbine_x[turbine_id] - turbine_x[turb]

        if x > tol

            # extract state and design info for current upstream turbine
            d_ust = rotor_diameter[turb]
            h_ust = hub_height[turb]
            yaw_ust = turbine_yaw[turb]
            ti_ust = turbine_local_ti[turb]
            ct_ust = turbine_ct[turb]

            # determine the far-wake onset location
            astar = ti_model.astar
            bstar = ti_model.bstar
            x0 = _gauss_yaw_potential_core(d_ust, yaw_ust, ct_ust, astar, ti_ust, bstar)

            # calculate wake spread rate for current upstream turbine
            kstar_ust = _k_star_func(ti_ust,ti_model.k1,ti_model.k2)

            # calculate horizontal and vertical spread standard deviations
            # println("sigma inputs: ", x, " ", x0, " ", kstar_ust, " ", d_ust, " ", yaw_ust)
            sigmay = sigmaz = _gauss_yaw_spread(d_ust, kstar_ust, x, x0, yaw_ust)

            # determine the initial wake angle at the onset of far wake
            theta0 = _bpa_theta_0(yaw_ust, ct_ust)

            # horizontal cross-wind wake displacement from hub
            println("wake offset: ", d_ust, " ", ct_ust, " ", yaw_ust, " ", kstar_ust, " ", kstar_ust, " ", sigmay, " ", sigmaz, " ", theta0, " ", x0)
            wake_offset = _bpa_deflection(d_ust, ct_ust, yaw_ust, kstar_ust, kstar_ust, sigmay, sigmaz, theta0, x0)

            # cross wind distance from point location to upstream turbine wake center
            delta_y = turbine_y[turbine_id]  - (turbine_y[turb] + wake_offset)

            # save ti_area_ratio and ti_dst to new memory locations to avoid
            # aliasing during differentiation
            ti_area_ratio_tmp = deepcopy(ti_area_ratio)

            # update local turbulence intensity
            ti_dst, ti_area_ratio = _niayifar_added_ti_function(x, d_dst, d_ust, h_ust, h_dst, ct_ust, kstar_ust, delta_y, ambient_ti, ti_ust, ti_dst, ti_area_ratio_tmp)

        end

    end

    if turbine_id == 10
        println("output")
        println(ambient_ti, " ", turbine_ct[turbine_id], " ", turbine_x[turbine_id], " ", rotor_diameter[turbine_id], " ", hub_height[turbine_id], " ",700, " ", ti_dst)
    end

    return ti_dst

end

"""
    GaussianTI(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, 
        sorted_turbine_index, ambient_ti; div_sigma=2.5, div_ti=1.2)


Calculate local turbulence intensity based on "On wake modeling, wind-farm gradients and AEP 
    predictions at the Anholt wind farm" by Pena Diaz, Alfredo; Hansen, Kurt Schaldemose; 
    Ott, Søren; van der Laan, Paul ??

# Arguments
- `loc::Array{Float,3}`: [x,y,z] location of point of interest in wind direction ref. frame
- `turbine_x::Array{Float,nTurbines}`: turbine wind direction locations in the wind direction 
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine cross wind locations in the wind direction 
    reference frame
- `rotor_diameter::Array{Float,nTurbines}`: rotor diameters of all turbines
- `hub_height::Array{Float,nTurbines}`: hub heights of all turbines relative to the ground
- `turbine_ct::Array{Float,nTurbines}`: thrust coefficient of each turbine for the given state
- `sorted_turbine_index::Array{Float,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `ambient_ti::Float`: ambient turbulence intensity
- `div_sigma::Float`: ?
- `div_ti::Float`: ?
"""
function GaussianTI(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_ti; div_sigma=2.5, div_ti=1.2)

    added_ti = 0.0
    e = 1.0*ambient_ti^0.1
    nturbines = length(turbine_x)

    for u=1:nturbines

        # get index of upstream turbine
        turb = sorted_turbine_index[u]

        # calculate downstream distance between wind turbines
        dx = loc[1] - turbine_x[turb]

        if dx > 1e-6
            dy = loc[2] - turbine_y[turb]
            dz = loc[3] - hub_height[turb]
            r = norm([dy, dz])
            ct = turbine_ct[turb]

            kstar = 0.11*ct^1.07*ambient_ti^0.2
            epsilon = 0.23*ct^-0.25*ambient_ti^0.17
            d = 2.3*ct^-1.2
            f = 0.7*ct^-3.2*ambient_ti^-0.45

            dist = 0.5
            if r/rotor_diameter[turb] <= dist
                k1 = cos(pi/2.0*(r/rotor_diameter[turb]-dist))^2
                k2 = cos(pi/2.0*(r/rotor_diameter[turb]+dist))^2
            else
                k1 = 1.0
                k2 = 0.0
            end

            sigma = kstar*dx + epsilon*rotor_diameter[turb]
            if dz >= 0.0
                delta = 0.0
            else
                delta = ambient_ti*sin(pi*dz/hub_height[turb])^2
            end

            #2.5 for low TI 2.0 for high TI
            sigma = sigma/div_sigma

            new_ex = -2.0 #orig -2.0
            p1 = 1.0/(d + e*dx/rotor_diameter[turb] + f*(1.0+dx/rotor_diameter[turb])^new_ex)
            p2 = k1*exp(-(r-rotor_diameter[turb]/2.0)^2/(2.0*sigma^2)) + k2*exp(-(r+rotor_diameter[turb]/2.0)^2/(2.0*sigma^2))
            dI = p1*p2 - delta
            #1.2 for low TI 2.0 for high TI
            added_ti += dI/div_ti
        end
    end
    return ambient_ti + added_ti
end


# function GaussianTI_stanley(loc,windfarm,windfarmstate,ambient_ti)
#
#     added_ti = 0.0
#     e = 1.0*ambient_ti^0.1
#     nturbines = length(windfarm.turbine_x)
#
#     for u=1:nturbines
#
#         # get index of upstream turbine
#         turb = windfarmstate.sorted_turbine_index[u]
#
#         # calculate downstream distance between wind turbines
#         dx = loc[1] - windfarmstate.turbine_x[turb]
#
#         if dx > 1e-6
#             turbine_type = windfarm.turbine_definition_ids[turb]
#             rotor_diameter = windfarm.turbine_definitions[turbine_type].rotor_diameter[1]
#             hub_height = windfarm.turbine_definitions[turbine_type].hub_height[1]
#             dy = loc[2] - windfarmstate.turbine_y[turb]
#             dz = loc[3] - hub_height
#             r = sqrt(dy^2 + dz^2)
#
#             ct = windfarmstate.turbine_ct[turb]
#             kstar = 0.11*ct^1.07*ambient_ti^0.2
#             epsilon = 0.23*ct^-0.25*ambient_ti^0.17
#             d = 2.3*ct^-1.2
#             f = 0.7*ct^-3.2*ambient_ti^-0.45
#
#             # dist = 0.0
#             # if r/rotor_diameter <= dist
#             #     k1 = cos(pi/2.0*(r/rotor_diameter-dist))^2
#             #     k2 = cos(pi/2.0*(r/rotor_diameter+dist))^2
#             # else
#             #     k1 = 1.0
#             #     k2 = 0.0
#             # end
#             k1 = 1.0
#             k2 = 0.0
#
#             sigma = kstar*dx + epsilon*rotor_diameter
#             if dz >= 0.0
#                 delta = 0.0
#             else
#                 delta = ambient_ti*sin(pi*dz/hub_height)^2
#             end
#
#
#             sigma = sigma*0.5
#             p1 = 1.0/(d + e*dx/rotor_diameter + f*(1.0+dx/rotor_diameter)^-2.0)
#             # p2 = k1*exp(-(r-rotor_diameter/2.0)^2/(2.0*sigma^2)) + k2*exp(-(r+rotor_diameter/2.0)^2/(2.0*sigma^2))
#             p2 = exp(-(r/2.0)^2/(2.0*sigma^2))
#             dI = p1*p2 - delta
#             # if r < rotor_diameter*4.0/5.0
#             added_ti += dI
#             # end
#         end
#     end
#     return ambient_ti + added_ti
# end
