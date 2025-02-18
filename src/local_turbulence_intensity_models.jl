export LocalTIModelNoLocalTI, LocalTIModelMaxTI, LocalTIModelGaussTI
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
struct LocalTIModelMaxTI{T1,T2,T3,T4} <: AbstractLocalTurbulenceIntensityModel
    astar::T1
    bstar::T2
    k1::T3
    k2::T4
end
LocalTIModelMaxTI(x, y) = LocalTIModelMaxTI(x, y, 0.3837, 0.003678)
LocalTIModelMaxTI() = LocalTIModelMaxTI(2.32, 0.154, 0.3837, 0.003678)

"""
    LocalTIModelGaussTI()

Calculate local turbulence intensity using the model presented in Qian and
Ishihara (2018)

"""
struct LocalTIModelGaussTI{} <: AbstractLocalTurbulenceIntensityModel

end

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
    ti_area_ratio = ti_area_ratio_in

    # only include turbines with area overlap in the softmax
    if wake_overlap > 0.0
        # Calculate the turbulence added to the inflow of the downstream turbine by the
        # wake of the upstream turbine based on Crespo, A.; Hernandez, J. Turbulence
        # characteristics in wind-turbine wakes. J. Wind Eng. Ind. Aerodyn. 1996,
        # 61, 71–85.
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
                    turbine_inflow_velocities, turbine_ct, ti_model::LocalTIModelMaxTI; turbine_id=1, tol=1E-6)

    # calculate local turbulence intensity at turbI

    # initialize the ti_dst and ti_area_ratio to 0.0 for current turbine
    ti_area_ratio = 0.0
    ti_dst = copy(ambient_ti)

    # extract downstream turbine information
    d_dst = rotor_diameter[turbine_id]
    h_dst = hub_height[turbine_id]

    # extract the number of turbines
    nturbines = length(turbine_x)

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

            # calculate the discontinuity point of the gauss yaw model
            xd = _gauss_yaw_discontinuity(d_ust, x0, kstar_ust, kstar_ust, yaw_ust, ct_ust)

            # calculate horizontal wake spread
            sigmay = _gauss_yaw_spread_interpolated(d_ust, kstar_ust, x, x0, yaw_ust, xd)

            # calculate vertical wake spread
            sigmaz = _gauss_yaw_spread_interpolated(d_ust, kstar_ust, x, x0, 0.0, xd)

            # determine the initial wake angle at the onset of far wake
            theta0 = _bpa_theta_0(yaw_ust, ct_ust)

            # horizontal cross-wind wake displacement from hub
            # println("wake offset: ", ct_ust, " ", kstar_ust, " ", kstar_ust, " ", sigmay, " ", sigmaz, " ", theta0, " ", x0)

            # println("wake offset: ", ct, " ", ky, " ", kz, " ", sigmay, " ", sigmaz, " ", theta0, " ", x0)
            wake_offset = _bpa_deflection(d_ust, ct_ust, yaw_ust, kstar_ust, kstar_ust, sigmay, sigmaz, theta0, x0)

            # cross wind distance from point location to upstream turbine wake center
            delta_y = turbine_y[turbine_id]  - (turbine_y[turb] + wake_offset)

            # save ti_area_ratio and ti_dst to new memory locations to avoid
            # aliasing during differentiation
            ti_area_ratio_tmp = deepcopy(ti_area_ratio)

            # update local turbulence intensity
            ti_dst, ti_area_ratio = _niayifar_added_ti_function(x, d_dst, d_ust, h_ust, h_dst, ct_ust, kstar_ust, delta_y, ambient_ti, ti_ust, ti_dst, ti_area_ratio_tmp)
            # println("ti output: ", turbine_id, " ", turb, " ", ti_area_ratio, " ", ti_dst)
        end

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

"""
    calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
    turbine_inflow_velcities, turbine_ct, ti_model::LocalTIModelGaussTI; turbine_id=1, tol=1E-6)

Returns local turbulence intensity calculated using methods in Qian 2018 from the Journal of Wind Energy https://doi.org/10.1016/j.jweia.2018.04.010
with modification to account for yaw coming from Qian 2018 from Energies doi:10.3390/en11030665

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
- `ti_model::LocalTIModelGaussTI`: contains a struct defining the desired turbulence intensity model
- `turbine_id::Int`: index of wind turbine of interest. Provide 1 as default.
- `tol::Float`: How far upstream a turbine should be before being included in TI calculations
"""
function calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_inflow_velocities, turbine_ct, ti_model::LocalTIModelGaussTI; turbine_id=1, tol=1E-6)

    nturbines = length(turbine_x)
    x_loc = turbine_x[turbine_id]
    y_loc = turbine_y[turbine_id]
    z_loc = hub_height[turbine_id]

    intensity = 0
    Ia = ambient_ti

    for i = 1:nturbines
        # get index of upstream turbine
        upstream_turbine = sorted_turbine_index[i]

        # put yaw in counter-clockwise notation as descirbed in Qian 2018 from Energies
        yaw = turbine_yaw[upstream_turbine]
        yaw *= -1

        # get modified ct from Qian 2018 Energies
        ct = turbine_ct[upstream_turbine] * cos(yaw)^3

        # get downstream distance between turbines
        dx = x_loc - turbine_x[upstream_turbine]

        if dx > tol
            dy = y_loc - turbine_y[upstream_turbine]
            dz = z_loc - hub_height[upstream_turbine]

            # Modifications for yaw
            theta_0 = 0.3 * yaw/cos(yaw) * (1 - sqrt(1 - ct))
            yd = theta_0 * dx
            r = sqrt(dz^2 + (dy + yd)^2)

            # calculate constants
            k_star = 0.11 * ct^1.07 * Ia^0.20
            epsilon = 0.23 * ct^(-0.25) * Ia^0.17
            d = 2.3 * ct^(-1.2)
            e = Ia^0.1
            f = 0.7 * ct^(-3.2) * Ia^(-0.45)
            k1 = 1
            k2 = 0
            D = rotor_diameter[upstream_turbine]
            if r/D <= 0.5
                k1 = (cos(pi/2 * (r/D - 0.5)))^2
                k2 = (cos(pi/2 * (r/D + 0.5)))^2
            end
            sigma = k_star * dx + epsilon * D

            if dz >= 0 || hub_height[upstream_turbine] == 0 #to avoid divide by zero
                delta = 0
            else
                delta = Ia * sin(pi*(-dz/hub_height[upstream_turbine]))^2
            end
            intensity += ((1 / (d + e*dx/D + f*(1+dx/D)^-2)) * (k1*exp(-((r-D/2)^2/(2*sigma^2))) + k2*exp(-((r+D/2)^2/(2*sigma^2)))) - delta)^10
        else
            break
        end
    end

    turbine_local_ti[turbine_id] = intensity^(1/10) + ambient_ti

    return turbine_local_ti[turbine_id]
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
