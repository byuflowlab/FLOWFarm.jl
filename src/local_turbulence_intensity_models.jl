abstract type AbstractLocalTurbulenceIntensityModel end

struct LocalTIModelNoLocalTI{} <: AbstractLocalTurbulenceIntensityModel

end

struct LocalTIModelMaxTI{TF} <: AbstractLocalTurbulenceIntensityModel
    astar::TF
    bstar::TF
end
LocalTIModelMaxTI() = LocalTIModelMaxTI(2.32, 0.154)


function calculate_local_ti(turbine_x, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_inflow_velcities, ct_model, ti_model::LocalTIModelNoLocalTI; turbine_id=1, tol=1E-6, k1=0.2, k2=0.003)
    return ambient_ti
end

# compute wake spread parameter based on local turbulence intensity
function _k_star_func(ti_ust;k1=0.2,k2=0.003)

    # calculate wake spread parameter from Niayifar and Porte Agel (2015, 2016)
    k_star_ust = k1*ti_ust + k2

    return k_star_ust

end

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
    wake_overlap = overlap_area_func(0.0, h_dst, d_dst, delta_y, h_ust, d_w)

    # initialize the wake/rotor area overlap ratio
    ti_area_ratio = 0.0

    # only include turbines with area overlap in the softmax
    if wake_overlap > 0.0

        # Calculate the turbulence added to the inflow of the downstream turbine by the
        # wake of the upstream turbine
        ti_added = 0.73*(axial_induction_ust^0.8325)*(ti_ust^0.0325)*((x/d_ust)^(-0.32))

        rotor_area_dst = 0.25*pi*d_dst^2.0
        ti_area_ratio_tmp = ti_added*(wake_overlap/rotor_area_dst)

        # Run through the smooth max to get an approximation of the true max TI area ratio
        ti_area_ratio = smooth_max(ti_area_ratio_in, ti_area_ratio_tmp, s=s)

        # Calculate the total turbulence intensity at the downstream turbine based on
        # the result of the smooth max function
        ti_dst = sqrt(ti_amb^2.0 + ti_area_ratio^2.0)

    end

    return ti_dst, ti_area_ratio

end


function calculate_local_ti(turbine_x, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_inflow_velcities, ct_model, ti_model::LocalTIModelMaxTI; turbine_id=1, tol=1E-6, k1=0.2, k2=0.003)

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

            # calculate ct at the current upstream turbine
            ct_model_ind = ct_model[turb]
            ct_ust = calculate_ct(turbine_inflow_velcities[turb], ct_model_ind)

            # determine the far-wake onset location
            astar = ti_model.astar
            bstar = ti_model.bstar
            x0 = _gauss_yaw_potential_core(d_ust, yaw_ust, ct_ust, astar, ti_ust, bstar)

            # calculate the distance from the onset of far-wake
            deltax0 = x - x0

            # calculate wake spread rate for current upstream turbine
            kstar_ust = _k_star_func(ti_ust,k1=k1,k2=k2)

            # calculate horizontal and vertical spread standard deviations
            sigmay = sigmaz = _gauss_yaw_spread(x, x0, kstar_ust, d_ust, yaw_ust)

            # determine the initial wake angle at the onset of far wake
            theta0 = _bpa_theta_0(yaw_ust, ct_ust)

            # horizontal cross-wind wake displacement from hub
            wake_offset = _bpa_deflection(d_ust, ct_ust, yaw_ust, kstar_ust, kstar_ust, sigmay, sigmaz, theta0, x0)

            # cross wind distance from point location to upstream turbine wake center
            delta_y = windfarmstate.turbine_y[turbine_id]  - (windfarmstate.turbine_y[turb] + wake_offset)

            # save ti_area_ratio and ti_dst to new memory locations to avoid
            # aliasing during differentiation
            ti_area_ratio_tmp = deepcopy(ti_area_ratio)

            # update local turbulence intensity
            ti_dst, ti_area_ratio = _niayifar_added_ti_function(x, d_dst, d_ust, h_ust, h_dst, ct_ust, kstar_ust, delta_y, ambient_ti, ti_ust, ti_dst, ti_area_ratio_tmp)

        end


    end

    return ti_dst

end


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

            turbine_type = windfarm.turbine_definition_ids[turb]
            rotor_diameter = rotor_diameter[turb]
            hub_height = hub_height[turb]
            dy = loc[2] - turbine_y[turb]
            dz = loc[3] - hub_height
            r = sqrt(dy^2 + dz^2)
            ct = turbine_ct[turb]

            kstar = 0.11*ct^1.07*ambient_ti^0.2
            epsilon = 0.23*ct^-0.25*ambient_ti^0.17
            d = 2.3*ct^-1.2
            f = 0.7*ct^-3.2*ambient_ti^-0.45

            dist = 0.5
            if r/rotor_diameter <= dist
                k1 = cos(pi/2.0*(r/rotor_diameter-dist))^2
                k2 = cos(pi/2.0*(r/rotor_diameter+dist))^2
            else
                k1 = 1.0
                k2 = 0.0
            end

            sigma = kstar*dx + epsilon*rotor_diameter
            if dz >= 0.0
                delta = 0.0
            else
                delta = ambient_ti*sin(pi*dz/hub_height)^2
            end

            #2.5 for low TI 2.0 for high TI
            sigma = sigma/div_sigma

            new_ex = -2.0 #orig -2.0
            p1 = 1.0/(d + e*dx/rotor_diameter + f*(1.0+dx/rotor_diameter)^new_ex)
            p2 = k1*exp(-(r-rotor_diameter/2.0)^2/(2.0*sigma^2)) + k2*exp(-(r+rotor_diameter/2.0)^2/(2.0*sigma^2))
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
