export JensenTopHat, JensenCosine, MultiZone, GaussOriginal, GaussYaw, GaussYawVariableSpread, GaussSimple, CumulativeCurl, NoWakeDeficit
abstract type AbstractWakeDeficitModel end

"""
    JensenTopHat(alpha)

Container for parameters related to the Jensen Top Hat deficit model

# Arguments
- `alpha::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.1
"""
struct JensenTopHat{TF} <: AbstractWakeDeficitModel
    alpha::TF
end
JensenTopHat() = JensenTopHat(0.2)

"""
    JensenCosine(alpha)

Container for parameters related to the Jensen Cosine deficit model

# Arguments
- `alpha::Float`: parameter controlling the wake deficit decay rate. Default value is 0.1
- `beta::Float`: parameter controlling the width of the cosine function. Default value is 20.0 deg., given in radians.
"""
struct JensenCosine{TF,ATF} <: AbstractWakeDeficitModel
    alpha::TF
    beta::TF
    wec_factor::ATF
end
JensenCosine() = JensenCosine(0.1, 20.0*pi/180.0, [1.0])
JensenCosine(x) = JensenCosine(x, 20.0*pi/180.0, [1.0])
JensenCosine(x, y) = JensenCosine(x, y, [1.0])

"""
    MultiZone(me, ke, MU, aU, bU)

Container for parameters related to the MultiZone deficit model

# Arguments
- `me::Float`: parameter controlling general wake expansion. Default value is 0.065
- `ke::Array{Float}(3)`: parameters controlling the wake expansion of each zone respectively. Default values are [-0.5 0.22 1.0].
- `MU::Array{Float}(3)`: parameters controlling the wake deficit decay of each zone respectively. Default values are [0.5 1.0 5.5].
- `aU::Float`: parameter impacting the wake deficit decay for a constant wake deflection. Default value is 5.0.
- `bU::Float`: parameter changing the wake deficit decay under yawed conditions. Default value is 1.66.
"""
struct MultiZone{ATF, TF} <: AbstractWakeDeficitModel
    me::ATF
    ke::TF
    MU::ATF
    aU::TF
    bU::TF
end
MultiZone() = MultiZone([-0.5 0.22 1.0], 0.065, [0.5 1.0 5.5], 5.0, 1.66)

"""
    GaussOriginal(k_star)

Container for parameters related to the origina Gaussian deficit model presented by Bastankhah and Porte-Agel 2014

# Arguments
- `k_star::Float`: parameter controlling the wake spreading rate and deficit decay. Default value is 0.075
"""
struct GaussOriginal{TF} <: AbstractWakeDeficitModel
    k_star::TF
end
GaussOriginal() = GaussOriginal(0.075)

"""
    GaussYaw(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star, interpolation)

Container for parameters related to the Gaussian deficit model with yaw presented by Bastankhah and Porte-Agel 2016

# Arguments
- `horizontal_spread_rate::Float`: parameter controlling the horizontal spread of the deficit model. Default value is 0.022.
- `vertical_spread_rate::Float`: parameter controlling the vertical spread of the deficit model. Default value is 0.022.
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
- `interpolation::Bool`: boolean stating if the the near wake should be interpolated. Default value is true.
"""
struct GaussYaw{TF, ATF, BO} <: AbstractWakeDeficitModel
    horizontal_spread_rate::TF
    vertical_spread_rate::TF
    alpha_star::TF
    beta_star::TF
    wec_factor::ATF
    interpolate_sigma::BO
end
GaussYaw() = GaussYaw(0.022, 0.022, 2.32, 0.154, [1.0], true)
GaussYaw(interp) = GaussYaw(0.022, 0.022, 2.32, 0.154, [1.0], interp)
GaussYaw(a, b, c, d) = GaussYaw(a, b, c, d, [1.0], true)
GaussYaw(a, b, c, d, interp) = GaussYaw(a, b, c, d, [1.0], interp)

"""
    GaussYawVariableSpread(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star, interpolation)

Container for parameters related to the Gaussian deficit model with yaw presented by Bastankhah and Porte-Agel 2016
    and the farm model presented by Niayifar and Porte-Agel in 2016.

# Arguments
- `alpha_star::Float`: parameter controlling the impact of turbulence intensity on the length of the near wake. Default value is 2.32.
- `beta_star::Float`: parameter controlling the impact of the thrust coefficient on the length of the near wake. Default value is 0.154.
- `interpolation::Bool`: boolean stating if the the near wake should be interpolated. Default value is true.
"""
struct GaussYawVariableSpread{TF, ATF, BO} <: AbstractWakeDeficitModel
    alpha_star::TF
    beta_star::TF
    k1::TF
    k2::TF
    wec_factor::ATF
    interpolate_sigma::BO
end
GaussYawVariableSpread() = GaussYawVariableSpread(2.32, 0.154, 0.3837, 0.003678, [1.0], true)
GaussYawVariableSpread(interp) = GaussYawVariableSpread(2.32, 0.154, 0.3837, 0.003678, [1.0], interp)
GaussYawVariableSpread(x, y, z) = GaussYawVariableSpread(x, y, 0.3837, 0.003678, z, true)
GaussYawVariableSpread(x, y, z, interp) = GaussYawVariableSpread(x, y, 0.3837, 0.003678, z, interp)
GaussYawVariableSpread(x, y) = GaussYawVariableSpread(x, y, 0.3837, 0.003678, [1.0], true)

"""
    GaussSimple(k, wec_factor)

Container for parameters related to the Gaussian deficit model with yaw presented by Bastankhah and Porte-Agel 2016

# Arguments
- `k::Float`: parameter controlling the spread of the wake
- `wec_factor::Array{Float}`: parameter artificial wake spreading for wake expansion continuation (WEC) optimization
"""
struct GaussSimple{TF, ATF} <: AbstractWakeDeficitModel
    k::TF
    wec_factor::ATF
end
GaussSimple(k) = GaussSimple(k, [1.0])
GaussSimple() = GaussSimple(0.0324555, [1.0])

"""
    CumulativeCurl(a_s, b_s, c_s1, c_s2, a_f, b_f, c_f)

Container for parameters related to the Cumulative Curl model used in FLORIS v3 as shown in Bay 2022 (https://doi.org/10.5194/wes-2022-17)

# Arguments
- `a_s::Float`: parameter relating turbulence intensity to wake spread. Default value is 0.179367259
- `b_s::Float`: parameter controlling the default wake spread. Default value is 0.0118889215
- `c_s1::Float`: parameter relating Ct to initial wake width. Default value is 0.0563691592
- `c_s2::Float`: parameter affecting initial wake width. Default value is 0.13290157
- `a_f::Float`: Default value is 3.11
- `b_f::Float`: Default value is -0.68
- `c_f::Float`: Default value is 2.41
- 'wec_factor': paramter for the wake expansion continuation Default is [1.0]
"""
struct CumulativeCurl{TF,ATF} <: AbstractWakeDeficitModel
    a_s::TF
    b_s::TF
    c_s1::TF
    c_s2::TF
    a_f::TF
    b_f::TF
    c_f::TF
    wec_factor::ATF
end
CumulativeCurl() = CumulativeCurl(0.179367259, 0.0118889215, 0.0563691592, 0.13290157, 3.11, -0.68, 2.41, [1.0])

"""
    NoWakeDeficit()
"""
struct NoWakeDeficit <: AbstractWakeDeficitModel
end


"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model::JensenTopHat)

Computes the wake deficit according to the original Jensen top hat wake model, from the paper:
"A Note on Wind Generator Interaction" by N.O. Jensen (1983)

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::JensenTopHat`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::JensenTopHat)


    r0 = rotor_diameter[upstream_turbine_id]/2.0 #turbine rotor radius


    if downstream_turbine_id == 0

        # find delta x, y, and z. dx is the downstream distance from the turbine to
        # the point of interest. dy and dz are the distances from the point of interest
        # and the wake center (in y and z)
        dx = locx-turbine_x[upstream_turbine_id]
        dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)
        dz = locz-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

        del = norm([dy, dz]) #distance from wake center to the point of interest
        r = model.alpha*dx + r0 #figure (1) from the paper

        if (dx < 0.0) || (del > r)
            totloss = 0.0
        else
            loss = 2.0*turbine_ai[upstream_turbine_id]*(r0/r)^2 #equation (2) from the paper
            totloss = loss
        end
    else
        r1 = rotor_diameter[downstream_turbine_id]/2.0
        # find delta x, y, and z. dx is the downstream distance from the upstream turbine hub to
        # the downstream turbine hub of interest. dy and dz are the distances from the point of interest
        # and the wake center (in y and z)
        dxt = turbine_x[downstream_turbine_id]-turbine_x[upstream_turbine_id]
        dyt = turbine_y[downstream_turbine_id]-(turbine_y[upstream_turbine_id]+deflection_y)
        dzt = turbine_z[downstream_turbine_id]-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

        r = model.alpha*dxt + r0 #figure (1) from the paper

        posz = turbine_z[downstream_turbine_id] + hub_height[downstream_turbine_id]

        # finds the wakediameter
        wake_center_y = (turbine_y[upstream_turbine_id]+deflection_y)
        wake_center_z = (turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)
        wake_diameter = 2.0*r

        loss = 2.0*turbine_ai[upstream_turbine_id]*(r0/r)^2 #equation (2) from the paper

        overlap = overlap_area_func(turbine_y[downstream_turbine_id], posz, 2.0*r1, wake_center_y, wake_center_z, wake_diameter; tol=1E-6)

        ovlpercentage = overlap/(pi*r1^2)

        totloss = loss*ovlpercentage

    end

    return totloss
end


"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model::JensenCosine)

Computes the wake deficit according to the original Jensen cosine wake model, from the paper:
"A Note on Wind Generator Interaction" by N.O. Jensen (1983)

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::JensenCosine`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::JensenCosine)
    """the original Jensen cosine wake model, from the paper: "A Note on Wind
    Generator Interaction" by N.O. Jensen (1983)"""

    # get wec factor (See Thomas and Ning 2018)
    wec_factor = model.wec_factor[1]

    # find delta x, y, and z. dx is the downstream distance from the turbine to
    # the point of interest. dy and dz are the distances from the point of interest
    # and the wake center (in y and z)
    dx = locx-turbine_x[upstream_turbine_id]
    dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)
    dz = locz-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

    r0 = rotor_diameter[upstream_turbine_id]/2.0 #turbine rotor radius
    del = norm([dy, dz]) #distance from wake center to the point of interest

    if dx < 0.
        loss = 0.0 # no loss outside the wake
    else
        d = wec_factor*r0/tan(model.beta) # distance from fulcrum of wake cone to wind turbine hub
        theta = atan(del/(dx+d)) # angle from center of wake to point of interest
        if theta > model.beta # if you're outside the wake
            loss = 0.0
        else # if you're inside the wake
            n = pi / model.beta # see Jensen 1983 eq. 3. Value used for n in paper was 9, corresponding to beta = 20.0 deg.
            ftheta = (1.0 + cos(n * theta)) / 2.0 # cosine term to be applied to loss equation as per Jensen 1983
            loss = 2.0*turbine_ai[upstream_turbine_id]*(ftheta*r0/(r0+model.alpha*dx))^2 #equation (2) from the paper
        end
    end

    return loss
end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model::MultiZone)

Computes the wake deficit at a given location using the original MultiZone "FLORIS" wake model, from the paper:
"Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study" by Gebraad et al. (2014)

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::MultiZone`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::MultiZone)

    dt = rotor_diameter[upstream_turbine_id]



    # extract model parameters
    ke = model.ke
    me = model.me
    MU = model.MU
    aU = model.aU
    bU = model.bU


    # aU is given in degrees when it should be radians
    aU = deg2rad(aU)


    if downstream_turbine_id == 0             # Finds the velocity at that point if not looking for a turbine

        # find delta x, y, and z. dx is the downstream distance from the turbine to
        # the point of interest. dy and dz are the distances from the point of interest
        # and the wake center (in y and z)

        dx = locx-turbine_x[upstream_turbine_id]
        dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)
        dz = locz-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)


        if dx < 0.
            c = 0.0
        else
            del = norm([dy, dz]) #distance from wake center to the point of interest

            # calculate the diameter of the wake in each of the three zones (at the specified dx)
            Dw = zeros(3)
            for i = 1:3
                Dw[i] = max(dt+2*ke*me[i]*dx,0) # equation (13) from the paper
            end
            Rw = Dw./2 # radius of the wake zones

            # calculate the coefficient c of the appropriate wake zone:
            if del > Rw[3]
                c = 0.0
            else
                # equations (15, 16, and 17) from the paper
                if del < Rw[1]
                    MUi = MU[1]
                elseif del < Rw[2]
                    MUi = MU[2]
                elseif del < Rw[3]
                    MUi = MU[3]
                end

                mU = MUi/(cos(aU+bU*turbine_yaw[upstream_turbine_id]))
                c = (dt/(dt+2.0*ke*mU*dx))^2
            end


        end

        loss = 2.0*turbine_ai[upstream_turbine_id]*c # calculate the wake loss. Equation (14) fro mthe paper

    else            # this section will calculate the wind effect on an entire turbine.

        # find delta x, y, and z. dx is the downstream distance from the turbine to
        # the point of interest. dy and dz are the distances from the point of interest
        # and the wake center (in y and z)
        dxt = turbine_x[downstream_turbine_id]-turbine_x[upstream_turbine_id]
        dyt = turbine_y[downstream_turbine_id]-(turbine_y[upstream_turbine_id]+deflection_y)
        dzt = (turbine_z[downstream_turbine_id]+hub_height[downstream_turbine_id])-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

        posz = turbine_z[downstream_turbine_id] + hub_height[downstream_turbine_id] # finds the z coordinate of the turbine hub

        if dxt < 0.
            c1,c2,c3 = 0,0,0
        else
            del = norm([dxt, dzt]) #distance from wake center to the point of interest

            wake_center_y = (turbine_y[upstream_turbine_id]+deflection_y)
            wake_center_z = (turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

            # calculate the diameter of the wake in each of the three zones (at the specified dx)
            # after calculating the size of the wake the model will then find its overlap with the turbine

            Dw = zeros(3)
            ovlp = zeros(3)
            for i = 1:3
                Dw[i] = max(dt+2*ke*me[i]*dxt,0) # equation (13) from the paper
                overlap = overlap_area_func(turbine_y[downstream_turbine_id], posz,
                rotor_diameter[downstream_turbine_id], wake_center_y, wake_center_z, Dw[i]; tol=1E-6)
                area = pi*(rotor_diameter[downstream_turbine_id]^2)/4
                ovlp[i] = overlap/area
            end

            Rw = Dw./2 # radius of the wake zones

            # equations (15, 16, and 17) from the paper calculated for all 3 zones
            # losses computed for each zone and the subsequently added together


            mU1 = MU[1]/(cos(aU+bU*turbine_yaw[upstream_turbine_id]))
            c1 = (dt/(dt+2.0*ke*mU1*dxt))^2

            mU2 = MU[2]/(cos(aU+bU*turbine_yaw[upstream_turbine_id]))
            c2 = (dt/(dt+2.0*ke*mU2*dxt))^2

            mU3 = MU[3]/(cos(aU+bU*turbine_yaw[upstream_turbine_id]))
            c3 = (dt/(dt+2.0*ke*mU3*dxt))^2

        end

        # losses calculated for all 3 zones with equation [14] from the paper
        loss1 = 2.0*turbine_ai[upstream_turbine_id]*c1*ovlp[1]
        loss2 = 2.0*turbine_ai[upstream_turbine_id]*c2*ovlp[2]
        loss3 = 2.0*turbine_ai[upstream_turbine_id]*c3*ovlp[3]

        # losses of 3 zones are summed to get total loss
        loss = loss1 + loss2 + loss3
    end


    return loss
end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model::GaussOriginal)

Computes the wake deficit at a given location using the Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "A new analytical model for wind-turbine wakes" (2014)

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::GaussOriginal`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussOriginal)

    dx = locx-turbine_x[upstream_turbine_id]
    dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)
    dz = locz-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

    # extract turbine properties
    dt = rotor_diameter[upstream_turbine_id]
    ct = turbine_ct[upstream_turbine_id]

    # extract model parameters
    ks = model.k_star       # wake spread rate (k* in 2014 paper)

    # calculate beta (paper eq: 6)
    beta = 0.5*(1.0 + sqrt(1.0 - ct))/sqrt(1.0 - ct)

    # use 2 as length of potential core #TODO calculate this
    x0 = 2.0*dt

    # calculate loss (paper eq: 23)
    num = (dz/dt)^2 + (dy/dt)^2
    if dx > x0
        denom = (ks*dx/dt + 0.2*sqrt(beta))^2
    else
        denom = (ks*x0/dt + 0.2*sqrt(beta))^2
    end

    loss = (1.0 - sqrt(1.0 - ct/(8.0*denom)))*exp(-num/(2.0*denom))

    return loss

end

"""
    _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)

Helper function for wake_deficit_model when using the GaussYaw model. Computes the length of the near wake potential core.
"""
function _gauss_yaw_potential_core(d, yaw, ct, as, ti, bs)
    # from Bastankhah and Porte-Agel 2016 eqn 7.3
    x0 = d*(cos(yaw)*(1.0+sqrt(1.0-ct)))/(sqrt(2.0)*(as*ti+bs*(1.0-sqrt(1.0-ct))))
    return x0
end

"""
    _gauss_yaw_spread(dt, k, dx, x0, yaw)

Helper function for wake_deficit_model when using the GaussYaw model. Computes the standard deviation of the wake.
"""
function _gauss_yaw_spread(dt, k, dx, x0, yaw)
    # from Bastankhah and Porte-Agel 2016 eqn 7.2
    sigma = k*(dx - x0) + dt*cos(yaw)/sqrt(8.0)

    return sigma

end

"""
    _gauss_yaw_spread_interpolated(dt, k, dx, x0, yaw)

Helper function for wake_deficit_model when using the GaussYaw model. Computes the standard deviation of the wake.
with an interpolation on the near wake.

"""
function _gauss_yaw_spread_interpolated(dt, k, dx, x0, yaw, xd; interpolate=true)
    # calculate wake spread

    if interpolate
        if dx > x0 # far wake
            sigma = _gauss_yaw_spread(dt, k, dx, x0, yaw)
        else # linear interpolation in the near wakes
            dx_interpolate = xd+((x0-xd)/(x0))*(dx)
            sigma = _gauss_yaw_spread(dt, k, dx_interpolate, x0, yaw)
        end
    else
        if dx > xd # far wake
            sigma = _gauss_yaw_spread(dt, k, dx, x0, yaw)
        else # use sigma at xd for undefined values
            sigma = _gauss_yaw_spread(dt, k, xd, x0, yaw)
        end
    end

    return sigma
end

function _gauss_yaw_discontinuity(dt, x0, ky, kz, yaw, ct)
    # for clarity, break out the terms in the equation
    a = ky + kz*cos(yaw)
    b = 4.0 * ky * kz * cos(yaw)*(ct - 1.0)
    c = 2.0 * sqrt(8.0) * ky * kz

    # distance from rotor to the last point where the wake model is undefined
    discontinuity_point = x0 + dt * (a - sqrt(a^2 - b))/c

    return discontinuity_point
end

function _gauss_yaw_model_deficit(dx, dy, dz, dt, yaw, ct, ti, as, bs, ky, kz, wf)

    if dx > 1E-6 # loss in the wake

        # println("x0 inputs:", diam, " ", yaw, " ", ct, " ", as, " ", ti, " ", bs)
        # calculate the length of the potential core (paper eq: 7.3)
        x0 = _gauss_yaw_potential_core(dt, yaw, ct, as, ti, bs)

        # calculate the discontinuity point of the gauss yaw model
        xd = _gauss_yaw_discontinuity(dt, x0, ky, kz, yaw, ct)

        # calculate horizontal wake spread (paper eq: 7.2)
        sigma_y = _gauss_yaw_spread_interpolated(dt, ky, dx, x0, yaw, xd)

        # calculate vertical wake spread (paper eq: 7.2)
        sigma_z = _gauss_yaw_spread_interpolated(dt, kz, dx, x0, 0.0, xd)

        # calculate velocity deficit #check - infty when large input ~= 500
        ey = exp(-0.5*(dy/(wf*sigma_y))^2)
        ez = exp(-0.5*(dz/(wf*sigma_z))^2)

        sqrtterm = 1.0-ct*cos(yaw)/(8.0*(sigma_y*sigma_z/dt^2))
        if sqrtterm >= 1e-8 #check - could try increasing tolerance
            loss = (1.0-sqrt(sqrtterm))*ey*ez
        else
            loss = ey*ez
        end
    else # loss upstream of turbine
        loss = 0.0
    end


    return loss

end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model::GaussYaw)

Computes the wake deficit at a given location using the The Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "Experimental and theoretical study of wind turbine wakes in yawed conditions" (2016)

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::GaussYaw`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussYaw)

    dx = locx-turbine_x[upstream_turbine_id]
    dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)
    dz = locz-(turbine_z[upstream_turbine_id]+hub_height[upstream_turbine_id]+deflection_z)

    # extract turbine properties
    dt = rotor_diameter[upstream_turbine_id]
    yaw = turbine_yaw[upstream_turbine_id]
    ct = turbine_ct[upstream_turbine_id]

    # extract model parameters
    # ks = model.k_star       # wake spread rate (k* in 2014 paper)
    ti = turbine_local_ti[upstream_turbine_id]
    ky = model.horizontal_spread_rate
    kz = model.vertical_spread_rate
    as = model.alpha_star
    bs = model.beta_star
    wec_factor = model.wec_factor[1]

    loss = _gauss_yaw_model_deficit(dx, dy, dz, dt, yaw, ct, ti, as, bs, ky, kz, wec_factor)

    return loss

end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussYawVariableSpread)

Computes the wake deficit at a given location using the The Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "Experimental and theoretical study of wind turbine wakes in yawed conditions" (2016)
The spread rate is adjusted based on local turbulence intensity as in Niayifar and Porte-Agel 2016

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::GaussYawVariableSpread`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussYawVariableSpread)

    dx = locx - turbine_x[upstream_turbine_id]
    dy = locy - (turbine_y[upstream_turbine_id] + deflection_y)
    dz = locz - (turbine_z[upstream_turbine_id] + hub_height[upstream_turbine_id] + deflection_z)

    # extract turbine properties
    dt = rotor_diameter[upstream_turbine_id]
    yaw = turbine_yaw[upstream_turbine_id]
    ct = turbine_ct[upstream_turbine_id]

    # extract model parameters
    # ks = model.k_star       # wake spread rate (k* in 2014 paper)
    ti = turbine_local_ti[upstream_turbine_id]
    ky = kz = _k_star_func(ti, model.k1, model.k2)

    as = model.alpha_star
    bs = model.beta_star
    wec_factor = model.wec_factor[1]

    loss = _gauss_yaw_model_deficit(dx, dy, dz, dt, yaw, ct, ti, as, bs, ky, kz, wec_factor)

    return loss
end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussSimple)

Computes the wake deficit at a given location using the Gaussian wake model presented by Bastankhah and Porte-Agel in the paper: "A new analytical model for wind-turbine wakes" (2014)
    as modified for IEA Task 37 Case Studies 3 and 4

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `model::GaussSimple`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::GaussSimple)

    dx = locx-turbine_x[upstream_turbine_id]
    dy = locy-(turbine_y[upstream_turbine_id]+deflection_y)

    # extract turbine properties
    dt = rotor_diameter[upstream_turbine_id]
    ct = turbine_ct[upstream_turbine_id]

    # extract model properties
    k = model.k
    wf = model.wec_factor[1]

    # calculate loss
    sigmay = k*dx + dt/sqrt(8.0)
    radical = 1.0 - ct/(8.0*(sigmay^2)/(dt^2))
    exponent = -0.5*(dy/(wf*sigmay))^2
    loss = (1.0 - sqrt(radical))*exp(exponent)

    return loss

end

"""
    wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::CumulativeCurl)

Computes the wake deficit at a given location using the Cumulative Curl Model https://doi.org/10.5194/wes-2022-17

# Arguments
- `locx::Float`: x coordinate where wind speed is calculated
- `locy::Float`: y coordinate where wind speed is calculated
- `locz::Float`: z coordinate where wind speed is calculated
- `turbine_x::Array(Float)`: vector containing x coordinates for all turbines in farm
- `turbine_y::Array(Float)`: vector containing y coordinates for all turbines in farm
- `turbine_z::Array(Float)`: vector containing z coordinates for all turbines in farm
- `deflection_y::Float`: deflection in the y direction of downstream wake
- `deflection_z::Float`: deflection in the z direction of downstream wake
- `upstream_turbine_id::Int`: index of the upstream wind turbine creating the wake
- `downstream_turbine_id::Int`: index of the downstream turbine feeling the wake (if not referencing a turbine set to zero)
- `hub_height::Array(Float)`: vector containing hub heights for all turbines in farm
- `rotor_diameter::Array(Float)`: vector containing rotor diameters for all turbines in farm
- `turbine_ai::Array(Float)`: vector containing initial velocity deficits for all turbines in farm
- `turbine_local_ti::Array(Float)`: vector containing local turbulence intensities for all turbines in farm
- `turbine_ct::Array(Float)`: vector containing thrust coefficients for all turbines in farm
- `turbine_yaw::Array(Float)`: vector containing the yaw angle? for all turbines in farm
- `wake_deficits`: matrix containing the wake deficits from every turbine to every other turbine for use in sparsity codes
- `contribution_matrix`: matrix containing the contribution coefficients defined in the Cumlative Curl Model
- `model::CumulativeCurl`: indicates the wake model in use

"""
function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::CumulativeCurl)
    x_n = turbine_x[upstream_turbine_id]
    y_n = turbine_y[upstream_turbine_id]
    z_n = turbine_z[upstream_turbine_id] + hub_height[upstream_turbine_id]

    # unpack model constants
    a_f = model.a_f
    b_f = model.b_f
    c_f = model.c_f
    wec_factor = model.wec_factor[1]

    x_tilde_n = abs(locx - x_n) / rotor_diameter[upstream_turbine_id]
    m = a_f*exp(b_f*x_tilde_n)+c_f
    a1 = 2^(2/m - 1)
    a2 = a1^2

    sigma_squared[upstream_turbine_id,downstream_turbine_id] = wake_expansion(turbine_ct[upstream_turbine_id],turbine_local_ti[upstream_turbine_id],x_tilde_n,model)

    sigma_n = sigma_squared[upstream_turbine_id,downstream_turbine_id]
    sqrt_num = (y_n-locy-deflections[upstream_turbine_id,downstream_turbine_id])^2 + (z_n-locz)^2
    if sqrt_num == 0
        r_tilde = 0.0
    else
        r_tilde = sqrt(sqrt_num)/rotor_diameter[upstream_turbine_id]
    end
    if (-1 * (r_tilde^m)/(2.0*sigma_n*wec_factor)) > -750

        sum_C = 0.0

        @inbounds begin
            for i = 1:current_index_loop-1
                other_turbine_id = Int(sorted_turbine_index[i])

                if contribution_matrix[other_turbine_id,downstream_turbine_id] == 0.0
                    continue
                end

                y_i = turbine_y[other_turbine_id]
                z_i = turbine_z[other_turbine_id] + hub_height[other_turbine_id]
                sigma_i = sigma_squared[other_turbine_id,downstream_turbine_id]
                dy_i = deflections[other_turbine_id,downstream_turbine_id]

                @fastmath sigma_plus = 1/(sigma_n+sigma_i)
                @fastmath expo = -0.5*sigma_plus*((y_n-y_i-dy_i)^2 + (z_n-z_i)^2)
                if expo > -750
                    lambda = sigma_n*sigma_plus*exp(expo)
                    sum_C += lambda * contribution_matrix[other_turbine_id,downstream_turbine_id]
                end
            end

            calc = abs_smooth(a2 - (m*turbine_ct[upstream_turbine_id]*cos(turbine_yaw[upstream_turbine_id]))/(16.0*gamma(2.0/m)*(sigma_n^(2.0/m))*(1.0-sum_C/wind_speed_internal)^2),0.1)
            contribution_matrix[upstream_turbine_id,downstream_turbine_id] = (1-sum_C/wind_speed_internal) * (a1-sqrt(calc))
            wake_deficits[upstream_turbine_id,downstream_turbine_id] = contribution_matrix[upstream_turbine_id,downstream_turbine_id] * exp(-1 * (r_tilde^m)/(2.0*sigma_n*wec_factor))

            # yaw color plot
            # keep = [155 114 14 5 1 2 8 20 39 64 95 133 177 228 285]
            # cols = [317 155 114 14 5 1 2 8 20 39 64 95 133 177 228]
            # for j = eachindex(cols)
            #     temp_def = wake_deficits[cols[j],keep[j]]
            #     temp_c = contribution_matrix[cols[j],keep[j]]
            #     wake_deficits[cols[j],:] .= 0.0
            #     contribution_matrix[cols[j],:] .= 0.0
            #     wake_deficits[cols[j],keep[j]] = temp_def
            #     contribution_matrix[cols[j],keep[j]] = temp_c
            # end

            # # 317 -> 155 -> 114 -> 14 -> 5 -> 1 -> 2 -> 8 -> 20 -> 39 -> 64 -> 95 -> 133 -> 177 -> 228 -> 285

            return wake_deficits[upstream_turbine_id,downstream_turbine_id]
        end
    end
    wake_deficits[upstream_turbine_id,downstream_turbine_id] = 0.0
    return wake_deficits[upstream_turbine_id,downstream_turbine_id]
end

# Helper function for the Cumulative Curl model
function wake_expansion(Ct,TI,x_tilde,model)
    beta = 0.5*(1.0+sqrt(1.0-Ct))/sqrt(1.0-Ct)
    epsilon = (model.c_s1*Ct+model.c_s2)*sqrt(beta)
    k = (model.a_s*TI+model.b_s)
    sigma = k*x_tilde+epsilon
    return sigma^2
end

function wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wake_deficits, contribution_matrix, deflections, current_index_loop, wind_speed_internal, sigma_squared, wtvelocities, sorted_turbine_index, model::NoWakeDeficit)
    return 0.0
end
