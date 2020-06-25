using FlowFarm
using CCBlade
using PyPlot

const ff=FlowFarm


"""
    multiple_components_op(U, V, W, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

Return the operating points along the blade considering all the inflow velocity components.

# Arguments
- `U::Array{Float}`: u velocity component of the inflow at each r
- `V::Array{Float}`: v velocity component of the inflow at each r
- `W::Array{Float}`: w velocity component of the inflow at each r
- `Omega::Float`: rotor rotational speed
- `r::Array{Float}`: radial locations of interest
- `precone::Float`: rotor precone angle
- `yaw::Float`: rotor yaw angle
- `tilt::Float`: rotor tilt angle
- `azimuth::Float`: blade azimuth angle
- `rho::Float`: air density

# Keyword Arguments
- `mu::Float`: air viscocity (can usually use the default)
- `asound::Float`: speed of sound (can usually use the default)
"""
function multiple_components_op(U, V, W, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    sc = sin(precone)
    cc = cos(precone)

    magnitude = sqrt.(U.^2 + V.^2 + W.^2)
    x0 = U./magnitude
    y0 = V./magnitude
    z0 = W./magnitude

    #first
    x1 = x0.*cy + y0.*sy
    y1 = -x0.*sy + y0.*cy
    z1 = z0

    #second
    x2 = -z1.*st + x1.*ct
    y2 = y1
    z2 = z1.*ct + x1.*st

    #third
    x3 = x2
    y3 = y2.*ca + z2.*sa
    z3 = -y2.*sa + z2.*ca

    #fourth
    x4 = z3.*sc + x3.*cc
    y4 = y3
    z4 = z3.*cc - x3.*sc

    # transform wind to blade c.s.
    Vwind_x = magnitude.*x4
    Vwind_y = magnitude.*y4
    Vwind_z = magnitude.*z4

    # coordinate in azimuthal coordinate system
    z_az = r.*cos(precone)

    # wind from rotation to blade c.s.
    Vrot_x = -Omega*sc
    Vrot_y = Omega.*z_az

    # total velocity
    Vx = Vwind_x .+ Vrot_x
    Vy = Vwind_y + Vrot_y

    # operating point
    return OperatingPoint(Vx, Vy, rho, mu, asound)
end


"""
    distributed_velocity_op(V, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

Return the operating points along the blade considering varied inflow along the blade.

# Arguments
- `V::Array{Float}`: velocity inflow at each r
- `Omega::Float`: rotor rotational speed
- `r::Array{Float}`: radial locations of interest
- `precone::Float`: rotor precone angle
- `yaw::Float`: rotor yaw angle
- `tilt::Float`: rotor tilt angle
- `azimuth::Float`: blade azimuth angle
- `rho::Float`: air density

# Keyword Arguments
- `mu::Float`: air viscocity (can usually use the default)
- `asound::Float`: speed of sound (can usually use the default)
"""
function distributed_velocity_op(V, Omega, r, precone, yaw, tilt, azimuth, rho, mu=1.81206e-05, asound=1.0)

    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    sc = sin(precone)
    cc = cos(precone)

    # coordinate in azimuthal coordinate system
    x_az = -r.*sin(precone)
    z_az = r.*cos(precone)

    # get section heights in wind-aligned coordinate system
    heightFromHub = (z_az*ca)*ct - x_az*st

    # transform wind to blade c.s.
    Vwind_x = V * ((cy*st*ca + sy*sa)*sc + cy*ct*cc)
    Vwind_y = V * (cy*st*sa - sy*ca)

    # wind from rotation to blade c.s.
    Vrot_x = -Omega*sc
    Vrot_y = Omega*z_az

    # total velocity
    Vx = Vwind_x + Vrot_x
    Vy = Vwind_y + Vrot_y

    # operating point
    return OperatingPoint(Vx, Vy, rho, mu, asound)
end


"""
    find_xyz_simple(x_hub,y_hub,z_hub,r,yaw,azimuth)

Find the xyz locations of points along a blade given it's location and azimuth angle.
Currently doesn't consider precone or tilt.

# Arguments
- `x_hub::Float`: x location of hub
- `y_hub::Float`: y location of hub
- `z_hub::Float`: z location of hub (hub height if no topology)
- `r::Array{Float}`: radial locations of interest
- `precone::Float`: rotor precone angle
- `yaw::Float`: rotor yaw angle
- `azimuth::Float`: blade azimuth angle
"""
function find_xyz_simple(x_hub,y_hub,z_hub,r,yaw,azimuth)

        sy = sin(yaw)
        cy = cos(yaw)
        sa = sin(azimuth)
        ca = cos(azimuth)

        x_locs = x_hub .- r.*(sy*sa)
        y_locs = y_hub .- r.*(cy*sa)
        z_locs = z_hub .+ r.*(ca)

        return x_locs, y_locs, z_locs
end


"""
    calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

Calculates stresses from bending moments on a hollow cylinder

# Arguments
- `mx::Float`: x moment
- `my::Float`: y moment
- `dx::Float`: x distance to the location of interest
- `dy::Float`: y distance to the location of interest

# Keyword Arguments
- `Rcyl::Float`: radius of the cylinder
- `tcyl::Float`: thickenss of the cylinder
"""
function calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

        I = 0.25*pi*(Rcyl^4-(Rcyl-tcyl)^4)
        stress = mx*dy/I + my*dx/I

        return stress
end


function calc_moment_stress_resultant(edge,flap,dx,dy,Rcyl=1.771,tcyl=0.06)

        I = 0.25*pi*(Rcyl^4-(Rcyl-tcyl)^4)

        magnitude = sqrt(edge^2 + flap^2)
        distance = abs(dy*edge - dx*flap)/sqrt(flap^2 + edge^2)
        stress = magnitude*distance/I

        return stress
end


function get_speeds(turbine_x,turbine_y,turbine_z,turb_index,hub_height,r,turbine_yaw,azimuth,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
                    sorted_turbine_index,wtvelocities,wind_resource,model_set;wind_farm_state_id=1,downwind_turbine_id=0)
        x_locs, y_locs, z_locs = find_xyz_simple(turbine_x[turb_index],turbine_y[turb_index],hub_height[turb_index],r,turbine_yaw[turb_index],azimuth)
        npts = length(x_locs)
        arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                                typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))
        point_velocities = zeros(arr_type,npts)
            for i in 1:npts
                # loc = [x_locs[i], y_locs[i], z_locs[i]]
                # point_velocities[i] = FlowFarm.point_velocity(loc, windfarm, windfarmstate, windresource, wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, 0)

                # point_velocities[i] = FlowFarm.point_velocity(loc, model_set, problem_description)

                point_velocities[i] = FlowFarm.point_velocity(x_locs[i], y_locs[i], z_locs[i], turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                                    wind_resource, model_set,
                                    wind_farm_state_id=1, downwind_turbine_id=0)



                # loc = [x_locs[i], y_locs[i], z_locs[i]]
                # point_velocities[i] = FlowFarm.point_velocity(loc, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                #                     rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                #                     wind_resource, model_set;
                #                     wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)
            end
        return point_velocities
end


"""

    rainflow(array_ext,uc_mult=0.5)
Rainflow counting of a signal's turning points

# Arguments
        array_ext (numpy.ndarray): array of turning points

# Keyword Arguments
        uc_mult (float): partial-load scaling [opt, default=0.5]

# Returns
        array_out (numpy.ndarray): (3 x n_cycle) array of rainflow values:
                                    1) load range
                                    2) range mean
                                    3) cycle count
"""
function rainflow(array_ext,uc_mult=0.5)

    tot_num = length(array_ext)             # total size of input array
    array_out = zeros(typeof(array_ext[1]),(3, tot_num-1))        # initialize output array

    pr = 1                                  # index of input array
    po = 1                                  # index of output array
    j = 0                                   # index of temporary array "a"
    a  = zeros(typeof(array_ext[1]),tot_num)                     # temporary array for algorithm

    # loop through each turning point stored in input array
    for i = 1:tot_num
        j += 1                  # increment "a" counter
        a[j] = array_ext[pr]    # put turning point into temporary array
        pr += 1                 # increment input array pointer

        while j >= 3 && abs( a[j-1] - a[j-2]) <= abs(a[j] - a[j-1])

            lrange = abs( a[j-1] - a[j-2] )

            # partial range
            if j == 3
                mean      = ( a[1] + a[2] ) / 2.0
                a[1]=a[2]
                a[2]=a[3]
                j=2
                if lrange > 0
                    array_out[1,po] = lrange
                    array_out[2,po] = mean
                    array_out[3,po] = uc_mult
                    po += 1
                end

            # full range
            else
                mean      = ( a[j-1] + a[j-2] ) / 2.0
                a[j-2]=a[j]
                j=j-2
                if (lrange > 0)
                    array_out[1,po] = lrange
                    array_out[2,po] = mean
                    array_out[3,po] = 1.00
                    po += 1
                end
            end
        end
    end

    # partial range
    for i = 1:j-1
        lrange    = abs( a[i] - a[i+1] )
        mean      = ( a[i] + a[i+1] ) / 2.0
        if lrange > 0
            array_out[1,po] = lrange
            array_out[2,po] = mean
            array_out[3,po] = uc_mult
            po += 1
        end
    end

    # get rid of unused entries
    out = array_out[:,1:po-1]

    return out
end


"""
    get_peaks(array)

get the turning point values of a signal

# Arguments
- `array::Array{Float}`: the signal to find the turning points, or peaks
"""
function get_peaks(array)
    A = array[:]
    # get rid of any zero slope in the beginning
    while A[2] == A[1]
        A = A[2:length(A)]
    end
    peaks = [A[1]]
    if A[2] > A[1]
        slope = "p"
    elseif A[2] < A[1]
        slope = "m"
    end
    for i = 1:length(A)-2
        ind = i+1
        if slope == "p"
            if A[ind+1] < A[ind]
                peaks = append!(peaks,A[ind])
                slope = "m"
            end
        elseif slope == "m"
            if A[ind+1] > A[ind]
                peaks = append!(peaks,A[ind])
                slope = "p"
            end
        end
    end
    peaks = append!(peaks,A[length(A)])
    return peaks
end


"""
    get_peaks_indices(array)

return the indices of the signal peaks

# Arguments
- `array::Array{Float}`: the signal to find the turning points, or peaks
"""
function get_peaks_indices(array)
    A = array[:]
    # get rid of any zero slope in the beginning
    while A[2] == A[1]
        A = A[2:length(A)]
    end
    peaks = [1]
    if A[2] > A[1]
        slope = "p"
    elseif A[2] < A[1]
        slope = "m"
    end
    for i = 1:length(A)-2
        ind = i+1
        if slope == "p"
            if A[ind+1] < A[ind]
                peaks = append!(peaks,ind)
                slope = "m"
            end
        elseif slope == "m"
            if A[ind+1] > A[ind]
                peaks = append!(peaks,ind)
                slope = "p"
            end
        end
    end
    peaks = append!(peaks,length(A))
    return peaks
end


function calc_TI(constant,ai,TI_free,initial,sep,downstream)

    ti_calculation = constant * (1.0/3.0)^ai * TI_free^initial * sep^downstream

    TI = sqrt(ti_calculation^2 + TI_free^2)

    return TI
end


"""
    get_moments(out,Rhub,Rtip,r,az,precone,tilt)

Trapezoidal integration to find the blade root bending moment using the loads distribution

# Arguments
- `out::CCBlade dict: output from running CCBlade solve
- `Rhub::Float`: radius of the rotor hub
- `Rtip::Float`: radius of the blade tip
- `r::Array{Float}`: radial locations of interest
- `az::Float`: blade azimuth angle
- `precone::Float`: rotor precone angle
- `tilt::Float`: rotor tilt angle
"""
function get_moments(out,Rhub,Rtip,r,az,precone,tilt)
    nr = length(r)+2

    arr_type = promote_type(typeof(out[1].Np),typeof(out[1].Tp))
    loads_flap = zeros(arr_type,nr)
    loads_edge = zeros(arr_type,nr)
    r_arr = zeros(nr)

    counter = 2
    for i = 1:length(r)
        loads_flap[counter] = out[i].Np/1000.0
        loads_edge[counter] = out[i].Tp/1000.0
        counter += 1
    end


    # arr_type = promote_type(typeof(out.Np[1]),typeof(out.Tp[1]))
    # loads_flap = zeros(arr_type,nr)
    # loads_edge = zeros(arr_type,nr)
    # r_arr = zeros(nr)
    # loads_flap[2:end-1] = out.Np/1000.0
    # loads_edge[2:end-1] = out.Tp/1000.0

    # println("loads_flap: ", loads_flap[5:8])
    # println("loads_edge: ", loads_edge[5:8])

    r_arr[2:end-1] = r
    #
    #approximate loads at r = Rhub
    dL_flap = loads_flap[3]-loads_flap[2]
    dL_edge = loads_edge[3]-loads_edge[2]
    dr = r_arr[3]-r_arr[2]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[2] + m_flap*(Rhub-r_arr[2])
    Lhub_edge = loads_edge[2] + m_edge*(Rhub-r_arr[2])

    loads_flap[1] = Lhub_flap
    loads_edge[1] = Lhub_edge
    r_arr[1] = Rhub

    #approximate loads at r = Rtip
    dL_flap = loads_flap[end-1]-loads_flap[end-2]
    dL_edge = loads_edge[end-1]-loads_edge[end-2]
    dr = r_arr[end-1]-r_arr[end-2]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[end] + m_flap*(Rtip-r_arr[end-1])
    Lhub_edge = loads_edge[end] + m_edge*(Rtip-r_arr[end-1])

    loads_flap[end] = Lhub_flap
    loads_edge[end] = Lhub_edge
    r_arr[end] = Rtip

    M_flap = trapz(r_arr,loads_flap.*(r_arr.-Rhub))
    M_edge = trapz(r_arr,loads_edge.*(r_arr.-Rhub))

    # add gravity loads
    blade_mass=17536.617
    blade_cm=20.650
    grav=9.81
    M_edge += sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0
    return M_flap, M_edge
end


function get_moments_twist(out,Rhub,Rtip,r,az,precone,tilt,rotor,sections)

    pitch = rotor.pitch

    nr = length(r)+2

    arr_type = promote_type(typeof(out[1].Np),typeof(out[1].Tp))
    loads_hor = zeros(arr_type,nr)
    loads_vert = zeros(arr_type,nr)
    r_arr = zeros(nr)

    counter = 2
    for i = 1:length(r)
        twist = sections[counter-1].theta
        loads_hor[counter] = ((out[i].Np*cos(twist) + out[i].Tp*sin(twist))*cos(pitch) +
                            (out[i].Tp*cos(twist) + out[i].Np*sin(twist))*sin(pitch))/1000.0
        loads_vert[counter] = ((out[i].Tp*cos(twist) + out[i].Np*sin(twist))*cos(pitch) +
                            (out[i].Np*cos(twist) + out[i].Tp*sin(twist))*sin(pitch))/1000.0
        counter += 1
    end

    r_arr[2:end-1] = r
    #
    #approximate loads at r = Rhub
    dL_flap = loads_hor[3]-loads_hor[2]
    dL_edge = loads_vert[3]-loads_vert[2]
    dr = r_arr[3]-r_arr[2]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_hor[2] + m_flap*(Rhub-r_arr[2])
    Lhub_edge = loads_vert[2] + m_edge*(Rhub-r_arr[2])

    loads_hor[1] = Lhub_flap
    loads_vert[1] = Lhub_edge
    r_arr[1] = Rhub

    #approximate loads at r = Rtip
    dL_flap = loads_hor[end-1]-loads_hor[end-2]
    dL_edge = loads_vert[end-1]-loads_vert[end-2]
    dr = r_arr[end-1]-r_arr[end-2]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_hor[end] + m_flap*(Rtip-r_arr[end-1])
    Lhub_edge = loads_vert[end] + m_edge*(Rtip-r_arr[end-1])

    loads_hor[end] = Lhub_flap
    loads_vert[end] = Lhub_edge
    r_arr[end] = Rtip

    M_hor = trapz(r_arr,loads_hor.*(r_arr.-Rhub))
    M_vert = trapz(r_arr,loads_vert.*(r_arr.-Rhub))

    # add gravity loads
    blade_mass=17536.617
    blade_cm=20.650
    grav=9.81
    M_vert += sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0
    return M_hor, M_vert
end



# function get_single_damage(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,turbine_ID,state_ID,nCycles,az_arr,
#     turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set,
#     turbine_velocities, turbine_ct,turbine_local_ti;
#     Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2)
#
#
#         arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
#                             typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))
#
#         naz = length(az_arr)
#
#         nturbines = length(turbine_x)
#
#         naz = length(az_arr)
#         ws = wind_resource.wind_speeds[state_ID]
#         measurementheight = wind_resource.measurement_heights
#         air_density = wind_resource.air_density
#         wind_shear_model = wind_resource.wind_shear_model
#         ambient_tis = wind_resource.ambient_tis
#
#         flap = zeros(arr_type,nCycles*naz)
#         edge = zeros(arr_type,nCycles*naz)
#         oms = zeros(arr_type,nCycles*naz)
#
#         # turbine_velocities, turbine_ct, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
#         #                     turbine_ai, sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
#         #                     model_set; wind_farm_state_id=state_ID)
#
#         TI_inst_arr = zeros(arr_type,length(az_arr))
#         for i = 1:length(az_arr)
#             az = az_arr[i]
#
#             x_locs, y_locs, z_locs = ff.find_xyz_simple(turbine_x[turbine_ID],turbine_y[turbine_ID],turbine_z[turbine_ID].+hub_height[turbine_ID],r,turbine_yaw[turbine_ID],az)
#
#             TI_arr = zeros(arr_type,length(r))
#             for k = 1:length(r)
#                 loc = [x_locs[k],y_locs[k],z_locs[k]]
#                 TI_arr[k] = turbulence_func(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_tis[state_ID], div_sigma=div_sigma, div_ti=div_ti)
#             end
#
#             TI_inst = 0.0
#             for k = 1:length(r)
#                 if k == length(r)
#                     TI_inst += (Rtip-r[end])*TI_arr[end]/Rtip
#                 else
#                     TI_inst += (r[k+1]-r[k])*TI_arr[k]/Rtip
#                 end
#             end
#             TI_inst_arr[i] = TI_inst
#         end
#
#         U_arr = zeros(arr_type,(length(az_arr),length(r)))
#         for i = 1:length(az_arr)
#             az = az_arr[i]
#             U_arr[i,:] = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
#                                 sorted_turbine_index,turbine_velocities,wind_resource,model_set;wind_farm_state_id=state_ID,downwind_turbine_id=turbine_ID)
#         end
#
#         time1 = time()
#         for i = 1:nCycles*naz
#             az = az_arr[(i+1)%naz+1]
#
#             # x_locs, y_locs, z_locs = ff.find_xyz_simple(turbine_x[turbine_ID],turbine_y[turbine_ID],turbine_z[turbine_ID].+hub_height[turbine_ID],r,turbine_yaw[turbine_ID],az)
#             #
#             # TI_arr = zeros(length(r))
#             # for k = 1:length(r)
#             #     loc = [x_locs[k],y_locs[k],z_locs[k]]
#             #     TI_arr[k] = turbulence_func(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_tis[state_ID], div_sigma=div_sigma, div_ti=div_ti)
#             # end
#             #
#             # TI_inst = 0.0
#             # for k = 1:length(r)
#             #     if k == length(r)
#             #         TI_inst += (Rtip-r[end])*TI_arr[end]/Rtip
#             #     else
#             #         TI_inst += (r[k+1]-r[k])*TI_arr[k]/Rtip
#             #     end
#             # end
#             TI_inst = TI_inst_arr[(i+1)%naz+1]
#
#             turbine_velocity_with_turb = arr_type(turbine_velocities[turbine_ID] + turb_samples[i]*TI_inst*turbine_velocities[turbine_ID])
#             Omega_rpm = omega_func(turbine_velocity_with_turb) #in rpm
#             Omega = Omega_rpm*0.10471975512 #convert to rad/s
#             pitch_deg = pitch_func(turbine_velocity_with_turb) #in degrees
#             pitch = pitch_deg*pi/180.0 #convert to rad
#             # U = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
#             #                     sorted_turbine_index,turbine_velocities,wind_resource,model_set;wind_farm_state_id=state_ID,downwind_turbine_id=turbine_ID)
#             U = U_arr[(i+1)%naz+1,:]
#             U = U + turb_samples[i]*TI_inst.*U
#
#             rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
#             op = ff.distributed_velocity_op.(U, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)
#             out = CCBlade.solve.(Ref(rotor), sections, op)
#
#             flap[i],edge[i] = ff.get_moments(out,Rhub,Rtip,r,az,precone,tilt)
#             oms[i] = Omega
#         end
#         # println("loop 1: ", time()-time1)
#
#         # println("flap: ", flap)
#         # println("edge: ", edge)
#         xlocs = zeros(Nlocs)
#         ylocs = zeros(Nlocs)
#         angles = range(-pi/2.0,stop=pi/2.0,length=Nlocs)
#         root_rad=3.542/2.0
#         for i in 1:Nlocs
#             xlocs[i] = cos(angles[i])*root_rad
#             ylocs[i] = sin(angles[i])*root_rad
#         end
#
#         avg_omega = sum(oms)/length(oms)
#         total_time = nCycles/(avg_omega/(2.0*pi))
#
#         d_arr = zeros(arr_type,Nlocs)
#         damage = 0.0
#         ind = 0
#         time2 = time()
#         for i  = 1:Nlocs
#             sigma = calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
#             peaks = get_peaks(sigma)
#             out = rainflow(peaks)
#
#             alternate = out[1,:]/2.
#             mean = out[2,:]
#             count = out[3,:]
#
#             su = 70000.0
#             m = 10.0
#             years = 25.0
#             freq = 1.0
#
#             mar = alternate./(1.0.-mean./su)
#
#             npts = length(mar)
#             #damage calculations
#             d = 0.0
#             for i = 1:npts
#                 Nfail = ((su)/(mar[i]*fos))^m
#                 mult = years*365.25*24.0*3600.0*freq/total_time
#                 d += count[i]*mult/Nfail
#             end
#             d_arr[i] = d
#             if d > damage
#                 damage = d
#                 ind = i
#             end
#         end
#         # println("loop 2: ", time()-time2)
#         return damage
# end
#
#
#
#
# function get_single_state_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
#                 state_ID,nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,
#                 wind_resource,model_set;
#                 Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2)
#
#     arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
#                         typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))
#
#     nturbines = length(turbine_x)
#     damage = zeros(arr_type,nturbines)
#
#     turbine_velocities, turbine_ct, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
#                         turbine_ai, sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
#                         model_set; wind_farm_state_id=state_ID)
#     for k = 1:nturbines
#         damage[k] = ff.get_single_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,k,state_ID,nCycles,az_arr,
#             turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set,
#             turbine_velocities,turbine_ct,turbine_local_ti;
#             Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti)
#     end
#
#     return damage
# end
#
#
# function get_total_farm_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,ct_model,
#         nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set;
#         Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2)
#
#     arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
#                         typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))
#
#     nturbines = length(turbine_x)
#
#     frequencies = wind_resource.wind_probabilities
#     directions = wind_resource.wind_directions
#     ndirections = length(frequencies)
#
#     damage = zeros(arr_type,nturbines)
#     for k = 1:ndirections
#         rot_x, rot_y =  ff.rotate_to_wind_direction(turbine_x, turbine_y, directions[k])
#         sorted_turbine_index = sortperm(rot_x)
#         state_damage = ff.get_single_state_damage_surr(rot_x,rot_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
#                         k,nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,
#                         wind_resource,model_set,
#                         Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti)
#         damage = damage + state_damage.*frequencies[k]
#     end
#
#     return damage
# end



struct NpTp{T}
    Np::T
    Tp::T
end


function get_single_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,turbine_ID,state_ID,nCycles,az_arr,
    turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set,
    turbine_velocities, turbine_ct,turbine_local_ti;
    Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2,B=3)


        arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

        naz = length(az_arr)

        nturbines = length(turbine_x)

        naz = length(az_arr)
        ws = wind_resource.wind_speeds[state_ID]
        measurementheight = wind_resource.measurement_heights
        air_density = wind_resource.air_density
        wind_shear_model = wind_resource.wind_shear_model
        ambient_tis = wind_resource.ambient_tis

        flap = zeros(arr_type,nCycles*naz)
        edge = zeros(arr_type,nCycles*naz)
        oms = zeros(arr_type,nCycles*naz)

        TI_inst_arr = zeros(arr_type,length(az_arr))
        for i = 1:length(az_arr)
            az = az_arr[i]

            x_locs, y_locs, z_locs = ff.find_xyz_simple(turbine_x[turbine_ID],turbine_y[turbine_ID],turbine_z[turbine_ID].+hub_height[turbine_ID],r,turbine_yaw[turbine_ID],az)

            TI_arr = zeros(arr_type,length(r))
            for k = 1:length(r)
                loc = [x_locs[k],y_locs[k],z_locs[k]]
                TI_arr[k] = turbulence_func(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_tis[state_ID], div_sigma=div_sigma, div_ti=div_ti)
            end

            TI_inst = 0.0
            for k = 1:length(r)
                if k == length(r)
                    TI_inst += (Rtip-r[end])*TI_arr[end]/Rtip
                else
                    TI_inst += (r[k+1]-r[k])*TI_arr[k]/Rtip
                end
            end
            TI_inst_arr[i] = TI_inst
        end

        U_arr = zeros(arr_type,(length(az_arr),length(r)))
        for i = 1:length(az_arr)
            az = az_arr[i]
            U_arr[i,:] = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
                                sorted_turbine_index,turbine_velocities,wind_resource,model_set;wind_farm_state_id=state_ID,downwind_turbine_id=turbine_ID)
        end

        time1 = time()
        for i = 1:nCycles*naz
            az = az_arr[(i+1)%naz+1]

            TI_inst = TI_inst_arr[(i+1)%naz+1]

            turbine_velocity_with_turb = arr_type(turbine_velocities[turbine_ID] + turb_samples[i]*TI_inst*turbine_velocities[turbine_ID])
            Omega_rpm = omega_func(turbine_velocity_with_turb) #in rpm
            Omega = Omega_rpm*0.10471975512 #convert to rad/s
            pitch_deg = pitch_func(turbine_velocity_with_turb) #in degrees
            pitch = pitch_deg*pi/180.0 #convert to rad

            U = U_arr[(i+1)%naz+1,:]
            U = U + turb_samples[i]*TI_inst.*U

            # println("U: ", U)
            # println("turbine_velocity_with_turb: ", turbine_velocity_with_turb)
            # np = flap_func(U,turbine_velocity_with_turb)
            # tp = edge_func(U,turbine_velocity_with_turb)
            # # # tp = tp .* cos(az-3.1415926535897/2.0)
            # out = NpTp(np,tp)

            rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
            op = ff.distributed_velocity_op.(U, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)
            out = CCBlade.solve.(Ref(rotor), sections, op)

            flap[i],edge[i] = ff.get_moments(out,Rhub,Rtip,r,az,precone,tilt)
            oms[i] = Omega
        end

        xlocs = zeros(Nlocs)
        ylocs = zeros(Nlocs)
        angles = range(-pi/2.0,stop=pi/2.0,length=Nlocs)
        root_rad=3.542/2.0
        for i in 1:Nlocs
            xlocs[i] = cos(angles[i])*root_rad
            ylocs[i] = sin(angles[i])*root_rad
        end

        avg_omega = sum(oms)/length(oms)
        total_time = nCycles/(avg_omega/(2.0*pi))

        d_arr = zeros(arr_type,Nlocs)
        damage = 0.0
        ind = 0
        time2 = time()
        for i  = 1:Nlocs
            sigma = calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
            peaks = get_peaks(sigma)
            out = rainflow(peaks)

            alternate = out[1,:]/2.
            mean = out[2,:]
            count = out[3,:]

            su = 70000.0
            m = 10.0
            years = 25.0
            freq = 1.0

            mar = alternate./(1.0.-mean./su)

            npts = length(mar)
            #damage calculations
            d = 0.0
            for i = 1:npts
                Nfail = ((su)/(mar[i]*fos))^m
                mult = years*365.25*24.0*3600.0*freq/total_time
                d += count[i]*mult/Nfail
            end
            d_arr[i] = d
            if d > damage
                damage = d
                ind = i
            end
        end
        return damage
end




function get_single_state_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
                state_ID,nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,
                wind_resource,model_set;
                Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                        typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

    nturbines = length(turbine_x)
    damage = zeros(arr_type,nturbines)

    turbine_velocities, turbine_ct, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                        turbine_ai, sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                        model_set; wind_farm_state_id=state_ID)
    for k = 1:nturbines
        damage[k] = ff.get_single_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,k,state_ID,nCycles,az_arr,
            turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set,
            turbine_velocities,turbine_ct,turbine_local_ti;
            Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti)
    end

    return damage
end


function get_total_farm_damage_surr(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,ct_model,
        nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,wind_resource,model_set;
        Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                        typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

    nturbines = length(turbine_x)

    frequencies = wind_resource.wind_probabilities
    directions = wind_resource.wind_directions
    ndirections = length(frequencies)

    damage = zeros(arr_type,nturbines)
    for k = 1:ndirections
        rot_x, rot_y =  ff.rotate_to_wind_direction(turbine_x, turbine_y, directions[k])
        sorted_turbine_index = sortperm(rot_x)
        state_damage = ff.get_single_state_damage_surr(rot_x,rot_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
                        k,nCycles,az_arr,turb_samples,flap_func,edge_func,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,precone,tilt,rho,
                        wind_resource,model_set,
                        Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti)
        damage = damage + state_damage.*frequencies[k]
    end

    return damage
end




function get_single_damage_super(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,turbine_ID,state_ID,nCycles,az_arr,
    turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,wind_resource,model_set,
    turbine_velocities, turbine_ct,turbine_local_ti;
    Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2,B=3)


        arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

        naz = length(az_arr)

        nturbines = length(turbine_x)

        naz = length(az_arr)
        ws = wind_resource.wind_speeds[state_ID]
        measurementheight = wind_resource.measurement_heights
        air_density = wind_resource.air_density
        wind_shear_model = wind_resource.wind_shear_model
        ambient_tis = wind_resource.ambient_tis

        flap = zeros(arr_type,nCycles*naz)
        edge = zeros(arr_type,nCycles*naz)
        oms = zeros(arr_type,nCycles*naz)

        TI_inst_arr = zeros(arr_type,length(az_arr))
        for i = 1:length(az_arr)
            az = az_arr[i]

            x_locs, y_locs, z_locs = ff.find_xyz_simple(turbine_x[turbine_ID],turbine_y[turbine_ID],turbine_z[turbine_ID].+hub_height[turbine_ID],r,turbine_yaw[turbine_ID],az)

            TI_arr = zeros(arr_type,length(r))
            for k = 1:length(r)
                loc = [x_locs[k],y_locs[k],z_locs[k]]
                TI_arr[k] = turbulence_func(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_tis[state_ID], div_sigma=div_sigma, div_ti=div_ti)
            end

            TI_inst = 0.0
            for k = 1:length(r)
                if k == length(r)
                    TI_inst += (Rtip-r[end])*TI_arr[end]/Rtip
                else
                    TI_inst += (r[k+1]-r[k])*TI_arr[k]/Rtip
                end
            end
            TI_inst_arr[i] = TI_inst
        end

        # loc_arr_y = [0.69,0.69,-0.69,-0.69]
        # loc_arr_z = [0.69,-0.69,-0.69,0.69]
        # eff_TI = 0.0
        # for i = 1:4
        #     loc = [turbine_x[turbine_ID],turbine_y[turbine_ID]+loc_arr_y[i]*rotor_diameter[turbine_ID],turbine_z[turbine_ID]+hub_height[turbine_ID]+loc_arr_z[i]*rotor_diameter[turbine_ID]]
        #     eff_TI += turbulence_func(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_tis[state_ID], div_sigma=div_sigma, div_ti=div_ti)^2
        # end
        # eff_TI = sqrt(eff_TI)/4.0 + 0.0
        eff_TI = sum(TI_inst_arr)/length(TI_inst_arr)

        # Omega_rpm = omega_func(turbine_velocities[turbine_ID]) #in rpm
        # Omega_rpm = 22.0
        # Omega = Omega_rpm*0.10471975512 #convert to rad/s
        Omega = 7.571322070955497*turbine_velocities[turbine_ID]/(rotor_diameter[turbine_ID]/2.0)
        if Omega > 1.267109036952
            Omega = 1.267109036952
        end
        pitch_deg = pitch_func(turbine_velocities[turbine_ID]) #in degrees
        # pitch_deg = pitch_deg/1.5
        # println("pitch_deg: ", pitch_deg)
        # pitch_deg = 1.979211195533589
        # pitch_deg = 5.0
        # pitch_deg = 0.0
        pitch = pitch_deg*3.1415926535897/180.0 #convert to rad

        m_arr_flap = zeros(arr_type,length(az_arr))
        m_arr_edge = zeros(arr_type,length(az_arr))
        for i = 1:length(az_arr)
            az = az_arr[i]
            U = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
                                sorted_turbine_index,turbine_velocities,wind_resource,model_set;wind_farm_state_id=state_ID,downwind_turbine_id=turbine_ID)

            rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
            op = ff.distributed_velocity_op.(U, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)

            # V = zeros(length(U)).+10.0*cos(az)
            # W = zeros(length(U)).-10.0*sin(az)
            # op = ff.multiple_components_op.(U, V, W, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)
            out = CCBlade.solve.(Ref(rotor), sections, op)

            """testing the correct theta"""
            # m_arr_flap[i],m_arr_edge[i] = ff.get_moments(out,Rhub,Rtip,r,az,precone,tilt)
            m_arr_flap[i],m_arr_edge[i] = ff.get_moments_twist(out,Rhub,Rtip,r,az,precone,tilt,rotor,sections)
        end

        m_arr_flap_0 = zeros(arr_type,length(az_arr))
        m_arr_edge_0 = zeros(arr_type,length(az_arr))
        wr = ff.DiscretizedWindResource(wind_resource.wind_directions, [11.4], wind_resource.wind_probabilities, wind_resource.measurement_heights, wind_resource.air_density, wind_resource.ambient_tis, wind_resource.wind_shear_model)
        for i = 1:length(az_arr)
            az = az_arr[i]
            U = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
                                sorted_turbine_index,turbine_velocities,wr,model_set;wind_farm_state_id=1,downwind_turbine_id=turbine_ID)

            pitch = 0.0
            rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
            op = ff.distributed_velocity_op.(U, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)
            out = CCBlade.solve.(Ref(rotor), sections, op)

            m_arr_flap_0[i],m_arr_edge_0[i] = ff.get_moments(out,Rhub,Rtip,r,az,precone,tilt)
        end


        m_arr_flap_P = zeros(arr_type,length(az_arr))
        m_arr_edge_P = zeros(arr_type,length(az_arr))
        wr = ff.DiscretizedWindResource(wind_resource.wind_directions, [20.0], wind_resource.wind_probabilities, wind_resource.measurement_heights, wind_resource.air_density, wind_resource.ambient_tis, wind_resource.wind_shear_model)
        for i = 1:length(az_arr)
            az = az_arr[i]
            # U = get_speeds(turbine_x,turbine_y,turbine_z,turbine_ID,hub_height,r,turbine_yaw,az,turbine_ct,turbine_ai,rotor_diameter,turbine_local_ti,
            #                     sorted_turbine_index,turbine_velocities,wr,model_set;wind_farm_state_id=1,downwind_turbine_id=turbine_ID)
            U = ones(length(r)) .* 20.0
            pitch_deg = pitch_func(20.0) #in degrees
            pitch = pitch_deg*3.1415926535897/180.0 #convert to rad
            rotor = CCBlade.Rotor(Rhub, Rtip, B, true, pitch, precone)
            op = ff.distributed_velocity_op.(U, Omega, r, precone, turbine_yaw[turbine_ID], tilt, az, air_density)
            out = CCBlade.solve.(Ref(rotor), sections, op)

            m_arr_flap_P[i],m_arr_edge_P[i] = ff.get_moments(out,Rhub,Rtip,r,az,precone,tilt)
        end



        switch_num = 0
        for i = 1:nCycles*naz
            az = az_arr[(i+1)%naz+1]

            TI_inst = TI_inst_arr[(i+1)%naz+1]

            # turbine_velocity_with_turb = arr_type(turbine_velocities[turbine_ID] + turb_samples[i]*TI_inst*turbine_velocities[turbine_ID])
            turbine_velocity_with_turb = arr_type(turbine_velocities[turbine_ID] + turb_samples[i]*eff_TI*turbine_velocities[turbine_ID])
            # Omega_rpm = omega_func(turbine_velocity_with_turb) #in rpm
            # Omega_rpm = 22.0
            # Omega = Omega_rpm*0.10471975512 #convert to rad/s
            Omega = 7.571322070955497*turbine_velocity_with_turb/(rotor_diameter[turbine_ID]/2.0)
            if Omega > 1.267109036952
                Omega = 1.267109036952
            end

            # if turbine_velocity_with_turb > 11.4
            # flap[i] = m_arr_flap[(i+1)%naz+1] * (1.0+turb_samples[i]*TI_inst)^2
            # edge[i] = m_arr_edge[(i+1)%naz+1] * (1.0+turb_samples[i]*TI_inst)

            blade_mass=17536.617
            blade_cm=20.650
            grav=9.81
            # aero_edge = m_arr_edge[(i+1)%naz+1] - sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0

            if turbine_velocity_with_turb < 11.4
                flap[i] = m_arr_flap[(i+1)%naz+1] * (1.0+turb_samples[i]*TI_inst)^2
                # aero_edge = aero_edge * (1.0+turb_samples[i]*TI_inst)^2
                aero_edge = m_arr_edge_0[(i+1)%naz+1]
                edge[i] = aero_edge
            else
                # flap[i] = m_arr_flap[(i+1)%naz+1] * (1.0+turb_samples[i]*TI_inst)^-1
                # aero_edge = aero_edge * (1.0+turb_samples[i]*TI_inst)^-1
                aero_edge = m_arr_edge[(i+1)%naz+1] - sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0
                slope_edge = -604.12
                aero_edge = aero_edge + (turb_samples[i]*TI_inst)*slope_edge*turbine_velocities[turbine_ID]
                edge[i] = aero_edge + sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0

                slope_flap = -1781.08
                flap[i] = m_arr_flap[(i+1)%naz+1] + (turb_samples[i]*TI_inst)*slope_flap*turbine_velocities[turbine_ID]
            end
            # edge[i] = aero_edge + sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0



                # flap[i] = m_arr_flap[(i+1)%naz+1]
                # edge[i] = m_arr_edge[(i+1)%naz+1]

            # else
            #     flap[i] = m_arr_flap_0[(i+1)%naz+1] * (turbine_velocity_with_turb/11.4)
            #     # edge[i] = m_arr_edge[(i+1)%naz+1] * (1.0+turb_samples[i]*TI_inst)
            #
            #     blade_mass=17536.617
            #     blade_cm=20.650
            #     grav=9.81
            #     aero_edge = m_arr_edge_0[(i+1)%naz+1] - sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0
            #     aero_edge = aero_edge * (1.0+turb_samples[i]*TI_inst)
            #     edge[i] = aero_edge + sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.0
            #     # flap[i] = m_arr_flap[(i+1)%naz+1]
            #     # edge[i] = m_arr_edge[(i+1)%naz+1]
            # end

            oms[i] = Omega

            if i == 1
                global sign
                if turbine_velocity_with_turb > 11.4
                    sign = "above"
                else
                    sign = "below"
                end
            else
                global sign
                if sign == "above"
                    if turbine_velocity_with_turb < 11.4
                        sign = "below"
                        switch_num += 1
                    end
                elseif sign == "below"
                    if turbine_velocity_with_turb > 11.4
                        sign = "above"
                        switch_num += 1
                    end
                end
            end
        end

        # plot(edge)
        # plot(flap)
        # ylim(-4000,13000)

        # println("switch_num: ", switch_num)
        add_edge = zeros(switch_num*2)
        add_flap = zeros(switch_num*2)
        for i = 1:switch_num
            if i%2 == 1
                add_edge[i*2-1:i*2] = m_arr_edge_0[1:2]
                add_flap[i*2-1:i*2] = m_arr_flap_0[1:2]
            else
                add_edge[i*2-1:i*2] = m_arr_edge_P[1:2]
                add_flap[i*2-1:i*2] = m_arr_flap_P[1:2]
            end
        end

        # append!(edge,add_edge)
        # append!(flap,add_flap)

        xlocs = zeros(Nlocs)
        ylocs = zeros(Nlocs)
        angles = range(-3.1415926535897/2.0,stop=3.1415926535897/2.0,length=Nlocs)
        root_rad=3.542/2.0
        for i in 1:Nlocs
            xlocs[i] = cos(angles[i])*root_rad
            ylocs[i] = sin(angles[i])*root_rad
        end

        avg_omega = sum(oms)/length(oms)
        total_time = (nCycles)/(avg_omega/(2.0*3.1415926535897))

        d_arr = zeros(arr_type,Nlocs)
        damage = 0.0
        ind = 0
        time2 = time()
        for i  = 1:Nlocs
            sigma = calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
            # sigma = calc_moment_stress_resultant.(edge,flap,xlocs[i],ylocs[i])
            peaks = get_peaks(sigma)
            out = rainflow(peaks)

            alternate = out[1,:]/2.
            mean = out[2,:]
            count = out[3,:]

            su = 70000.0
            m = 10.0
            years = 25.0
            freq = 1.0

            mar = alternate./(1.0.-mean./su)

            npts = length(mar)
            #damage calculations
            d = 0.0
            for i = 1:npts
                Nfail = ((su)/(mar[i]*fos))^m
                mult = years*365.25*24.0*3600.0*freq/total_time
                d += count[i]*mult/Nfail
            end
            d_arr[i] = d
            if d > damage
                damage = d
                ind = i
            end
        end
        return damage
end



function get_single_state_damage_super(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
                state_ID,nCycles,az_arr,turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,
                wind_resource,model_set;
                Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2,B=3)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                        typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

    nturbines = length(turbine_x)
    damage = zeros(arr_type,nturbines)

    turbine_velocities, turbine_ct, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                        turbine_ai, sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                        model_set; wind_farm_state_id=state_ID)
    for k = 1:nturbines
        damage[k] = ff.get_single_damage_super(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,k,state_ID,nCycles,az_arr,
            turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,wind_resource,model_set,
            turbine_velocities,turbine_ct,turbine_local_ti;
            Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti,B=B)
    end

    return damage
end


function get_total_farm_damage_super(turbine_x,turbine_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,ct_model,
        nCycles,az_arr,turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,wind_resource,model_set;
        Nlocs=20,fos=1.15,rotor_sample_points_y=[0.69,0.69,-0.69,-0.69],rotor_sample_points_z=[0.69,-0.69,0.69,-0.69],div_sigma=2.5,div_ti=1.2,B=3)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                        typeof(hub_height[1]),typeof(turbine_yaw[1]),typeof(turbine_ai[1]))

    nturbines = length(turbine_x)

    frequencies = wind_resource.wind_probabilities
    directions = wind_resource.wind_directions
    ndirections = length(frequencies)

    damage = zeros(arr_type,nturbines)

    for k = 1:ndirections
        # t = time()
        rot_x, rot_y =  ff.rotate_to_wind_direction(turbine_x, turbine_y, directions[k])
        sorted_turbine_index = sortperm(rot_x)
        state_damage = ff.get_single_state_damage_super(rot_x,rot_y,turbine_z,rotor_diameter,hub_height,turbine_yaw,turbine_ai,sorted_turbine_index,ct_model,
                        k,nCycles,az_arr,turb_samples,omega_func,pitch_func,turbulence_func,r,sections,Rhub,Rtip,precone,tilt,
                        wind_resource,model_set,
                        Nlocs=Nlocs,fos=fos,rotor_sample_points_y=rotor_sample_points_y,rotor_sample_points_z=rotor_sample_points_z,div_sigma=div_sigma,div_ti=div_ti,B=B)
        damage = damage + state_damage.*frequencies[k]
        # println("state $k: ", time()-t)
    end

    return damage
end
