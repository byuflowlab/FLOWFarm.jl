using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using Random

const ff=FlowFarm


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


function calc_moment_stress(mx,my,dx,dy,Rcyl=1.771,tcyl=0.06)

        I = 0.25*pi*(Rcyl^4-(Rcyl-tcyl)^4)
        stress = mx*dy/I + my*dx/I

        return stress
end


function get_speeds(turbineX,turbineY,turb_index,hubHt,r,yaw,azimuth,model_set,problem_description)
        x_locs, y_locs, z_locs = find_xyz_simple(turbineX[turb_index],turbineY[turb_index],hubHt,r,yaw,azimuth)
        npts = length(x_locs)
        point_velocities = zeros(npts)
            for i in 1:npts
                loc = [x_locs[i], y_locs[i], z_locs[i]]
                # point_velocities[i] = FlowFarm.point_velocity(loc, windfarm, windfarmstate, windresource, wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, 0)
                point_velocities[i] = FlowFarm.point_velocity(loc, model_set, problem_description)
            end
        return point_velocities
end


function rainflow(array_ext;uc_mult=0.5)

    # """ Rainflow counting of a signal's turning points with Goodman correction
    #
    #     Args:
    #         array_ext (numpy.ndarray): array of turning points
    #
    #     Keyword Args:
    #         uc_mult (float): partial-load scaling [opt, default=0.5]
    #
    #     Returns:
    #         array_out (numpy.ndarray): (3 x n_cycle) array of rainflow values:
    #                                     1) load range
    #                                     2) range mean
    #                                     3) cycle count
    # """

    tot_num = length(array_ext)             # total size of input array
    array_out = zeros(3, tot_num-1)         # initialize output array

    pr = 1                                  # index of input array
    po = 1                                  # index of output array
    j = 0                                   # index of temporary array "a"
    a  = zeros(tot_num)                     # temporary array for algorithm

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


function get_moments(out,Rhub,Rtip,r,az)
    loads_flap = out.Np/1000.
    loads_edge = out.Tp/1000.
    #
    #approximate loads at r = Rhub
    dL_flap = loads_flap[2]-loads_flap[1]
    dL_edge = loads_edge[2]-loads_edge[1]
    dr = r[2]-r[1]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[1] + m_flap*(Rhub-r[1])
    Lhub_edge = loads_edge[1] + m_edge*(Rhub-r[1])

    loads_flap = append!([Lhub_flap],loads_flap)
    loads_edge = append!([Lhub_edge],loads_edge)
    r = append!([Rhub],r)

    #approximate loads at r = Rtip
    dL_flap = loads_flap[end]-loads_flap[end-1]
    dL_edge = loads_edge[end]-loads_edge[end-1]
    dr = r[end]-r[end-1]

    m_flap = dL_flap/dr
    m_edge = dL_edge/dr

    Lhub_flap = loads_flap[end] + m_flap*(Rtip-r[end])
    Lhub_edge = loads_edge[end] + m_edge*(Rtip-r[end])

    loads_flap = append!(loads_flap,[Lhub_flap])
    loads_edge = append!(loads_edge,[Lhub_edge])
    r = append!(r,[Rtip])

    M_flap = trapz(r,loads_flap.*(r.-Rhub))
    M_edge = trapz(r,loads_edge.*(r.-Rhub))

    # add gravity loads
    blade_mass=17536.617
    blade_cm=20.650
    grav=9.81
    M_edge += sin(az)*cos(precone)*cos(tilt)*blade_mass*grav*blade_cm/1000.
    return M_flap, M_edge
end


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


function get_single_damage(model_set,problem_description,turbine_ID,state_ID,nCycles,az_arr,
    turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip;
    Nlocs=20)

        flap = []
        edge = []
        oms = []
        naz = length(az_arr)

        turbine_x = problem_description.wind_farm_states[state_ID].turbine_x
        turbine_y = problem_description.wind_farm_states[state_ID].turbine_y
        nturbines = length(turbine_x)

        naz = length(az_arr)
        ws = problem_description.wind_resource.wind_speeds[state_ID]
        measurementheight = problem_description.wind_resource.measurement_heights
        air_density = problem_description.wind_resource.air_density
        wind_shear_model = problem_description.wind_resource.wind_shear_model
        ambient_tis = problem_description.wind_resource.ambient_tis
        windfarm = problem_description.wind_farm
        windfarmstate = problem_description.wind_farm_states[state_ID]
        turb_type = problem_description.wind_farm.turbine_definition_ids[turbine_ID]
        hub_height = problem_description.wind_farm.turbine_definitions[turb_type].hub_height
        yaw = problem_description.wind_farm_states[state_ID].turbine_yaw[turbine_ID]


        for i = 1:nCycles*naz
            az = az_arr[(i+1)%naz+1]

            """need to figure this out"""
            x_locs, y_locs, z_locs = find_xyz_simple(turbine_x[turbine_ID],turbine_y[turbine_ID],turbine_z[turbine_ID].+hub_height,r,yaw,az)
            TI_arr = zeros(length(r))
            for k = 1:length(r)
                loc = [x_locs[k],y_locs[k],z_locs[k]]
                TI_arr[k] = turbulence_func(loc)
            end
            TI_inst = sum(TI_arr)/length(TI_arr)
            # TI_inst = 0.16
            """"""

            windspeeds = ws + turb_samples[i]*TI_inst*ws
            turbine_inflow_velcities = zeros(nturbines) .+ windspeeds
            temp_resource = ff.DiscretizedWindResource([3*pi/2], [windspeeds], [1.0], measurementheight, air_density,ambient_tis, wind_shear_model)
            temp_pd = ff.WindFarmProblemDescription(windfarm, temp_resource, [windfarmstate])

            ff.turbine_velocities_one_direction!(points_x, points_y, model_set, temp_pd)
            turbine_inflow_velcities = temp_pd.wind_farm_states[state_ID].turbine_inflow_velcities
            Omega_rpm = omega_func(turbine_inflow_velcities[turb_index]) #in rpm
            Omega = Omega_rpm*0.10471975512 #convert to rad/s
            pitch_deg = pitch_func(turbine_inflow_velcities[turb_index]) #in degrees
            pitch = pitch_deg*pi/180.0 #convert to rad
            U = get_speeds(turbine_x,turbine_y,turbine_ID,hub_height,r,yaw,az,model_set,temp_pd)
            op = distributed_velocity_op.(U, Omega, r, precone, yaw, tilt, az, rho)
            out = CCBlade.solve.(Ref(rotor), sections, op)
            # out = CCBlade.solve(rotor, sections[1], op[1])
            flap_temp,edge_temp = get_moments(out,Rhub,Rtip,r,az)
            flap = append!(flap,flap_temp)
            edge = append!(edge,edge_temp)
            oms = append!(oms,Omega)

        end

        xlocs = zeros(Nlocs)
        ylocs = zeros(Nlocs)
        angles = range(-pi/2.,stop=pi/2.,length=Nlocs)
        root_rad=3.542/2.
        for i in 1:Nlocs
            xlocs[i] = cos(angles[i])*root_rad
            ylocs[i] = sin(angles[i])*root_rad
        end

        avg_omega = sum(oms)/length(oms)
        total_time = nCycles/(avg_omega/(2.0*pi))

        damage = 0.0
        for i  = 1:Nlocs
            sigma = calc_moment_stress.(edge,flap,xlocs[i],ylocs[i])
            peaks = get_peaks(sigma)
            out = rainflow(peaks)

            alternate = out[1,:]/2.
            mean = out[2,:]
            count = out[3,:]

            su = 70000.
            m = 10.
            years = 25.
            freq = 1.0

            mar = alternate./(1.0.-mean./su)

            npts = length(mar)
            #damage calculations
            d = 0.0
            d_arr = []
            for i = 1:npts
                Nfail = ((su)/(mar[i]*fos))^m
                mult = years*365.25*24.0*3600.0*freq/total_time
                d += count[i]*mult/Nfail
                d_arr = append!(d_arr,count[i]*mult/Nfail)
            end
            if d > damage
                damage = d
            end
        end

        return damage
end


function get_single_state_damage(model_set,problem_description,state_ID,nCycles,az_arr,
    turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip;
    Nlocs=20)

    turbine_x = problem_description.wind_farm_states[state_ID].turbine_x
    nturbines = length(turbine_x)
    damage = zeros(nturbines)
    for k = 1:nturbines
        damage[k] = get_single_damage(model_set,problem_description,k,state_ID,nCycles,az_arr,
            turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,
            Nlocs=Nlocs)
    end

    return damage
end


function get_total_farm_damage(model_set,problem_description,nCycles,az_arr,
    turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip;
    Nlocs=20)

    turbine_x = problem_description.wind_farm_states[1].turbine_x
    nturbines = length(turbine_x)

    frequencies = problem_description.wind_resource.wind_probabilities
    ndirections = length(frequencies)

    damage = zeros(nturbines)
    for k = 1:ndirections

        problem_description.wind_farm_states[i].turbine_x[:],problem_description.wind_farm_states[i].turbine_y[:] =
                rotate_to_wind_direction(wind_farm.turbine_x, wind_farm.turbine_y, wind_resource.wind_directions[i])

        problem_description.wind_farm_states[i].sorted_turbine_index[:] = sortperm(problem_description.wind_farm_states[i].turbine_x)

        state_damage = get_single_state_damage(model_set,problem_description,k,nCycles,az_arr,
            turb_samples,points_x,points_y,omega_func,pitch_func,turbulence_func,r,rotor,sections,Rhub,Rtip,
            Nlocs=Nlocs)
        damage = damage + state_damage.*frequencies[k]
    end

    return damage
end


for i = 1:nstates
    problem_description.wind_farm_states[i].turbine_x[:],problem_description.wind_farm_states[i].turbine_y[:] =
            rotate_to_wind_direction(wind_farm.turbine_x, wind_farm.turbine_y, wind_resource.wind_directions[i])

    problem_description.wind_farm_states[i].sorted_turbine_index[:] = sortperm(problem_description.wind_farm_states[i].turbine_x)

    turbine_velocities_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
        model_set, problem_description, wind_farm_state_id=i)

    turbine_powers_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
        problem_description, wind_farm_state_id=i)

    state_power = sum(problem_description.wind_farm_states[i].turbine_generators_powers)
    state_energy[i] = state_power*seconds_per_year*wind_probabilities[i]
end
