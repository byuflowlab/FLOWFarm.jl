abstract type AbstractModelSet end
# using CSV
# using DataFrames

"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, local_ti_model)

Container for objects defining models to use in wind farm calculations

# Arguments
- `wake_defiict_model::AbstractWakeDeficitModel`: contains a struct defining the desired wake deficit model
- `wake_deflection_model::AbstractWakeDeflectionModel`: contains a struct defining the desired wake deflection model
- `wake_combination_model::AbstractWakeCombinationModel`: contains a struct defining the desired wake combination model
- `local_ti_model::AbstractTurbulenceIntensityModel`: contains a struct defining the desired turbulence intensity model
"""
struct WindFarmModelSet{DTM,DNM,CM,TIM} <: AbstractModelSet

    wake_deficit_model::DTM
    wake_deflection_model::DNM
    wake_combination_model::CM
    local_ti_model::TIM

end

"""
    point_velocity(loc, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the wind speed at a given point for a given state

# Arguments
- `loc::Array{TF,3}`: Location of interest
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: turbine thrust coefficients for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: array containing indices of wind turbines 
    from most upwind to most downwind turbine in the given state
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: contains wind resource discreption (directions,
    speeds, frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `downwind_turbine_id::Int`: index of wind turbine of interest (if any). If not a point for
    calculating effective wind speed of a turbine, then provide 0 (default)
"""
function point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set::AbstractModelSet;
                    wind_farm_state_id=1, downwind_turbine_id=0)

    wakedeficitmodel = model_set.wake_deficit_model
    wakedeflectionmodel = model_set.wake_deflection_model
    wakecombinationmodel = model_set.wake_combination_model
    
    # extract flow information
    wind_speed = wind_resource.wind_speeds[wind_farm_state_id]
    reference_height = wind_resource.measurement_heights[wind_farm_state_id]

    # set ground height 
    ground_height = wind_resource.wind_shear_model.ground_height    # TODO: allow topology to be given

    # find order for wind shear and deficit calculations
    shear_order = wind_resource.wind_shear_model.shear_order

    # adjust wind speed for wind shear
    if shear_order == "nothing"
        wind_speed_internal = wind_speed
    elseif shear_order == "first"
        wind_speed_internal = adjust_for_wind_shear(locz, wind_speed, reference_height, ground_height, wind_resource.wind_shear_model)
    else
        wind_speed_internal = wind_speed
    end

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = Int(sorted_turbine_index[u])

        # don't allow turbine to impact itself
        if upwind_turb_id == downwind_turbine_id; continue; end

        # downstream distance between upstream turbine and point
        x = locx - turbine_x[upwind_turb_id]

        # check turbine relative locations
        if x > 1E-6
            # skip this loop if it would include a turbine's impact on itself)
            if upwind_turb_id==downwind_turbine_id; continue; end

            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct,
                            upwind_turb_id, rotor_diameter, turbine_local_ti, wakedeflectionmodel)

            vertical_deflection = 0.0

            # velocity difference in the wake
            deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                            upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                            turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel)

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, wakecombinationmodel)
            # println(deficit_sum, " ", downwind_turbine_id, " ", upwind_turb_id)
            
        end
    end

    # find velocity at point without shear
    point_velocity = wind_speed_internal - deficit_sum

    if shear_order == "nothing"
        point_velocity_out = point_velocity
    elseif shear_order == "first"
        point_velocity_out = point_velocity
    else
        point_velocity_out = adjust_for_wind_shear(locz, point_velocity, reference_height, ground_height, wind_resource.wind_shear_model)        
    end

    return point_velocity_out

end

"""
    calculate_transverse_velocity(U_i, U_inf, dx, dy, z, rotor_diameter, HH, tilt, ct, TSR, turbine_ai, shearexponent)

Calculates the vertical and horizontal spanwise velocities induced by upstream turbines on
the downstream turbine being looked at

# Arguments
- `U_inf::Float`: streamwise wind velocity at turbine being compared to
- `W::Float`: vertical spanwise wind velocity at turbine being compared to
- `U_inf_initial::Float`: streamwise wind velocity of wind farm
- `deltay::Float`: distance in y-direction between turbine of interest and turbine being compared to
- `z_i::Float`: vertical location of center of turbine being compared to
- `rotor_diameter::Float`: rotor diameter of turbine being comapared to
- `hub_height::Float`: hub height of turbine being compared to
- `cT::Float`: coefficient of thrust for turbine being compared to
- `axial_induction::Float`: axial induction factor for turbine being compared to
- `shearexponent::Float`: exponent used for determing ambient wind speed at different heights
"""

function calculate_transverse_velocity(U_i, U_inf, dx, dy, z, rotor_diameter, HH, tilt, cT, TSR, turbine_ai, shearexponent)

    # turbine parameters
    D = rotor_diameter
    ai = turbine_ai

    # flow parameters
    Uinf = U_inf

    # epsilon gain
    eps_gain = 0.2
    eps = eps_gain * D

    # Setting yaw to zero, still need to create way to combine tilt and yaw
    yaw = 0

    # find velocity at top and bottom of rotor swept area
    vel_top = ((HH+D/2)/HH)^shearexponent
    vel_bottom = ((HH-D/2)/HH)^shearexponent

    # find Gamma at the top and bottom of the rotor swept area
    # Gamma top and bottom may only be for yaw, tilt is gamma left and right
    # Gamma_top = gamma(D, vel_top, Uinf, cT)
    Gamma_top = sin(yaw)*cos(yaw)*(pi/8)*D*Uinf*cT*vel_top                # why do we use Uinf here, should we use U_inf?
    # Gamma_bot = -1*gamma(D, vel_bottom, Uinf, cT)
    Gamma_bot = sin(yaw)*cos(yaw)*(pi/8)*D*Uinf*cT*vel_bottom*-1.0            # does Ct already contain cos(tilt)? yes


    # Does this need to be multiplied by Uinf????
    Gamma_left = sin(tilt)*cos(tilt)*(pi/8)*D*Uinf*cT                       # why do we use Uinf here, should we use U_inf?
    # Gamma_bot = -1*gamma(D, vel_bottom, Uinf, cT)
    Gamma_right = sin(tilt)*cos(tilt)*(pi/8)*D*Uinf*cT*-1.0                   # does Ct already contain cos(tilt)? yes

    turbine_average_velocity = U_i

    Gamma_wr = 2*0.25*pi*D*(ai-ai^2)*turbine_average_velocity/TSR

    ### compute the spanwise and vertical velocities induced by tilt

    # decay the vortices as they move downstream using mixing length
    lambda = D/8.0
    kappa = 0.41
    lm = (kappa*z)/(1+(kappa*z/lambda))
    ### TODO: assign dudz_initial as an input, how to calculate this????  ###
    dudz = abs(Uinf*(vel_top-vel_bottom))/(D)
    turbulent_visc = (lm^2)*abs(dudz)

    decay = (eps^2)/((4*turbulent_visc*dx/Uinf) + eps^2)

    # top vortex
    dz_top = z - (HH+D/2)             
    r_top = dy^2 + dz_top^2             
    core_shape = 1-exp(-r_top/(eps^2))          
    v_1 = ((Gamma_top*dz_top)/(2*pi*r_top))*core_shape*decay
    w_1 = ((-1*Gamma_top*dy)/(2*pi*r_top))*core_shape*decay

    # bottom vortex
    dz_bot = z - (HH - D/2)
    r_bot = dy^2 + dz_bot^2
    core_shape = 1-exp(-r_bot/(eps^2))
    v_2 = ((Gamma_bot*dz_bot)/(2*pi*r_bot))*core_shape*decay
    w_2 = ((-1*Gamma_bot*dy)/(2*pi*r_bot))*core_shape*decay

    # left vortex
    dz_left = z - HH
    dy_left = dy-(D/2)
    r_l = dy_left^2 + dz_left^2
    core_shape = 1-exp(-r_l/(eps^2))     
    v_3 = ((Gamma_left*dz_left)/(2*pi*r_l))*core_shape*decay
    w_3 = ((-1*Gamma_left*dy_left)/(2*pi*r_l))*core_shape*decay

    # right vortex
    dz_right = z - HH
    dy_right = dy+(D/2)
    r_r = dy_right^2 + dz_right^2
    core_shape = 1-exp(-r_r/(eps^2))     
    v_4 = ((Gamma_right*dz_right)/(2*pi*r_r))*core_shape*decay
    w_4 = ((-1*Gamma_right*dy_right)/(2*pi*r_r))*core_shape*decay

    # wake rotation vortex
    dz_center = z - HH
    r_center = dy^2 + dz_center^2
    core_shape = 1 - exp(-r_center/(eps^2))
    v_5 = ((Gamma_wr*dz_center)/(2*pi*r_center))*core_shape*decay
    w_5 = ((-1*Gamma_wr*dy)/(2*pi*r_center))*core_shape*decay

    ### Boundary condition - ground mirror vortex
    # top vortex
    dz_top = z + (HH+D/2)             
    r_top = dy^2 + dz_top^2             
    core_shape = 1-exp(-r_top/(eps^2))          
    v_6 = ((-1*Gamma_top*dz_top)/(2*pi*r_top))*core_shape*decay
    w_6 = ((Gamma_top*dy)/(2*pi*r_top))*core_shape*decay

    # bottom vortex
    dz_bot = z + (HH - D/2)
    r_bot = dy^2 + dz_bot^2
    core_shape = 1-exp(-r_bot/(eps^2))
    v_7 = ((-1*Gamma_bot*dz_bot)/(2*pi*r_bot))*core_shape*decay
    w_7 = ((Gamma_bot*dy)/(2*pi*r_bot))*core_shape*decay

    # wake rotation vortex
    dz_center = z + HH
    r_center = dy^2 + dz_center^2
    core_shape = 1 - exp(-r_center/(eps^2))
    v_8 = ((-1*Gamma_wr*dz_center)/(2*pi*r_center))*core_shape*decay
    w_8 = ((Gamma_wr*dy)/(2*pi*r_center))*core_shape*decay

    # left vortex
    dz_left = z + HH
    dy_left = dy-(D/2)
    r_l = dy_left^2 + dz_left^2
    core_shape = 1-exp(-r_l/(eps^2))     
    v_9 = ((-1*Gamma_left*dz_left)/(2*pi*r_l))*core_shape*decay
    w_9 = ((Gamma_left*dy_left)/(2*pi*r_l))*core_shape*decay

    # right vortex
    dz_right = z + HH
    dy_right = dy+(D/2)
    r_r = dy_right^2 + dz_right^2
    core_shape = 1-exp(-r_r/(eps^2))     
    v_10 = ((-1*Gamma_right*dz_right)/(2*pi*r_r))*core_shape*decay
    w_10 = ((Gamma_right*dy_right)/(2*pi*r_r))*core_shape*decay

    # total spanwise velocities
    V = v_1 + v_2 + v_3 + v_4 + v_5 + v_6 + v_7 + v_8 + v_9 + v_10
    W = w_1 + w_2 + w_3 + w_4 + w_5 + w_6 + w_7 + w_8 + w_9 + w_10

    if dx < 0.0
        V = 0
        W = 0
    end

    # Why can't W be negative?
    if W < 0
        W = 0 
    end

    return V, W

end

"""
    tilt_added_turbulence_intensity(u_i, w_i, I_i, v_i, turb_v_i, turb_w_i)

Calculates the added tilt due to secondary wake steering

# Arguments
- `u_i::Float`: 
- `w_i::Float`: 
- `I_i::Float`: 
- `v_i::Float`: 
- `z_i::Float`: 
- `turb_v_i::Float`:
- `turb_w_i::Float`: 
"""

function tilt_added_turbulence_intensity(u_i, w_i, I_i, v_i, turb_v_i, turb_w_i)
    # Convert ambient TI to TKE
    k = ((u_i*I_i)^2)/(2/3)
    u_term = sqrt(2*k)
    v_term = v_i + turb_v_i
    w_term = w_i + turb_w_i

    # Compute new TKE
    k_total = 0.5 * (u_term^2 + v_term^2 + w_term^2)

    # Convert TKE back to TI
    I_total = sqrt((2/3)*k_total)/u_i

    # solve for TI due to mixing
    I_mixing = I_total - I_i

    return I_mixing
end

"""
    brent(f, a, b; args=(), atol=2e-12, rtol=4*eps(), maxiter=100)
1D root finding using Brent's method.  Based off the brentq implementation in scipy.
**Arguments**
- `f`: scalar function, that optionally takes additional arguments
- `a`::Float, b::Float`: bracketing interval for a root - sign changes sign between: (f(a) * f(b) < 0)
- `args::Tuple`: tuple of additional arguments to pass to f
- `atol::Float`: absolute tolerance (positive) for root
- `rtol::Float`: relative tolerance for root
- `maxiter::Int`: maximum number of iterations allowed
**Returns**
- `xstar::Float`: a root of f
- `info::Tuple`: A named tuple containing:
    - `iter::Int`: number of iterations
    - 'fcalls::Int`: number of function calls
    - 'flag::String`: a convergence/error message.
"""
function brent(f, a, b; args=(), atol=2e-12, rtol=4*eps(), maxiter=100)

    xpre = a; xcur = b
    # xblk = 0.0; fblk = 0.0; spre = 0.0; scur = 0.0
    error_num = "INPROGRESS"

    fpre = f(xpre, args...)
    fcur = f(xcur, args...)

    xblk = zero(fpre); fblk = zero(fpre); spre = zero(fpre); scur = zero(fpre)
    funcalls = 2
    iterations = 0
    
    if fpre*fcur > 0
        error_num = "SIGNERR"
        return 0.0, (iter=iterations, fcalls=funcalls, flag=error_num)
    end
    if fpre == zero(fpre)
        error_num = "CONVERGED"
        return xpre, (iter=iterations, fcalls=funcalls, flag=error_num)
    end
    if fcur == zero(fcur)
        error_num = "CONVERGED"
        return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
    end

    for i = 1:maxiter
        iterations += 1
        if fpre*fcur < 0
            xblk = xpre
            fblk = fpre
            spre = scur = xcur - xpre
        end
        if abs(fblk) < abs(fcur)
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre
        end

        delta = (atol + rtol*abs(xcur))/2
        sbis = (xblk - xcur)/2
        if fcur == zero(fcur) || abs(sbis) < delta
            error_num = "CONVERGED"
            return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
        end

        if abs(spre) > delta && abs(fcur) < abs(fpre)
            if xpre == xblk
                # interpolate
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                # extrapolate
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre))
            end
            if 2*abs(stry) < min(abs(spre), 3*abs(sbis) - delta)
                # good short step
                spre = scur
                scur = stry
            else
                # bisect
                spre = sbis
                scur = sbis
            end
        else 
            # bisect
            spre = sbis
            scur = sbis
        end

        xpre = xcur; fpre = fcur
        if abs(scur) > delta
            xcur += scur
        else
            xcur += (sbis > 0 ? delta : -delta)
        end

        fcur = f(xcur, args...)
        funcalls += 1
    end
    error_num = "CONVERR"
    return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
end


"""
    tilt_eval(tilt, args...)

calculates the difference between W at the turbine and W induced by a tilt angle to find the added tilt due to secondary wake steering

# Arguements
- `tilt::Float`: tilt angle of turbine in radians
- `args::Tuple`: tuple of additional arguments to pass to W calculations
"""
function tilt_eval(tilt, args...)

    D = args[1]
    Uinf = args[2]
    cT = args[3]
    ai = args[4]
    TSR = args[5]
    z_i = args[6]
    HH = args[7]
    deltay = args[8]
    eps = args[9]
    avgW = args[10]
    U_inf = args[11]

    Gamma_left = (pi/8)*D*Uinf*cT*sin(tilt)*cos(tilt)
    Gamma_right = (pi/8)*D*Uinf*cT*sin(tilt)*cos(tilt)*-1

    # Use turbine average velocity to find Gamma due to wake rotation
    turbine_average_velocity = U_inf
    Gamma_wake_rotation = 0.25*2.0*pi*D*(ai-ai^2)*turbine_average_velocity/TSR

    # Left vortex
    dz_left = z_i- HH
    dy_left = deltay-(D/2)
    r_l = dy_left^2 + dz_left^2
    core_shape = 1-exp(-r_l/(eps^2))     
    w_left = ((-1*Gamma_left*dy_left)/(2*pi*r_l))*core_shape

    # right vortex
    dz_right = z_i - HH
    dy_right = deltay+(D/2)
    r_r = dy_right^2 + dz_right^2
    core_shape = 1-exp(-r_r/(eps^2))     
    w_right = ((-1*Gamma_right*dy_right)/(2*pi*r_r))*core_shape

    # wake rotation vortex
    dz_center = z_i - HH
    r_center = deltay^2 + dz_center^2
    core_shape = 1 - exp(-r_center/(eps^2))
    w_center = ((-1*Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape

    ### Boundary condition - ground mirror vortex
    # Left vortex
    dz_left = z_i + HH
    dy_left = deltay-(D/2)
    r_l = dy_left^2 + dz_left^2
    core_shape = 1-exp(-r_l/(eps^2))     
    w_left_g = ((Gamma_left*dy_left)/(2*pi*r_l))*core_shape

    # right vortex
    dz_right = z_i + HH
    dy_right = deltay+(D/2)
    r_r = dy_right^2 + dz_right^2
    core_shape = 1-exp(-r_r/(eps^2))     
    w_right_g = ((Gamma_right*dy_right)/(2*pi*r_r))*core_shape

    # wake rotation vortex
    dz_center = z_i + HH
    r_center = deltay^2 + dz_center^2
    core_shape = 1 - exp(-r_center/(eps^2))
    w_center_g = ((Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape

    W_eff = w_left+w_right+w_center + w_left_g+w_right_g+w_center_g
    # W_eff = w_left+w_right
    # W_eff = w_left+w_right + w_left_g+w_right_g
    
    residual = (avgW - 2*W_eff)
    # print("residual: ", residual, "\n")

    return residual
end


"""
    wake_added_tilt(U_inf, W, U_inf_initial, deltay, z_i, rotor_diameter, hub_height, 
    cT, TSR, axial_induction)

Calculates the added tilt due to secondary wake steering

# Arguments
- `U_inf::Float`: streamwise wind velocity at turbine being compared to
- `W::Float`: vertical spanwise wind velocity at turbine being compared to
- `U_inf_initial::Float`: streamwise wind velocity of wind farm
- `deltay::Float`: distance in y-direction between turbine of interest and turbine being compared to
- `z_i::Float`: vertical location of center of turbine being compared to
- `rotor_diameter::Float`: rotor diameter of turbine being comapared to
- `hub_height::Float`: hub height of turbine being compared to
- `cT::Float`: coefficient of thrust for turbine being compared to
- `axial_induction::Float`: axial induction factor for turbine being compared to
"""
function wake_added_tilt(U_inf, W, U_inf_initial, deltay, z_i, rotor_diameter, hub_height, cT, TSR, axial_induction)

    # turbine parameters
    D = rotor_diameter
    HH = hub_height
    
    ai = axial_induction
    avgW = W
    # print("avgW: ", avgW, "\n")
    # flow parameters
    Uinf = U_inf_initial
    # print("U_inf_initial: ", U_inf_initial)
    # print("U_inf: ", U_inf)

    # epsilon gain
    eps_gain = 0.2
    eps = eps_gain * D

    # # find velocity at top and bottom of rotor swept area for yaw
    # vel_top = ((HH+D/2)/HH)^0.12
    # vel_bottom = ((HH-D/2)/HH)^0.12

    

    # findzeros(tilt_eval, -45*pi/180, 45*pi/180)
    added_tilt, info = brent(tilt_eval, -45.0*pi/180, 45.0*pi/180, atol=0.00000000001, args=(D, Uinf, cT, ai, TSR, z_i, HH, deltay, eps, avgW, U_inf))
    # print("info: ", info, "\n")
    # # Define range of test values for tilt
    # tilt = -45.0:0.1:45
    # tilt = tilt.*pi/180
    # minTilt = 0.001
    # target_tilt_index = -10000.0
    # for i=1:length(tilt)
    #     Gamma_left = (pi/8)*D*Uinf*cT*sin(tilt[i])*cos(tilt[i])
    #     Gamma_right = (pi/8)*D*Uinf*cT*sin(tilt[i])*cos(tilt[i])*-1

    #     # Use turbine average velocity to find Gamma due to wake rotation
    #     turbine_average_velocity = U_inf
    #     Gamma_wake_rotation = 0.25*2.0*pi*D*(ai-ai^2)*turbine_average_velocity/TSR

    #     # Left vortex
    #     dz_left = z_i- HH
    #     dy_left = deltay-(D/2)
    #     r_l = dy_left^2 + dz_left^2
    #     core_shape = 1-exp(-r_l/(eps^2))     
    #     w_left = ((-1*Gamma_left*dy_left)/(2*pi*r_l))*core_shape

    #     # right vortex
    #     dz_right = z_i - HH
    #     dy_right = deltay+(D/2)
    #     r_r = dy_right^2 + dz_right^2
    #     core_shape = 1-exp(-r_r/(eps^2))     
    #     w_right = ((-1*Gamma_right*dy_right)/(2*pi*r_r))*core_shape

    #     # wake rotation vortex
    #     dz_center = z_i - HH
    #     r_center = deltay^2 + dz_center^2
    #     core_shape = 1 - exp(-r_center/(eps^2))
    #     w_center = ((-1*Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape

    #     ### Boundary condition - ground mirror vortex
    #     # Left vortex
    #     dz_left = z_i + HH
    #     dy_left = deltay-(D/2)
    #     r_l = dy_left^2 + dz_left^2
    #     core_shape = 1-exp(-r_l/(eps^2))     
    #     w_left_g = ((Gamma_left*dy_left)/(2*pi*r_l))*core_shape

    #     # right vortex
    #     dz_right = z_i + HH
    #     dy_right = deltay+(D/2)
    #     r_r = dy_right^2 + dz_right^2
    #     core_shape = 1-exp(-r_r/(eps^2))     
    #     w_right_g = ((Gamma_right*dy_right)/(2*pi*r_r))*core_shape

    #     # wake rotation vortex
    #     dz_center = z_i + HH
    #     r_center = deltay^2 + dz_center^2
    #     core_shape = 1 - exp(-r_center/(eps^2))
    #     w_center_g = ((Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape

    #     W_eff = w_left+w_right+w_center + w_left_g+w_right_g+w_center_g
    #     # W_eff = w_left+w_right+w_center 
        
    #     tmp = abs(avgW - W_eff)
    #     if tmp < minTilt
    #         target_tilt_index = i
    #         # print("tilt[i]: ", tilt[i], "\n")
    #     end
    # end

    # if target_tilt_index == -10000.0
    #     # print("Error: no effective tilt found, setting it to 0")
    #     added_tilt = 0.0
    # else
    #     added_tilt = tilt[target_tilt_index]
    # end

    # print("added_tilt: ", added_tilt, "\n")
    # # # find Gamma at the top and bottom of the rotor swept area
    # # # Gamma_top = gamma(D, vel_top, Uinf, cT)
    # # Gamma_top = (pi/8)*D*Uinf*cT*vel_top*sin(tilt)*cos(tilt)                # why do we use Uinf here, should we use U_inf?
    # # # Gamma_bot = -1*gamma(D, vel_bottom, Uinf, cT)
    # # Gamma_bot = (pi/8)*D*Uinf*cT*vel_bottom*-1.0*sin(tilt)*cos(tilt)

    # Gamma_left = (pi/8)*D*Uinf*cT*sin(tilt)*cos(tilt)
    # Gamma_right = (pi/8)*D*Uinf*cT*sin(tilt)*cos(tilt)*-1

    # # Use turbine average velocity to find Gamma due to wake rotation
    # turbine_average_velocity = U_inf
    # Gamma_wake_rotation = 0.25*2.0*pi*D*(ai-ai^2)*turbine_average_velocity/TSR           # why in FLORIS do they multiply this by 0.25*2?
    # # Gamma_wake_rotation = pi*D*(ai-ai^2)*turbine_average_velocity/TSR

    # ### Calculate the spanwise and vertical velocities induced by tilt ###
    # # top vortex
    # dz_top = z_i - (HH+D/2)             
    # r_top = deltay^2 + dz_top^2             
    # core_shape = 1-exp(-r_top/(eps^2))          
    # w_top = ((-1*Gamma_top*deltay)/(2*pi*r_top))*core_shape
    # # w_top = ((Gamma_top*deltay)/(2*pi*r_top))*core_shape


    # # bottom vortex
    # dz_bot = z_i - (HH - D/2)
    # r_bot = deltay^2 + dz_bot^2
    # core_shape = 1-exp(-r_bot/(eps^2))
    # w_bot = ((-1*Gamma_bot*deltay)/(2*pi*r_bot))*core_shape
    # # w_bot = ((Gamma_bot*deltay)/(2*pi*r_bot))*core_shape


    # # wake rotation vortex
    # dz_center = z_i - HH
    # r_center = deltay^2 + dz_center^2
    # core_shape = 1 - exp(-r_center/(eps^2))
    # # w_center = ((-1*Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape
    # w_center = ((Gamma_wake_rotation*deltay)/(2*pi*r_center))*core_shape


    # # print("deltay: ", deltay, "\n")
    # # print("Gamma_wake_rotation: ", Gamma_wake_rotation, "\n")
    # val = 2*(avgW-w_center)/(w_top+w_bot)           # why is this multiplied by 2?
    # # val = (avgW-w_center)/(w_top+w_bot)
    # # print("val: ", val, "\n")
    # # cap the added_tilt to be between -45 and 45, otherwise asin() will through an error
    # if val > 1.0
    #     val = 1.0
    # elseif val < -1.0
    #     val = -1.0
    # end

    # # print("val: ", val, "\n")
    # added_tilt = 0.5*asin(val)              # why is this multiplied by 0.5?

    # added_tilt = asin(val) 
    return added_tilt
end

"""
    point_velocity_tilt(loc, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities, W_sorted, V_sorted,
    wind_resource, nrotorpoints, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the wind speed at a given point for a given state (with tilt)

# Arguments
- `loc::Array{TF,3}`: Location of interest
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: turbine thrust coefficients for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: array containing indices of wind turbines 
    from most upwind to most downwind turbine in the given state
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `W_sorted::Array{TF,nTurbines}`: vertical spanwise velocities for each turbine
- `V_sorted::Array{TF,nTurbines}`: horizontal spanwise velocities for each turbine
- `shearexponent::Float`: shearexponent for determing wind speed at different heights
- `wind_resource::DiscretizedWindResource`: contains wind resource discreption (directions,
    speeds, frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `downwind_turbine_id::Int`: index of wind turbine of interest (if any). If not a point for
    calculating effective wind speed of a turbine, then provide 0 (default)
- `nrotorpoints::Int`: number of rotor points used to evaluate wind speed in rotor swept area
"""
function point_velocity_tilt(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities, W_sorted, V_sorted, shearexponent,
                    wind_resource, nrotorpoints, model_set::AbstractModelSet;
                    wind_farm_state_id=1, downwind_turbine_id=0)

    wakedeficitmodel = model_set.wake_deficit_model
    wakedeflectionmodel = model_set.wake_deflection_model
    wakecombinationmodel = model_set.wake_combination_model

    # The streamwise vertical and horizonal velocities must be found for each turbine
    # There is the effect of a turbine on itself
    # There is the effect of upstream turbines on a downstream turbine
    
    # extract flow information
    wind_speed = wind_resource.wind_speeds[wind_farm_state_id]
    reference_height = wind_resource.measurement_heights[wind_farm_state_id]

    # set ground height 
    ground_height = wind_resource.wind_shear_model.ground_height    # TODO: allow topology to be given

    # find order for wind shear and deficit calculations
    shear_order = wind_resource.wind_shear_model.shear_order

    # adjust wind speed for wind shear
    if shear_order == "nothing"
        wind_speed_internal = wind_speed
    elseif shear_order == "first"
        wind_speed_internal = adjust_for_wind_shear(locz, wind_speed, reference_height, ground_height, wind_resource.wind_shear_model)
    else
        wind_speed_internal = wind_speed
    end

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize the vertical and horizontal spanwise velocities and deficit summation term to zero
    deficit_sum = 0.0
    V = 0.0
    W = 0.0

    # print("NEW TURBINE", downwind_turbine_id, "\n")

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = Int(sorted_turbine_index[u])

        # don't allow turbine to impact itself
        # currently the code may be calculating the added vorticity from upstream
        # turbines, but it should also calculate the W and V at the turbine itself
        # if upwind_turb_id == downwind_turbine_id; continue; end

        # downstream distance between upstream turbine and point
        x = locx - turbine_x[upwind_turb_id]

        # check turbine relative locations
        # include if statement here for when we are looking at most upstream turbine
        # find W for this turbine to be updated ouside this point_velocity_tilt
        # Since Ct, ai, and ti aren't updated until dV is found, the v_wake and w_wake
        # of a turbine on itself will need to be found now that we are looking at a new
        # turbine
        # print("V_sorted before: ", V_sorted, "\n")
        if upwind_turb_id == downwind_turbine_id-1
            dx = 0.0
            dy = 0.0
            v_wake, w_wake = calculate_transverse_velocity(wtvelocities[upwind_turb_id], wind_speed, dx, dy, locz, rotor_diameter[upwind_turb_id], 
            hub_height[upwind_turb_id], turbine_tilt[upwind_turb_id], turbine_ct[upwind_turb_id], TSR, turbine_ai[upwind_turb_id], shearexponent)
            
            ### SOMETHING IS WRONG HERE ###
            ### SHOULDN'T BE JUST ADDING, NEED TO TAKE AVERAGE OF V_WAKE OVER ROTOR POINTS ###
            V_sorted[upwind_turb_id] = ((V_sorted[upwind_turb_id]*nrotorpoints) + v_wake)/(nrotorpoints + 1)
            W_sorted[upwind_turb_id] = ((W_sorted[upwind_turb_id]*nrotorpoints) + w_wake)/(nrotorpoints + 1)
            print("V_sorted: ", V_sorted[upwind_turb_id])
            print("W_sorted: ", W_sorted[upwind_turb_id])
        end
        # print("V_sorted after: ", V_sorted)

        if x > 1E-6
            # skip this loop if it would include a turbine's impact on itself)
            if upwind_turb_id==downwind_turbine_id; continue; end

            # find variables for wake_added_tilt function
            # spanwise distance between upstream turbine and point
            deltay = locy-turbine_y[upwind_turb_id]

            # find the added tilt angle due to the vortices
            # TODO: update wtvelocities and turbine tilt after all the comparisons
            # print("W_sorted: ", W_sorted, "\n")
            # Something is wrong with W_sorted not being updated correctly
            # This should be W_sorted from the upstream turbine to the downwind_turbine_id
            # This means that outside of this if statement W from the upstream turbine to all
            # downstream turbines should be calculated and updated in W_sorted.


            # Acutally since FLOWFarm steps through a sorted list of three turbines
            # and then compares only effects of upstream turbine on the turbine of interest
            # we can compute w and v in a seperate for loop before this for loop and 
            # use that to find the added_tilt. This computed w and v would go in the _gauss_yaw_potential_core
            # of W_sorted.

            # Remember that this is the added_tilt in the upwind_turb_id
            # That means we only need W of the upstream turbines, not the downstream turbine (turbine being looked at)
            # This means that the W for the most upstream turbine should be non-zero

            added_tilt = wake_added_tilt(wtvelocities[upwind_turb_id], W_sorted[upwind_turb_id], wind_speed, deltay, 
            locz, rotor_diameter[upwind_turb_id], hub_height[upwind_turb_id], turbine_ct[upwind_turb_id], TSR, turbine_ai[upwind_turb_id])
            
            

            # added_tilt only applies to calculating the wake deflection, no velocity deficit
            # the added ti and actual tilt angle are what are used to find the added wake recovery
            # in the velocity deficit equation. Therefore, set the added_tilt as a seperate variable.
            effective_tilt = deepcopy(turbine_tilt)
            effective_tilt[upwind_turb_id] = turbine_tilt[upwind_turb_id] + added_tilt
            # turbine_tilt[upwind_turb_id] += added_tilt
            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = 0.0

            vertical_deflection = wake_deflection_model(locx, locy, locz, turbine_x, effective_tilt, turbine_ct,
            upwind_turb_id, rotor_diameter, turbine_local_ti, wakedeflectionmodel)

            # find the vertical and horizontal spanwise velocities
            # Why is the adjusted turbulence not calculated before vertical_deflection is calculated???
                    # Potential reason is that the added_tilt accounts for the added deflection
            
            # Does TSR need to be updated based on new velocities?
            # This v_wake and w_wake are the vertical and horizontal spanwise velocities
            # induced on the turbine of interest by the upstream turbines
            dx = locx - turbine_x[upwind_turb_id]
            dy = deltay

            v_wake, w_wake = calculate_transverse_velocity(wtvelocities[upwind_turb_id], wind_speed, dx, dy, locz, rotor_diameter[upwind_turb_id], 
            hub_height[upwind_turb_id], turbine_tilt[upwind_turb_id], turbine_ct[upwind_turb_id], TSR, turbine_ai[upwind_turb_id], shearexponent)

            # print("v_wake: ", v_wake, "\n")
            # print("w_wake: ", w_wake, "\n")
            V += v_wake
            W += w_wake
            # print("W: ", W, "\n")
            # find the TI mixing induced by the tilt
            # Instead of v_wake, w_wake, this should be W_sorted[upwind_turb_id] and V_sorted[upwind_turb_id]
            I_mixing = tilt_added_turbulence_intensity(wtvelocities[upwind_turb_id], W_sorted[upwind_turb_id], turbine_local_ti[upwind_turb_id], 
            V_sorted[upwind_turb_id], v_wake, w_wake)
            # print("I_mixing: ", I_mixing, "\n")

            # What to set gch_gain to?
            gch_gain = 2
            turbine_local_ti[upwind_turb_id] += gch_gain*(I_mixing)
            # velocity difference in the wake
            deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                            upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                            turbine_local_ti, turbine_ct, turbine_tilt, wakedeficitmodel)

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, wakecombinationmodel)
            # println(deficit_sum, " ", downwind_turbine_id, " ", upwind_turb_id)
            
        end
    end

    # find velocity at point without shear
    point_velocity = wind_speed_internal - deficit_sum

    if shear_order == "nothing"
        point_velocity_out = point_velocity
    elseif shear_order == "first"
        point_velocity_out = point_velocity
    else
        point_velocity_out = adjust_for_wind_shear(locz, point_velocity, reference_height, ground_height, wind_resource.wind_shear_model)        
    end
    
    # this will return W and V at the turbine being looked at
    # so that it can be updated when it is being looked at as an upstream
    # turbine
    return point_velocity_out, W_sorted, V_sorted

end

"""
    point_velocity(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id=1)

Calculates the wind speed at a given point for a given state

# Arguments
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine sorted order upstream to downstream 
    for given state
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes 
    with state etc
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across 
    the rotor swept area when calculating the effective wind speed for the wind turbine. 
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades) 
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the 
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)
    
    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    # initialize correct array types
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_yaw[1]))

    # initialize arrays
    turbine_velocities = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    # loop over all turbines
    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0
        # turbine_velocities[downwind_turbine_id] = 0.0

        # loop over all rotor sample points to approximate the effective inflow velocity
        for p=1:n_rotor_sample_points

            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            # put rotor sample points in wind direction coordinate system, and account for yaw
            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    wind_resource, model_set,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_yaw[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        
        # update local turbulence intensity for downstream turbine
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                            turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

        # println("local ti turb 9: ", turbine_local_ti[downwind_turbine_id])
    end

    # df = DataFrame(ID=1:n_turbines, V=turbine_velocities, TI=turbine_local_ti, CT=turbine_ct)
    # CSV.write("internaldata.txt", df)

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end
end

"""
    turbine_velocities_one_direction_tilt(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource, shearexponnent,
    model_set::AbstractModelSet; wind_farm_state_id=1)

Calculates the wind speed at a given point for a given state

# Arguments
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
    radians
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine sorted order upstream to downstream 
    for given state
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes 
    with state etc
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across 
    the rotor swept area when calculating the effective wind speed for the wind turbine. 
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades) 
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the 
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `shearexponent::Float`: used to determine wind speed at different height using wind power law
"""
function turbine_velocities_one_direction_tilt(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource, shearexponent,
                    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)
    
    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    # initialize correct array types
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_tilt[1]))

    # initialize arrays
    turbine_velocities = zeros(arr_type, n_turbines)
    W_sorted = zeros(arr_type, n_turbines)
    V_sorted = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    # loop over all turbines
    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0

        # initialize vertical and horizontal spanwise velocity to zero
        W_wake = zeros(arr_type, n_turbines)
        V_wake = zeros(arr_type, n_turbines)

        # initialize adjusted tilt to zero
        # Tilt = zeros(arr_type, n_turbines)
        # turbine_velocities[downwind_turbine_id] = 0.0

        # Calculate initial Ct, ai, and ti for most upstream turbine because
        # these values are needed in finding W and V in point_velocity_tilt
        # if d == 1

        #     ##### IS THE TURBINE TILT HERE THE EFFECTIVE TILT OR NORMAL TILT?   ######
        #     ##### PRETTY SURE THE EFFECTIVE TILT IS ONLY USED IN DEFLECTION CALC #####
        #     velocities = wind_resource.wind_speeds[wind_farm_state_id]
        #     # update thrust coefficient for downstream turbine
        #     turbine_ct[downwind_turbine_id] = calculate_ct(velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        #     # update axial induction for downstream turbine
        #     turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_tilt[downwind_turbine_id])

        #     # get local turbulence intensity for this wind state
        #     ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
            
        #     # update local turbulence intensity for downstream turbine
        #     turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_tilt, turbine_local_ti, sorted_turbine_index,
        #                         velocities, turbine_ct, turbine_tilt, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)
        # end

        # loop over all rotor sample points to approximate the effective inflow velocity
        for p=1:n_rotor_sample_points

            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            # put rotor sample points in wind direction coordinate system, and account for yaw
            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_tilt[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] .+ local_rotor_sample_point_z*cos(turbine_tilt[downwind_turbine_id])

            # calculate the velocity at given point
            point_velocity_with_shear, W_adjusted, V_adjusted = point_velocity_tilt(locx, locy, locz, turbine_x, 
                                    turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    W_sorted, V_sorted, shearexponent, wind_resource, n_rotor_sample_points, model_set,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear
            # do the same with vertical and horizontal spanwise velocities and tilt
            # print("W_adjusted: ", W_adjusted)
            # print("W_wake: ", W_wake)
            W_wake .+= W_adjusted
            V_wake .+= V_adjusted
            # Tilt .+= tilt_adjusted


        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points
        W_wake /= n_rotor_sample_points
        V_wake /= n_rotor_sample_points
        # Tilt /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = wind_turbine_velocity

        # update tilt of upstream/all turbines
        # turbine_tilt = Tilt

        # update vertical and horizontal spanwise velocities for downwind turbine
        W_sorted = W_wake
        V_sorted = V_wake
        # print("W_sorted_out: ", W_sorted, "\n")
        # print("V_sorted_out: ", V_sorted, "\n")

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_tilt[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        
        # update local turbulence intensity for downstream turbine
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_tilt, turbine_local_ti, sorted_turbine_index,
                            turbine_velocities, turbine_ct, turbine_tilt, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

        # println("local ti turb 9: ", turbine_local_ti[downwind_turbine_id])
    end

    # df = DataFrame(ID=1:n_turbines, V=turbine_velocities, TI=turbine_local_ti, CT=turbine_ct)
    # CSV.write("internaldata.txt", df)

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti, W_sorted, V_sorted, turbine_tilt
    end
end

function turbine_velocities_one_direction(x, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)

    n_turbines = Int(length(x)/2)
    # println(typeof(x), n_turbines)
    turbine_x = x[1:n_turbines] 
    turbine_y = x[n_turbines+1:end]
    # println(turbine_x)
    # println("turbine_x type ", typeof(turbine_x))
    # println("type of x ", typeof(x))

    # get number of turbines and rotor sample point
    # n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                typeof(hub_height[1]),typeof(turbine_yaw[1]))
    turbine_velocities = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0
        # turbine_velocities[downwind_turbine_id] = 0.0

        for p=1:n_rotor_sample_points


            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                wind_resource, model_set,
                                wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_yaw[downwind_turbine_id])

        # update local turbulence intensity for downstream turbine
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_velocities, turbine_ct, turbine_yaw, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

    end

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end

end

# turbine_velocities_one_direction!(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1) = turbine_velocities_one_direction!([0.0], [0.0],
# model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)

"""
calculate_flow_field(xrange, yrange, zrange, model_set::AbstractModelSet, turbine_x, 
    turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, rotor_diameter, hub_height, 
    turbine_local_ti, sorted_turbine_index, wtvelocities, wind_resource; wind_farm_state_id=1)

Generates a flow field for a given state and cross section

# Arguments
- `xrange::Range`: range defining east-west locations to sample in global reference frame
- `yrange::Range`: range defining north-west locations to sample in global reference frame
- `zrange::Range`: range defining vertical locations to sample in global reference frame
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the global reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: thrust coefficient of each turbine for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource; wind_farm_state_id=1)

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi] = point_velocity(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)

            end
        end
    end

    

    # if zlen == 1
    #     return point_velocities[1,1:ylen,1:xlen]
    # elseif ylen == 1
    #     return point_velocities[1:zlen,1,1:xlen]
    # elseif xlen == 1
    #     return point_velocities[1:zlen,1:ylen,1]
    # else
    return point_velocities[1:zlen,1:ylen,1:xlen]
    # end

end

function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw,
    rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = turbine_velocities_one_direction(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, wind_farm_state_id=wind_farm_state_id)

end

"""
calculate_flow_field_tilt(xrange, yrange, zrange, model_set::AbstractModelSet, turbine_x, 
    turbine_y, turbine_z, turbine_tilt, turbine_ct, turbine_ai, rotor_diameter, hub_height, 
    turbine_local_ti, sorted_turbine_index, wtvelocities, wind_resource; wind_farm_state_id=1)

Generates a flow field for a given state and cross section

# Arguments
- `xrange::Range`: range defining east-west locations to sample in global reference frame
- `yrange::Range`: range defining north-west locations to sample in global reference frame
- `zrange::Range`: range defining vertical locations to sample in global reference frame
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the global reference frame
- `turbine_tilt::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: thrust coefficient of each turbine for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `shearexponent::Float`: exponent used in wind power law
"""

function calculate_flow_field_tilt(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, W_sorted, V_sorted, shearexponent; wind_farm_state_id=1)

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    # variable not used
    nrotorsamplepoints = 1
    
    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi], W, V = point_velocity_tilt(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities, W_sorted, V_sorted, shearexponent,
                    wind_resource, nrotorsamplepoints, model_set,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)
                ### W_sorted and V_sorted would change here because we are looking at spaced between turbines
                ### instead of just turbines.
                    # point_velocity_tilt(locx, locy, locz, turbine_x, 
                    #                 turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                    #                 rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                    #                 W_sorted, V_sorted, wind_resource, model_set,
                    #                 wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            end
        end
    end

    

    # if zlen == 1
    #     return point_velocities[1,1:ylen,1:xlen]
    # elseif ylen == 1
    #     return point_velocities[1:zlen,1,1:xlen]
    # elseif xlen == 1
    #     return point_velocities[1:zlen,1:ylen,1]
    # else
    return point_velocities[1:zlen,1:ylen,1:xlen]
    # end

end

function calculate_flow_field_tilt(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_tilt,
    rotor_diameter, hub_height, ct_models, TSR, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource, shearexponent; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti, W_sorted, V_sorted, turbine_tilt = turbine_velocities_one_direction_tilt(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource, shearexponent,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field_tilt(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, W_sorted, V_sorted, shearexponent, wind_farm_state_id=wind_farm_state_id)

end

