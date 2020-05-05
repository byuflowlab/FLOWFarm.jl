abstract type AbstractLocalTurbulenceIntensityModel end

struct LocalTIModelNoLocalTI{} <: AbstractLocalTurbulenceIntensityModel

end

struct LocalTIModelMaxTI{TF} <: AbstractLocalTurbulenceIntensityModel
    ambient_ti::TF
end

struct LocalTIModelGaussTI{TF} <: AbstractLocalTurbulenceIntensityModel
    ti::TF
end

function calculate_local_ti(ambient_ti, ti_model::LocalTIModelNoLocalTI)
    return ambient_ti
end

function calculate_local_ti(ambient_ti, windfarmstate, ti_model::LocalTIModelMaxTI; turbine_id=1)

    # calculate local turbulence intensity at turbI
    
    # initialize the TI_area_ratio to 0.0 for each turbine
    TI_area_ratio = 0.0

    # initialize local ti tmp for downstream turbine
    ti_dst_tmp = windfarmstate.turbine_local_ti[turbine_id]

    # extract the number of turbines
    nturbines = length(ti_dst_tmp)
    
    # loop over upstream turbines
    for u=1:nturbines
    
        # get index of upstream turbine
        turb = windfarmstate.sorted_turbine_index[u]
        
        # skip turbine's influence on itself
        if (turb .eq. turbI) cycle
        
        ! calculate downstream distance between wind turbines
        x = turbineXw(turbI) - turbineXw(turb)
        
        if (x > tol) then
            ! determine the far-wake onset location 
            call x0_func(rotorDiameter(turb), yaw(turb), Ct_local(turb), alpha, & 
                        & TIturbs(turb), beta, x0)
        
            ! calculate the distance from the onset of far-wake
            deltax0 = x - x0
        
            ! horizontal spread 
            call sigmay_func(x, x0, ky_local(turb), rotorDiameter(turb), yaw(turb), sigmay)

            ! vertical spread 
            call sigmaz_func(x, x0, kz_local(turb), rotorDiameter(turb), sigmaz)
        
            ! determine the initial wake angle at the onset of far wake
            call theta_c_0_func(yaw(turb), Ct_local(turb), theta_c_0)
    
            ! horizontal cross-wind wake displacement from hub
            call wake_offset_func(x, rotorDiameter(turb), theta_c_0, x0, yaw(turb), &
                                    & ky_local(turb), kz_local(turb), &
                                    Ct_local(turb), sigmay, sigmaz, wake_offset)
    
            ! cross wind distance from point location to upstream turbine wake center
            deltay = turbineYw(turbI) - (turbineYw(turb) + wake_offset)   
        
            ! save ti_area_ratio and ti_dst to new memory locations to avoid 
            ! aliasing during differentiation
            TI_area_ratio_tmp = TI_area_ratio
            TI_dst_tmp = TIturbs(turbI)
            TI_ust_tmp = TIturbs(turb)
    
            ! update local turbulence intensity
            call added_ti_func(TI, Ct_local(turb), x, ky_local(turb), rotorDiameter(turb), & 
                                & rotorDiameter(turbI), deltay, turbineZ(turb), &
                                & turbineZ(turbI), sm_smoothing, TI_ust_tmp, &
                                & TI_calculation_method, TI_area_ratio_tmp, &
                                & TI_dst_tmp, TI_area_ratio, TIturbs(turbI))
        end if
    
    end do
    
    ! calculate wake spreading parameter at turbI based on local turbulence intensity
    if (calc_k_star .eqv. .true.) then

        call k_star_func(TIturbs(turbI), k_star)
        ky_local(turbI) = k_star
        kz_local(turbI) = k_star

    end if
        
end