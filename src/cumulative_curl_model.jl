using Statistics
function turbine_velocities_one_direction_CC!(turbine_x::Vector{T0}, turbine_y::Vector{T1}, turbine_z::Vector{T2}, rotor_diameter::Vector{T3}, hub_height::Vector{T4}, turbine_yaw::Vector{T5},
    sorted_turbine_index::Vector{Int}, ct_model::Vector{<:AbstractThrustCoefficientModel}, rotor_sample_points_y::Vector{T6}, rotor_sample_points_z::Vector{T6}, wind_resource,
    model_set::AbstractModelSet, turbine_velocities::Vector{T7},
    turbine_ct::Vector{T7}, turbine_ai::Vector{T7}, turbine_local_ti::Vector{T7}; wind_farm_state_id::Int=1, velocity_only::Bool=true) where {T0, T1, T2, T3, T4, T5, T6, T7}

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                                typeof(hub_height[1]),typeof(turbine_yaw[1]))

    U_inf = wind_resource.wind_speeds[wind_farm_state_id]

    # get number of turbines, rotor sample points, and initialize contribution vector (see CC model)
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)
    C = zeros(arr_type,n_turbines,n_turbines)
    point_velocities = zeros(arr_type,n_turbines,n_rotor_sample_points)
    point_velocities .= U_inf
    deficits = zeros(arr_type,n_turbines,n_rotor_sample_points)
    deflections = zeros(arr_type,n_turbines,n_turbines)
    sigma = zeros(arr_type,n_turbines,n_turbines)

    ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]

    a_f = model_set.wake_deficit_model.a_f
    b_f = model_set.wake_deficit_model.b_f
    c_f = model_set.wake_deficit_model.c_f

    # loop over all turbines (n)
    for n=1:n_turbines
        current_turbine_id = Int(sorted_turbine_index[n])
        x_n = turbine_x[current_turbine_id]
        y_n = turbine_y[current_turbine_id]
        z_n = turbine_z[current_turbine_id] + hub_height[current_turbine_id]

        # update current turbine velocity from sample points
        tot = 0
        for t = 1:n_rotor_sample_points
            tot += point_velocities[current_turbine_id,t]
        end
        turbine_velocities[current_turbine_id] = tot / n_rotor_sample_points

        # update coefficient of thrust for current turbine
        turbine_ct[current_turbine_id] = calculate_ct(turbine_velocities[current_turbine_id], ct_model[current_turbine_id])

        # update local TI for current turbine
        turbine_local_ti[current_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                                        turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=current_turbine_id, tol=1E-6)

        for d = n+1:n_turbines
            downwind_turbine_id = Int(sorted_turbine_index[d])
            for p = 1:n_rotor_sample_points
                # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
                local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
                local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

                # put rotor sample points in wind direction coordinate system, and account for yaw
                x = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
                y = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
                z = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

                x_tilde_n = (x - x_n) / rotor_diameter[current_turbine_id]

                m = a_f*exp(b_f*x_tilde_n)+c_f
                a1 = 2^(2/m - 1)
                a2 = 2^(4/m - 2)

                if p == 1
                    dy = wake_deflection_model(x, y, z, turbine_x, turbine_yaw, turbine_ct,
                        current_turbine_id, rotor_diameter, turbine_local_ti, model_set.wake_deflection_model)
                    deflections[current_turbine_id,downwind_turbine_id] = dy

                    sigma_n = wake_expansion(turbine_ct[current_turbine_id],turbine_local_ti[current_turbine_id],x_tilde_n,model_set.wake_deficit_model)
                    sigma[current_turbine_id,downwind_turbine_id] = sigma_n
                else
                    dy = deflections[current_turbine_id,downwind_turbine_id]
                    sigma_n = sigma[current_turbine_id,downwind_turbine_id]
                end

                sum_C = 0
                for i = 1:n-1
                    other_turbine_id = Int(sorted_turbine_index[i])

                    y_i = turbine_y[other_turbine_id]
                    z_i = turbine_z[other_turbine_id] + hub_height[other_turbine_id]

                    sigma_i = sigma[other_turbine_id,current_turbine_id]
                    
                    dy_i = deflections[other_turbine_id,current_turbine_id]

                    lambda_n_i = sigma_n^2/(sigma_n^2+sigma_i^2) * exp(-((y_n-y_i-dy_i)^2 + (z_n-z_i)^2)/(2*(sigma_n^2+sigma_i^2)))

                    sum_C += lambda_n_i * C[other_turbine_id,current_turbine_id]
                end
                calc = a2 - (m*turbine_ct[current_turbine_id]*cos(turbine_yaw[current_turbine_id]))/(16.0*gamma(2/m)*sign(sigma_n)*(abs(sigma_n)^(4/m))*(1-sum_C/U_inf)^2)
                if calc < 0
                    calc = 0
                end
                C_point = (1-sum_C/U_inf) * (a1-sqrt(calc))

                r_tilde = sqrt((y-y_n-dy)^2 + (z-z_n)^2)/rotor_diameter[current_turbine_id]

                velDef = C_point*exp(-1 * (r_tilde^m)/(2*sigma_n^2))
                deficits[downwind_turbine_id,p] += velDef * turbine_velocities[current_turbine_id]
                if d == n+1
                    point_velocities[downwind_turbine_id,p] = U_inf - deficits[downwind_turbine_id,p]
                end
                if p == 1
                    C[current_turbine_id,downwind_turbine_id] = C_point
                end
            end
        end
    end
end

function wake_expansion(Ct,TI,x_tilde,model)
    beta = 0.5*(1.0+sqrt(1.0-Ct))/(sqrt(1.0-Ct))
    epsilon = (model.c_s1*Ct+model.c_s2)*sqrt(beta)
    k = model.a_s*TI+model.b_s
    sigma = k*x_tilde+epsilon
    return sigma
end