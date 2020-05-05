abstract type AbstractThrustCoefficientModel end

"""
    ThrustModelConstantCt(ct::Float)

Stores a constant ct value for wake calculations

# Arguments
- `ct::Float`: a constant ct value for computation
"""
struct ThrustModelConstantCt{TF} <: AbstractThrustCoefficientModel
    ct::TF
end

"""
    ThrustModelCtPoints(vel_points, ct_points)

Stores the thrust coefficient curve in terms of corresponding velocity and thrust
coefficient points. ct and velocity points should be in the same order and ordered from
lowest wind speed to highest wind speed.

# Arguments
- `inflow_velocity::Float`: inflow velocity of the wind turbine
- `thrust_model::ThrustModelCtPoints`: Struct containing ct and velocity points for ct curve
"""
struct ThrustModelCtPoints{ATF} <: AbstractThrustCoefficientModel
    vel_points::ATF
    ct_points::ATF
end

"""
    calculate_ct(model::ThrustModelConstantCt)

Calculate the thrust coefficient for a wind turbine based on a pre-determined constant ct

# Arguments
- `inflow_velocity::Float`: inflow velocity of the wind turbine (unused for const. ct)
- `thrust_model::ThrustModelConstantCt`: struct containing a constant ct value for computation
"""
function calculate_ct(inflow_velocity, thrust_model::ThrustModelConstantCt)
    return thrust_model.ct
end

"""
    calculate_ct(inflow_velocity, thrust_model::ThrustModelCtPoints)

    Calculate the thrust coefficient for a wind turbine based on a pre-determined ct curve 
        with linear interpolation.

# Arguments
- `inflow_velocity::Float`: inflow velocity of the wind turbine
- `thrust_model::ThrustModelCtPoints`: Struct containing ct and velocity points for ct curve
"""
function calculate_ct(inflow_velocity, thrust_model::ThrustModelCtPoints)

    minv = thrust_model.vel_points[1]
    maxv = thrust_model.vel_points[end]
    minct = thrust_model.ct_points[1]
    maxct = thrust_model.ct_points[end]

    if inflow_velocity < minv
        ct = minct
    elseif inflow_velocity < maxv
        ct = linear(thrust_model.vel_points, thrust_model.ct_points, inflow_velocity)
    else
        ct = maxct
    end
    
    return ct
end