abstract type AbstractWakeCombinationModel end

#=
To add to wake combination models
Zong, H., & Fernando Porté-Agel. (2020). A momentum-conserving wake superposition method for wind farm power prediction. Journal of Fluid Mechanics, 889 doi:https://doi.org/10.1017/jfm.2020.77
=#

struct LinearFreestreamSuperposition <: AbstractWakeCombinationModel
    # Lissaman 1979
    # new_deficit_sum = old_deficit_sum + wind_speed*deltav

end

struct SumOfSquaresFreestreamSuperposition <: AbstractWakeCombinationModel
    # Katic et al. 1986
    # new_deficit_sum = sqrt(old_deficit_sum**2 + (wind_speed*deltav)**2)

end

struct SumOfSquaresLocalVelocitySuperposition <: AbstractWakeCombinationModel
    # Voutsinas 1990
    # new_deficit_sum = sqrt(old_deficit_sum**2 + (turb_inflow*deltav)**2)

end

struct LinearLocalVelocitySuperposition <: AbstractWakeCombinationModel
    # Niayifar and Porte Agel 2015, 2016
    # new_deficit_sum = old_deficit_sum + turb_inflow*deltav

end

function check_negative_deficits!(new_deficit_sum, wind_speed)
    if new_deficit_sum > wind_speed
        new_deficit_sum = wind_speed - eps(typeof(new_deficit_sum))
    end
end

function wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum,  model::LinearFreestreamSuperposition)
    # Lissaman 1979

    new_deficit_sum = old_deficit_sum + wind_speed*deltav

    check_negative_deficits!(new_deficit_sum, wind_speed)

    return new_deficit_sum

end

function wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model::SumOfSquaresFreestreamSuperposition)
    # Katic et al. 1986

    new_deficit_sum = nansafesqrt(old_deficit_sum^2 + (wind_speed*deltav)^2)

    check_negative_deficits!(new_deficit_sum, wind_speed)

    return new_deficit_sum

end

function wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model::SumOfSquaresLocalVelocitySuperposition)
    # Voutsinas 1990

    # square root of the sum of the squares
    new_deficit_sum = nansafesqrt(old_deficit_sum^2 + (turb_inflow*deltav)^2)
    # enforce non-negative wind speed
    check_negative_deficits!(new_deficit_sum, turb_inflow)
    return new_deficit_sum

end

function wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model::LinearLocalVelocitySuperposition)
    # Niayifar and Porte Agel 2015, 2016

    new_deficit_sum = old_deficit_sum + turb_inflow*deltav

    check_negative_deficits!(new_deficit_sum, turb_inflow)

    return new_deficit_sum

end
