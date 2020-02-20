abstract type AbstractWindResourceModel end

struct DiscretizedWindResource{AF, TF} <: AbstractWindResourceModel
    
    wind_directions::AF
    wind_speeds::AF
    wind_probabilities::AF
    measurement_heights::AF
    shear_exponent::TF

end