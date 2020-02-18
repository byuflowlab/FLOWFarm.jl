abstract type AbstractWindResourceModel end

struct DiscretizedWindResource{AF, TF} <: AbstractWindResourceModel
    
    winddirections::AF
    windspeeds::AF
    windprobabilities::AF
    measurementheight::AF
    shearexponent::TF

end