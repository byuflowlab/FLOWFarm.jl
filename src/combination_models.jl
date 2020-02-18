abstract type AbstractCombinationModel end

struct TwoNorm <: AbstractCombinationModel end

function combination_model(loss_vec, ::TwoNorm)

    loss = sqrt(sum(loss_vec.^2))
    return loss

end

# combination_model(loss_vec, TwoNorm())
