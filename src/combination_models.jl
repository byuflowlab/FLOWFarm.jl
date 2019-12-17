function combination_model(loss_vec)

    loss = sqrt(sum(loss_vec.^2))
    return loss

end
