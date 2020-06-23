using Snopt
using DiffResults, ForwardDiff

function rosenbrock(x)

    f = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2

    return f 
    
end

function opt_wrapper(x)
    result = DiffResults.GradientResult(x)
    result = ForwardDiff.gradient!(result, rosenbrock, x)
    f = DiffResults.value(result)
    g = DiffResults.gradient(result)

    fail = false

    c = []
    dcdx = []

    return f, c, g, dcdx, fail
end



x0 = [4.0; 4.0]
lb = [-5.0; -5.0]
ub = [5.0; 5.0]
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 1
options["Major optimality tolerance"] = 1e-6

xopt, fopt, info = snopt(opt_wrapper, x0, lb, ub, options)