# Wind Farm Struct

This tutorial covers the usage of the wind farm struct. The wind farm struct simplifies optimization set up and allows for the use of sparse methods in gradient calculation. The AEP value and the gradient are both stored in the farm struct

## Farm Struct
All of the setup from the previous tutorial sections 1 and 2 is required. Then the wind farm struct in defined.
```
function my_update_function(farm,x)
    farm.turbine_x .= x[1:length(x)÷2]
    farm.turbine_y .= x[length(x)÷2+1:end]
end

farm = ff.build_wind_farm_struct(x0, x, y, turbine_z, hub_height, turbine_yaw, rotor_diameter, ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set, my_update_function; AEP_scale=1.0, input_type="ForwardDiff", opt_x=true, opt_y=true)
```

build_wind_farm_struct has many inputs and options. While most are well defined in the previous tutorial, some are new. 
- `x0` is a vector of the design variables and should already be scaled. 
- `my_update_function` is a function that is of the form `f(wind_farm_struct, design_variables)` and updates the wind farm struct in place from the design variables.
- `AEP_scale` applies a scale factor to the AEP calculation and by default uses `1.0/ideal_AEP`.
- `input_type` is FowardDiff when using the farm struct for optimization. Otherwise, deafults to `nothing` and will use the type of `x0`
- `opt_x`, `opt_y`, `opt_hub`, `opt_yaw`, and `opt_diam` determine which variables are being optimized so that only those variables are stored as `input_type`. This allows for faster gradient computation.

## Constraint Structs
Constraint structs are also provided for both spacing and boundary constraints.
```
spacing = 2.0 * rotor_diameter[1] # 2 rotor diameters spacing between turbines
scaling = 1.0/rotor_diameter[1]^2 # scaling factor for spacing constraint
spacing_struct = ff.build_spacing_struct(x0, num_turbines, spacing, scaling, my_update_function)

boundary_function(a,x,y) = ff.circle_boundary!(a, boundary_center, boundary_radius, x, y) # boundary function, must be in-place and take the form f(a,x,y) where a is the constraint values
n_constraints = num_turbines # number of boundary constraints
scaling = 1.0/boundary_radius^2 # scaling factor for boundary constraint
boundary_struct = ff.build_boundary_struct(x0, num_turbines, n_constraints, scaling, boundary_function, my_update_function)
```

The resulting objective function then looks like:
```
function opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)

    # calculate spacing constraint value and jacobian
    ff.calculate_spacing_jacobian!(spacing_struct,x)

    # calculate boundary constraint and jacobian
    ff.calculate_boundary_jacobian!(boundary_struct,x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_struct.spacing_vec; boundary_struct.boundary_vec]
    g[:] .= ((c[:]))
    dg[:] = (([spacing_struct.jacobian; boundary_struct.jacobian]))

    # calculate AEP and gradient
    ff.calculate_aep_gradient!(farm,x)

    AEP = -farm.AEP[1]
    df[:] .= ((-farm.AEP_gradient))

    return AEP
end

# Force signature to allow for pass into SNOW.minimize
opt!(g,df,dg,x) = opt!(g,df,dg,x,farm,spacing_struct,boundary_struct)
```

# Sparsity Options
Using sparsity can greatly accelerate wind farm optimization. For detailed theory please see
```
Varela, B., and Ning, A., “Sparsity Applications for Gradient-Based Optimization of Wind Farms,” Oct. 2024, (in review).
```
The use of sparsity is built into FLOWFarm for user convenience.

## Farm Sparsity
There are two types of sparsity used in calculating the gradient of the AEP, stable and unstable sparsity. Stable sparsity refers to the situation where the sparsity pattern does not change during the optimization. Unstable sparsity refers to situations where the sparsity pattern may change during the optimzation. The sparsity pattern referred to depends on how each turbine affects the power produced by every other turbine due to its wake interference. 

For example, any optimization where the x or y positions of the turbines are design variables will result in changing sparsity patterns as the turbines move and interact with different turbines. In contrast a yaw optimization is unlikely to cause a change in the sparsity pattern during the optimization.

Use of sparsity in AEP optimization requires a sparse_struct that must be passed into the objective function in the same manner as constraint structs. The AEP and the AEP gradient are then stored in the farm struct just as with non-sparse methods

```
ff.calculate_aep_gradient!(farm,x0,sparse_struct)
println("AEP = ",farm.AEP[1])
println("AEP Gradient = ",farm.AEP_gradient)
```

### Stable Sparsity
Stable sparse farm structs are constructed in a similar method to regular farm structs
```
farm, sparse_struct = ff.build_stable_sparse_struct(x0, x, y, turbine_z, hub_height, turbine_yaw, rotor_diameter, ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set, my_update_function_yaw; AEP_scale=1.0, opt_yaw=true, tolerance=1E-16)
```

The only new parameter is `tolerance` which defaults to 1E-16. This tolerance is used when computing the sparsity pattern. When creating the `sparse_struct` the sparsity pattern is computed by calculating the partial derivatives of the power of each turbine with respect to each design variable. Then, all partial derivatives with a magnitude below the `tolerance` are ignored in the sparsity pattern. The default value is safe for almost all cases and the effect of the `tolerance` is explored more in the above paper. 

Important note: `build_stable_sparse_struct` attempts to calculate an accurate sparsity pattern by perturbing all of the design variables randomly to avoid poor starting locations where derivatives zero out, such as in yaw optimization when the turbine faces directly into the wind. While this does help get an accurate sparsity pattern it is advised to call `build_stable_sparse_struct` with `x0` set to a reasonable start point to avoid derivatives that are not representative of the influence of the design variable.

### Unstable Sparsity
Unstable sparse methods can always be used as the sparsity pattern is recomputed at each iteration using the wake deficits from the AEP calculation, but unstable sparse methods are not as efficient as stable sparse methods, again see the theory paper.

Unstable sparse structs are computed in the same way as stable sparse struct with a different function call.
```
farm, sparse_struct = ff.build_unstable_sparse_struct(x0, x, y, turbine_z, hub_height, turbine_yaw, rotor_diameter, ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set, my_update_function_yaw; AEP_scale=1.0, opt_yaw=true, tolerance=1E-16)
```

## Constraint Sparsity

### Boundary Constraints
The boundary constraint structs use sparsity by default as it has not effect on the solutions and little effect on the computation speed. This can be turned off with `using_sparsity=false`.

### Spacing Constraints
The use of sparsity in spacing constraints is more complicated than with the farm AEP and requires at least 2 optimizations. First an optimization is performed with no spacing constraints, then a second optimization is performed with select spacing constraints.

The function call for building a sprase spacing struct is
```
sparse_spacing_struct = ff.build_sparse_spacing_struct(x0, x, y, spacing, scaling, spacing_update_function!; first_opt=true, relevant_spacing_factor=2)

ff.calculate_spacing_jacobian!(sparse_spacing_struct,x)
```
`first_opt` defaults to true and causes no spacing constraints to be calculated, resulting in empty spacing jacobian and spacing constraints objects. This is done to allow the objective function to remain unchanged between the first and second optimizations. `first_opt` should be set to false when creating the spacing struct for the second optimization
`relevant_spacing_factor` defaults to 2 and determines which spacing constraints are relevant during the optimization by only including constraints where turbine pairs are within `spacing` * `relevant_spacing_factor` from eachother. These constraints are determined based on the x and y passed into the build function.

By using no spacing constraints in the first optimization the design space is freed from many constraints, causing the solution to be improved and reached faster. The second optimization is then short as an optimum has been reached and only a small number of spacing constraints must be satisfied. This works well, but depending on your optimizer the second optimization can be difficult if started from a position that violates constraints. To ease this sparse spacing structs can store the last position that satisfied all the spacing constraints during the first optimization. 
```
ff.update_safe_design_variables!(sparse_spacing_struct,x)
```
`spacing_struct.safe_design_variables` can then be used as the start point for the second optimization.

The computing of the spacing constraints during the first optimization to determine if any are violated is minimal as the speed up of the optimization is caused by the optimizer having less constraints to satisfy.

## Sparse Objective Function in SNOW
When giving constraints to SNOW.jl it is possible to give the optimizer a sparse set of constaints whether or not the sparse structs are used. This is shown in the following using `sp`. 

```
function opt!(g,df,dg,x,farm,spacing_struct,boundary_struct,sp,sparse_struct)

    ff.update_safe_design_variables!(spacing_struct,x)

    # calculate spacing constraint value and jacobian
    ff.calculate_spacing_jacobian!(spacing_struct,x)

    # calculate boundary constraint and jacobian
    ff.calculate_boundary_jacobian!(boundary_struct,x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_struct.spacing_vec; boundary_struct.boundary_vec]
    g[:] .= ((c[:]))

    dcdx = (sparse([spacing_struct.jacobian; boundary_struct.jacobian]))
    for i = 1:length(sp.rows)
        dg[i] = dcdx[sp.rows[i],sp.cols[i]]
    end

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    ff.calculate_aep_gradient!(farm,x,sparse_struct)

    AEP = -farm.AEP[1]

    df[:] .= ((-farm.AEP_gradient))

    return AEP
end
```

`sp` is calculated as follows assuming that the constraint jacobians have already been calculated.
```
sp = SNOW.SparsePattern(sparse(dropzeros(vcat(spacing_struct.jacobian,boundary_struct.jacobian))))
```