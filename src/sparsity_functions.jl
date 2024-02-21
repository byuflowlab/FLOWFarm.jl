"""file with functions used in wind farm optimization employing sparse methods
created January 26, 2024
author: Benjamin Varela
"""

abstract type AbstractSparseMethod end

#TODO: documentation and examples

# ∇AEP Optimization (single pattern) #######################################################


# ∇AEP Optimization (multi pattern) #######################################################


# Sparse spacing constraint methods #######################################################


# Sparse boundary constraint methods #######################################################
struct sparse_boundary_struct{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractSparseMethod
    turbine_x::T1
    turbine_y::T2
    jacobian::T3
    ad::T4
    cache::T5
    boundary_vec::T6
    boundary_function::T7
    update_function::T8
    boundary_scaling_factor::T9
end
