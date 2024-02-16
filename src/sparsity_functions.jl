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
struct sparse_boundary_struct{T1,T2,T3,T4,T5} <: AbstractSparseMethod
    boundary_jacobian::T1
    ad::T2
    cache::T3
    boundary_vec::T4
    boundary_function::T5
end
# sparse_jacobian!(jacobian, ad, cache, f, vec, x)

function build_sparse_boundary_struct(x,boundary,farm)
    jacobian = dropzeros(sparse(boundary.boundary_jacobian))
    ad = AutoSparseForwardDiff()
    sd = JacPrototypeSparsityDetection(; jac_prototype=jacobian)
    vec = boundary.boundary_vec
    f(a,x) = calculate_boundary!(a,x,farm,boundary)
    cache = sparse_jacobian_cache(ad, sd, f, vec, x)
    return sparse_boundary_struct(jacobian,ad,cache,vec,f)
end
