export normalize, unnormalize, lift_state_RBF, lift_data_RBF
export EDMD, TREDMD


"""
    EDMD(X::Matrix, en_level::Int64=6) → Ã, r

Extended Dynamic Mode Decomposition.

# Arguments
- `Xm::Matrix`: data matrix
- `Xp::Matrix`: data matrix (either x_{i+1} or dx/dt)
- `dt`: time step
- `en_level`: energy level for the reduced order (e.g. default 6 corresponds to 10^(-6)), if 0 we do not reduce the order

# Returns
- `Ã`: reduced order inferred state A operator from the lifted data
- `r`: reduced order

# References
J. L. Proctor, S. L. Brunton, and J. N. Kutz, “Dynamic Mode Decomposition with Control,” 
SIAM Journal on Applied Dynamical Systems, vol. 15, no. 1. Society for Industrial & 
Applied Mathematics (SIAM), pp. 142–161, Jan. 2016. doi: 10.1137/15m1013857.
"""
function EDMD(Xm::Matrix, Xp::Matrix, dt::Real, en_level::Int64=6)
    Nx = size(Xm, 1)
    Xm_svd = svd(Xm)
    Um = Xm_svd.U
    Σm = Xm_svd.S
    Vtm = Xm_svd.Vt

    if en_level != 0
        r_all = choose_ro(Σm)
        r = r_all[en_level] 
        Σm_inv = Diagonal(Σm) \ I(Nx)
        Ã = Um[:,1:r]' * Xp * Vtm[1:r,:]' * Σm_inv[1:r,1:r]

        # Compute the eigendecomposition of Ã        
        FOO = eigen(Ã)
        W = FOO.vectors
        Λ = FOO.values
        Ω = log.(Λ) / dt
        
        # Reconstruct the eigendecomposition of A 
        Φ = Xp * Vtm[1:r,:]' * Σm_inv[1:r,1:r] * W

        Ā = real.(Φ * Diagonal(Ω) / Φ)[1:r,1:r]

        return Φ, Ω, Ā, r
    else
        Σm_inv = Diagonal(Σm) \ I(Nx)
        Ā = Xp * Vtm' * Σm_inv * Um'
        FOO = eigen(Ā)
        W = FOO.vectors
        Λ = FOO.values
        Ω = log.(Λ) / dt
        Φ = Xp * Vtm' * Σm_inv * W
        Ā = real.(Φ * Diagonal(Ω) / Φ)

        return Φ, Ω, Ā, Nx
    end
end


"""
    TREDMD(Xm::Matrix, Xp::Matrix, δ::Real, en_level::Int64=6) → Ã, r

Tikhonov Regularized extended Dynamic Mode Decomposition with control.
The regularization is only applied to the states.

# Arguments
- `Xm`: data matrix
- `Xp`: data matrix (either x_{i+1} or dx/dt)
- `dt`: time step
- `δ`: the tikhonov regularization constant
- `en_level`: energy level for the reduced order (e.g. default 6 corresponds to 10^(-6)) when 0 no reduction

# Returns
- `Ã`: reduced order inferred state A operator from the lifted data
- `r`: reduced order
"""
function TREDMD(Xm::Matrix, Xp::Matrix, dt::Real, δ::Real, en_level::Int64=6)
    Nx = size(Xm, 1)
    Xm_svd = svd(Xm)
    Um = Xm_svd.U
    Σm = Xm_svd.S
    Vtm = Xm_svd.Vt

    if en_level != 0
        r_all = choose_ro(Σm)
        r = r_all[en_level] 
        Σm = Diagonal(Σm)
        Σm_inv = Diagonal(Σm) \ I(Nx)
        Reg = Σm[1:r,1:r].^2 + δ*I(r)
        Reg_inv = Reg \ I(r)
        Ã = Um[:,1:r]' * Xp * Vtm[1:r,:]' * Reg_inv

        # Compute the eigendecomposition of Ã        
        FOO = eigen(Ã)
        W = FOO.vectors
        Λ = FOO.values
        Ω = log.(Λ) / dt
        
        # Reconstruct the eigendecomposition of A 
        Φ = Xp * Vtm[1:r,:]' * Σm_inv[1:r,1:r] * W

        return Φ, Ω, Ã, r
    else
        Σm = Diagonal(Σm)
        Σm_inv = Diagonal(Σm) \ I(Nx)
        Reg = Σm.^2 + δ*I(Nx)
        Reg_inv = Reg \ I(Nx)
        Ā = Xp * Vtm' * Reg_inv * Um'
        FOO = eigen(Ā)
        W = FOO.vectors
        Λ = FOO.values
        Ω = log.(Λ) / dt
        Φ = Xp * Vtm' * Σm_inv * W

        return Φ, Ω, Ā, Nx
    end
end


"""
    choose_ro(Σ::Vector) → Vector

Choose reduced order (ro) that preserves an acceptable energy.

# Arguments
- `Σ`: Singular value vector from the SVD of some Hankel Matrix

# Returns
- `Vector`: ro values that preserve an acceptable energy
"""
function choose_ro(Σ::Vector)
    # Energy loss from truncation
    en = 1 .- sqrt.(cumsum(Σ .^ 2)) / norm(Σ)

    # loop through ROM sizes
    en_vals = map(x -> 10.0^x, -1.0:-1.0:-15.0)
    r_all = Vector{Float64}()
    for rr = axes(en_vals, 1)
        # determine # basis functions to retain based on energy lost
        en_thresh = en_vals[rr]
        push!(r_all, findfirst(x -> x < en_thresh, en))
    end

    return Int.(r_all)
end


"""
    normalize(val::Real, x_min::Real, x_max::Real, a::Real=-1, b::Real=1) → Real

Normalize `val` ∈ [x_min, x_max] to be in [a,b]

# Arguments
- `val::Real`: value to be normalized
- `x_min::Real`: minimum value of the range
- `x_max::Real`: maximum value of the range
- `a::Real`: minimum value of the normalized range
- `b::Real`: maximum value of the normalized range

# Returns
- `Real`: normalized value
"""
function normalize(val::Real, x_min::Real, x_max::Real, a::Real=-1, b::Real=1)
    return (b-a) * (val - x_min)/(x_max - x_min) - a
end


"""
    normalize(vec::Vector, xs_min::Vector, xs_max::Vector, a::Real=-1, b::Real=1) → Vector

Normalize values in vector `vec` ∈ [xs_min, xs_max] to be in [a,b]

# Arguments
- `vec::Vector`: vector to be normalized
- `xs_min::Vector`: minimum values of the range
- `xs_max::Vector`: maximum values of the range
- `a::Real`: minimum value of the normalized range
- `b::Real`: maximum value of the normalized range

# Returns
- `Vector`: normalized vector
"""
function normalize(vec::Vector, xs_min::Vector, xs_max::Vector, a::Real=-1, b::Real=1)
    nx = length(vec)
    return (b-a)*(vec - xs_min) ./ (xs_max - xs_min) + a*ones(nx)
end


"""
    normalize(mat::Matrix, lbs::Vector, ubs::Vector, a::Real=-1, b::Real=1) → Matrix

Normalize values in matrix `mat` ∈ [lbs, ubs] to be in [a,b]

# Arguments
- `mat::Matrix`: matrix to be normalized
- `lbs::Vector`: minimum values of the range
- `ubs::Vector`: maximum values of the range
- `a::Real`: minimum value of the normalized range
- `b::Real`: maximum value of the normalized range

# Returns
- `Matrix`: normalized matrix
"""
function normalize(mat::Matrix, lbs::Vector, ubs::Vector, a::Real=-1, b::Real=1)
    nx, N_sample = size(mat)
    return hcat(
        [normalize(
            mat[i,:],
            lbs[i]*ones(N_sample),
            ubs[i]*ones(N_sample),
            a, b
        ) for i = 1:nx]...
    )'
end


"""
    unnormalize(val::Real, x_min::Real, x_max::Real, a::Real=-1, b::Real=1) → Real

Get value back to its original range

# Arguments
- `val::Real`: value to be unnormalized
- `x_min::Real`: minimum value of the range
- `x_max::Real`: maximum value of the range
- `a::Real`: minimum value of the normalized range
- `b::Real`: maximum value of the normalized range  

# Returns
- `Real`: unnormalized value
"""
function unnormalize(val::Real, x_min::Real, x_max::Real, a::Real=-1,b::Real=1)
    return (val - a) * (x_max - x_min)/(b-a) + x_min
end


"""
    unnormalize(vec::Vector, xs_min::Vector, xs_max::Vector, a::Real=-1, b::Real=1) → Vector

Get vector back to its original range

# Arguments
- `vec::Vector`: vector to be unnormalized
- `xs_min::Vector`: minimum values of the range
- `xs_max::Vector`: maximum values of the range
- `a::Real`: minimum value of the normalized range
- `b::Real`: maximum value of the normalized range

# Returns
- `Vector`: unnormalized vector
"""
function unnormalize(vec::Vector, xs_min::Vector, xs_max::Vector, a::Real=-1,b::Real=1)
    return [unnormalize(vec[i], xs_min[i], xs_max[i], a, b) for i = 1:length(vec)]
end


"""
    lift_state_RBF(x, n_KO, Ms) → Vector

Compute lifted state using radial basis functions

# Arguments
- `x::Vector`: state vector
- `n_KO::Int`: number of radial basis functions
- `Ms::Matrix`: centers of radial basis functions

# Returns
- `Vector`: lifted state
"""
function lift_state_RBF(x::Vector, n_KO::Int, Ms::Matrix)
    ϕ = vcat(x, [norm(x - Ms[:,i])^2 * log(norm(x - Ms[:,i])) for i = 1:n_KO])
end


"""
    lift_data_RBF(Xm, Xp, n_KO, Ms, verbose=true) → Xm_lifted, Xp_lifted, Ms

Lift input data using RBF, implementation with known centers Ms

# Arguments
- `Xm::AbstractArray`: initial states, `n`-by-`N_sample` where `n` is the state dimension
- `Xp::AbstractArray`: final states (either x_{i+1} or dx/dt), `n`-by-`N_sample`
- `n_KO`: number of RBF basis to use
- `Ms::Matrix`: centers of radial basis functions
- `verbose::Bool`: verbosity flag for progress monitoring

# Returns
- `Xm_lifted::Matrix`: lifted initial states, `(n+n_KO)`-by-`N_sample`
- `Xp_lifted::Matrix`: lifted final states, `(n+n_KO)`-by-`N_sample`
- `Ms::Matrix`: centers of radial basis functions
"""
function lift_data_RBF(
    Xm::AbstractArray,
    Xp::AbstractArray,
    n_KO::Int,
    Ms::Matrix;
    verbose::Bool=true
)
    # lift states
    if verbose
        @info "Lifting states..."
    end
    nx, n_data = size(Xm)
    Xm_lifted = zeros(n_KO+nx, n_data)
    Xp_lifted = zeros(n_KO+nx, n_data)
    @showprogress for i = 1:n_data
        Xm_lifted[:,i] = lift_state_RBF(Xm[:,i], n_KO, Ms)
        Xp_lifted[:,i] = lift_state_RBF(Xp[:,i], n_KO, Ms)
    end
    return Xm_lifted, Xp_lifted, Ms
end


"""
    lift_data_RBF(Xm, Xp, n_KO, display=:iter, verbose=true) → Xm_lifted, Xp_lifted, Ms

Lift input data using RBF

# Arguments
- `Xm::AbstractArray`: initial states, `n`-by-`N_sample` where `n` is the state dimension
- `Xp::AbstractArray`: final states (either x_{i+1} or dx/dt), `n`-by-`N_sample`
- `n_KO`: number of RBF basis to use
- `display::Symbol`: verbosity for kmeans clustering
- `verbose::Bool`: verbosity flag for progress monitoring

# Returns
- `Xm_lifted::Matrix`: lifted initial states, `(n+n_KO)`-by-`N_sample`
- `Xp_lifted::Matrix`: lifted final states, `(n+n_KO)`-by-`N_sample`
- `Ms::Matrix`: centers of radial basis functions
"""
function lift_data_RBF(Xm::AbstractArray, Xp::AbstractArray, n_KO::Int; 
            display::Symbol=:iter, verbose::Bool=true)
    # Get clusters
    R = Clustering.kmeans(hcat(Xm,Xp), n_KO; maxiter=2000, display=display);
    Ms = R.centers    # get the cluster centers
    return lift_data_RBF(Xm, Xp, n_KO, Ms; verbose=verbose)
end
