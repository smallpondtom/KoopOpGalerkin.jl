export KoopmanSolve

"""
    KoopmanSolve(ord::Integer, fx::Array{Float64,3}, x0::AbstractArray{Float64}, 
        tspan::AbstractArray{Float64}) → AbstractArray{Float64}

Solve a nonlinear system using Koopman Operator with Orthogonal Polynomial basis.

# Arguments
- `ord::Integer`: Chosen order of the polynomial.
- `fx::Array{Float64,3}`: Function array including coefficients and exponent info.
- `x0::AbstractArray{Float64}`: State initial conditions.
- `tspan::AbstractArray{Float64}`: Time span of the simulation.

# Returns
`AbstractArray{Float64}`: The solution as a function of time.
"""
function KoopmanSolve(ord::Integer, fx::Array{Float64,3}, x0::AbstractArray{Float64}, 
            tspan::AbstractArray{Float64})::AbstractArray{Float64}
    Nx, Nt, _ = size(fx)
    @assert ord >= 3 "Order of basis should be greater than 3 for this application."
    
    # Number of basis functions
    ns = numOfBasis(ord, Nx)
    
    # Legendre basis indexes
    ind = basisIndex(ord, ns, Nx)
    
    # Definition of the Legendre polynomials
    LPC = legendrePoly(ord)
    
    # Legendre polynomials in multiple dimensions
    NLPC = normalizedLegendrePoly(ord, LPC)
    MNLPC = multiLegendrePoly(ns, Nx, ind, NLPC)
    
    # Derivative of Normalized Legendre polynomials
    DNLPC = legendreDiff(ord, NLPC)
    
    # Operator Matrix
    K = KOM(ns, Nx, Nt, fx, ind, NLPC, DNLPC, MNLPC)
    
    # Observables
    H = observables(ns, Nx, ind, MNLPC)
    
    # Eigenvalue decomposition
    D, V = eigen(K)  # eigenvalues and right eigenvectors
    iV = V \ 1.0I(ns)

    # Initial condition mapped onto Legendre
    L0 = funcSpaceIC(ns, x0, ind, MNLPC)
    
    # Computation of the solution as a function of time
    nt = length(tspan)
    Sol = zeros(Float64, Nx, nt)

    for k in 1:nt
        Sol[:, k] = real.(H * V * exp(Diagonal(D) * tspan[k]) * iV * L0)
    end

    return Sol
end



"""
    KOM(ns::Integer, Nx::Integer, Nt::Integer, fx::Array{Float64,3}, ind::Matrix{<:Integer}, 
             NLPC::AbstractArray{Float64}, DNLPC::AbstractArray{Float64}, 
             MNLPC::AbstractArray{Float64}) → AbstractArray{Float64}

Compute the Koopman Operator Matrix (KOM).

# Arguments
- `ns::Integer`: Total number of basis.
- `Nx::Integer`: Total number of variables.
- `Nt::Integer`: Total number of terms in f(x).
- `fx::Array{Float64,3}`: Function array including coefficients and exponent info.
- `ind::Matrix{<:Integer}`: Legendre basis index.
- `NLPC::AbstractArray{Float64}`: Coefficients of normalized Legendre polynomials.
- `DNLPC::AbstractArray{Float64}`: Coefficients of derivative of normalized Legendre polynomials.
- `MNLPC::AbstractArray{Float64}`: Coefficients of multivariate normalized Legendre polynomials.

# Returns
`AbstractArray{Float64}`: The Koopman Operator Matrix.
"""
function KOM(ns::Integer, Nx::Integer, Nt::Integer, fx::Array{Float64,3}, ind::Matrix{<:Integer}, 
                  NLPC::AbstractArray{Float64}, DNLPC::AbstractArray{Float64}, 
                  MNLPC::AbstractArray{Float64})::AbstractArray{Float64}
    K = zeros(Float64, ns, ns)
    MDLP = zeros(Float64, Nx, ns)
    par = zeros(Integer, Nx)

    for i in 1:ns
        for j in 1:ns
            for k in 1:Nx
                MDLP[k, j] = DNLPC[ind[i, k] + 1, ind[j, k] + 1]
                for l in 1:Nx
                    if k != l
                        MDLP[k, j] *= NLPC[ind[i, l] + 1, ind[j, l] + 1]
                    end
                end
            end
        end

        DB = zeros(Float64, ns * 2, Nx + 1)
        s = 0
        for dim in 1:Nx
            for ifx in 1:Nt
                if fx[dim, ifx, 1] == 0
                    break
                else
                    for j in 1:ns
                        if MDLP[dim, j] != 0
                            s += 1
                            DB[s, 1] = fx[dim, ifx, 1] * MDLP[dim, j]
                            for k in 2:(Nx + 1)
                                DB[s, k] = fx[dim, ifx, k] + ind[j, k - 1]
                            end
                        end
                    end
                end
            end
        end

        for k in 1:ns
            for j in 1:ns
                if MNLPC[k, j] != 0
                    for ifx in 1:s
                        flag = 1
                        for iN in 1:Nx
                            par[iN] = round(Int, ind[j, iN] + DB[ifx, iN + 1] + 1)
                            if par[iN] % 2 == 0
                                flag = 0
                                break
                            end
                        end
                        if flag == 1
                            K[i, k] += (2^Nx) * MNLPC[k, j] * DB[ifx, 1] / prod(par)
                        end
                    end
                end
            end
        end
    end

    return K
end


"""
    observables(ns::Integer, Nx::Integer, ind::Matrix{<:Integer}, 
                    MLPC::AbstractArray{Float64}) → AbstractArray{Float64}

Calculate the observables for a given set of basis functions and multi-variate Legendre polynomial coefficients.

# Arguments
- `ns::Integer`: Total number of basis functions.
- `Nx::Integer`: Total number of variables.
- `ind::Matrix{<:Integer}`: Legendre basis index.
- `MNLPC::AbstractArray{Float64}`: Coefficients of multivariate Legendre polynomial.

# Returns
`AbstractArray{Float64}`: The observables matrix.
"""
function observables(ns::Integer, Nx::Integer, ind::Matrix{<:Integer}, 
        MNLPC::AbstractArray{Float64})::AbstractArray{Float64}

    H = zeros(Float64, Nx, ns)
    par = zeros(Integer, Nx)

    for i in 1:Nx
        # Observable Polynomials
        Obs = zeros(Float64, 1, Nx + 1)
        Obs[1, 1] = 1
        Obs[1, i + 1] = 1
    
        # Matrix integration
        for k in 1:ns
            for j in 1:ns
                if MNLPC[k, j] != 0
                    flag = 1
                    for iN in 1:Nx
                        par[iN] = round(Integer, ind[j, iN] + Obs[1, iN + 1] + 1)
                        if par[iN] % 2 == 0
                            flag = 0
                            break
                        end
                    end
                    if flag == 1
                        H[i, k] += (2^Nx) * MNLPC[k, j] * Obs[1, 1] / prod(par)
                    end
                end
            end
        end
    end

    return H
end


"""
    funcSpaceIC(ns::Integer, x0::AbstractArray{Float64}, ind::Matrix{<:Integer}, 
        MNLPC::AbstractArray{Float64}) → AbstractArray{Float64}

Generate the functional space initial conditions.

# Arguments
- `ns::Integer`: Total number of basis functions.
- `x0::AbstractArray{Float64}`: Initial condition of states, reshaped into a row vector if necessary.
- `ind::Matrix{<:Integer}`: Legendre polynomial index.
- `MNLPC::AbstractArray{Float64}`: Coefficients of multivariate Legendre polynomial.

# Returns
`AbstractArray{Float64}`: The functional space initial conditions.
"""
function funcSpaceIC(ns::Integer, x0::AbstractArray{Float64}, ind::Matrix{<:Integer}, 
            MNLPC::AbstractArray{Float64})::AbstractArray{Float64}
    g0 = zeros(Float64, ns)
    for i in 1:ns, j in 1:ns
        g0[i] += MNLPC[i, j] * prod(x0 .^ ind[j, :])
    end
    return g0
end
