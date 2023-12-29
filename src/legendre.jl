export legendreDiff, multiLegendrePoly, normalizedLegendrePoly, legendrePoly

"""
    DNLP(c::Integer, NLPC::AbstractArray{Float64}) -> AbstractArray{Float64}

Compute the coefficients of the Derivative of the Normalized Legendre Polynomials (DNLP).

# Arguments
- `c::Integer`: Chosen order of the polynomial.
- `NLPC::AbstractArray{Float64}`: Coefficients of normalized Legendre polynomial.

# Returns
`Matrix{Float64}`: Coefficients of the derivative of the normalized Legendre polynomials.
"""
function legendreDiff(c::Integer, NLPC::AbstractArray{Float64})::AbstractArray{Float64}
    DNLPC = zeros(Float64, c + 1, c + 1)

    for i = 2:(c + 1)
        for j = 1:c
            DNLPC[i, j] = j * NLPC[i, j + 1]
        end
    end

    return DNLPC
end


"""
    MLP(ns::Int, Nx::Int, ind::Matrix{<:Integer}, NLPC::Matrix{Float64}) → Matrix{Float64}

Compute the coefficients of the Multivariate Normalized Legendre polynomials (MNLP).

# Arguments
- `ns::Int`: Total number of basis.
- `Nx::Int`: Total number of variables.
- `ind::Matrix{<:Integer}`: Legendre Polynomial index.
- `NLPC::AbstractArray{Float64}`: Coefficients of normalized Legendre polynomial.

# Returns
`AbstractArray{Float64}`: Coefficients of the multivariate Legendre polynomials.
"""
function multiLegendrePoly(ns::Integer, Nx::Integer, 
        ind::Matrix{<:Integer}, NLPC::AbstractArray{Float64})::AbstractArray{Float64}
    MNLPC = ones(Float64, ns, ns)

    for i = 1:ns
        for j = 1:ns
            for k = 1:Nx
                MNLPC[i,j] *= NLPC[ind[i,k]+1,ind[j,k]+1]
            end
        end
    end

    return MNLPC
end


"""
    normalizedLegendrePoly(c::Integer, LPC::AbstractArray{Float64}) → AbstractArray{Float64}

Return the normalized Legendre polynomials up to degree `c`.

# Arguments
- `c::Integer`: The degree of the polynomial.
- `LPC::AbstractArray{Float64}`: The Legendre polynomial coefficients.

# Returns
- `AbstractArray{Float64}`: The normalized Legendre polynomials.
"""
function normalizedLegendrePoly(c::Integer, LPC::AbstractArray{Float64})::AbstractArray{Float64}
    NLPC = zeros(c + 1, c + 1)

    for i in 1:(c + 1)
        for j in 1:(c + 1)
            NLPC[i,j] = sqrt((2 * (i - 1) + 1) / 2) * LPC[i,j]
        end
    end

    return NLPC
end


"""
    legendrePoly(c::Integer) → AbstractArray{Float64}

Return the Legendre polynomials up to degree `c`.

# Arguments
- `c::Integer`: The degree of the polynomial.

# Returns
- `AbstractArray{Float64}`: The Legendre polynomials.
"""
function legendrePoly(c::Integer)::AbstractArray{Float64}
    # Initialize the Legendre polynomial coefficient matrix
    LPC = zeros(c + 1, c + 1)
    LPC[1, 1] = 1
    LPC[2, 2] = 1

    for i in 3:(c + 1)
        for j in 1:(i - 1)
            LPC[i, j + 1] += (2 * (i - 2) + 1) / (i - 1) * LPC[i - 1, j]
            LPC[i, j] -= (i - 2) / (i - 1) * LPC[i - 2, j]
        end
    end

    return LPC
end


"""
    basisIndex(c::Integer, ns::Integer, Nx::Integer) → Matrix{<:Integer}

Return the indices of the basis functions for a polynomial of degree `c` in `n` dimensions.

# Arguments
- `c::Integer`: The degree of the polynomial.
- `n::Integer`: The number of dimensions.
- `Nx::Integer`: The number of points in each dimension.

# Returns
- `Matrix{<:Integer}`: The indices of the basis functions.
"""
function basisIndex(c::Integer, ns::Integer, Nx::Integer)::Matrix{<:Integer}
    ind = zeros(Int, ns, Nx)
    s = 1
    for ord in 0:c
        ind, s = iterateDims(ord, zeros(Int, 1, Nx), 1, c, Nx, ind, s)
    end
    return ind
end


"""
    iterateDims(ord::Integer, currentComb::AbstractArray, currentDim::Integer, 
        c::Integer, Nx::Integer, ind::Matrix{<:Integer}, s::Integer) → Tuple{AbstractArray{Integer}, Integer}

Iterate over the dimensions of the polynomial.

# Arguments
- `ord::Integer`: The order of the polynomial.
- `currentComb::Array`: The current combination of indices.
- `currentDim::Integer`: The current dimension.
- `c::Integer`: The degree of the polynomial.
- `Nx::Integer`: The number of points in each dimension.
- `ind::Matrix{<:Integer}`: The indices of the basis functions.
- `s::Integer`: The current index of the basis function.

# Returns
- `AbstractArray{Integer}`: The indices of the basis functions.
- `Integer`: The current index of the basis function.
"""
function iterateDims(ord::Integer, currentComb::Array, currentDim::Integer, 
    c::Integer, Nx::Integer, ind::Matrix{<:Integer}, s::Integer)::Tuple{Matrix{<:Integer}, Integer}

    if currentDim == Nx
        currentComb[Nx] = ord - sum(currentComb[1:Nx-1])
        if sum(currentComb) == ord
            ind[s, :] = currentComb
            s += 1
        end
        return ind, s
    end

    for i in 0:ord
        newComb = deepcopy(currentComb)
        newComb[currentDim] = i
        if sum(newComb) <= ord
            ind, s = iterateDims(ord, newComb, currentDim + 1, c, Nx, ind, s)
        end
    end
    return ind, s
end


"""
    numOfBasis(c::Integer, n::Integer) → Integer

Return the number of basis functions for a polynomial of degree `c` in `n` dimensions.

# Examples
```julia-repl
julia> using KoopOpGalerkin

julia> numOfBasis(0,0)

julia> numOfBasis(4,2)

julia> numOfBasis(2,4)
```

# Arguments
- `c::Integer`: The degree of the polynomial.
- `n::Integer`: The number of dimensions.

# Returns
- `Integer`: The number of basis functions.
"""
function numOfBasis(c::Integer, n::Integer)
    return binomial(n + c, c)
end
