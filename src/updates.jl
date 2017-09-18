"""

    momentupdat(M::SymmetricTensor{Float, N}, X::Matrix, Xup::Matrix)

Returns SymmetricTensor{Float, N} updated moment, given original moment, original data and update
of data - dataup
"""

function momentupdat(M::SymmetricTensor{T, N}, X::Matrix{T}, Xup::Matrix{T}) where T<:AbstractFloat where N
  tup = size(Xup,1)
  M + tup/size(X, 1)*(moment(Xup, N, M.bls) - moment(X[1:tup,:], N, M.bls))
end

"""
  momentarray(X::Matrix{Float}, m::Int, b::Int)

Returns an array of Symmetric Tensors of moments given data and maximum moment order
"""

momentarray(X::Matrix{T}, m::Int = 4, b::Int = 2) where T <: AbstractFloat =
    [moment(X, i, b) for i in 1:m]

"""

  moms2cums!(M::Vector{SymmetricTensor})

Changes vector of Symmetric Tensors of moments to vector of Symmetric Tensors of cumulants
"""

function moms2cums1!(M::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  for i in 1:length(M)
    for sigma in 2:i
      @inbounds M[i] -= outerprodcum(i, sigma, M...; exclpartlen = 0)
    end
  end
end


function moms2cums!(M::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  m = length(M)
  for i in 1:m
    f(sigma::Int) = outerprodcum(i, sigma, M...; exclpartlen = 0)
    prods = pmap(f, [(2:m)...])
    for k in 2:m
      @inbounds M[i] -= prods[k-1]
    end
  end
end

"""

  cums2moms(cum::Vector{SymmetricTensor})

Returns vector of Symmetric Tensors of moments given vector of Symmetric Tensors
of cumulants
"""

function cums2moms1(cum::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  m = length(cum)
  Mvec = Array{SymmetricTensor{T}}(m)
  for i in 1:m
    Mvec[i] = cum[i]
    for sigma in 2:i
      Mvec[i] += outerprodcum(i, sigma, cum...; exclpartlen = 0)
    end
  end
  Mvec
end


function cums2moms(cum::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  m = length(cum)
  Mvec = Array{SymmetricTensor{T}}(m)
  for i in 1:m
    f(sigma::Int) = outerprodcum(i, sigma, cum...; exclpartlen = 0)
    @inbounds Mvec[i] = cum[i]
    prods = pmap(f, [(2:m)...])
    for k in 2:m
      @inbounds Mvec[i] += prods[k-1]
    end
  end
  Mvec
end


"""

    cumulantsupdat(cum::Vector{SymmetricTensor}, X::Matrix, Xup::Matrix)

Returns Vector{SymmetricTensor} of updated cumulants

"""

function cumulantsupdat(cum::Vector{SymmetricTensor{T}}, X::Matrix{T}, Xup::Matrix{T}) where T <: AbstractFloat
  M = cums2moms(cum)
  @inbounds Mup = [momentupdat(M[i], X, Xup) for i in 1:length(cum)]
  moms2cums!(Mup)
  Mup
end
