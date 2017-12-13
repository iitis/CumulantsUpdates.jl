"""
  dataupdat(X::Matrix{T}, Xup::Matrix{T}) where T<:AbstractFloat

Returns Matrix{Float} of size(X), first u = size(Xup, 1) rows of X are removed and
at the end the updat Xup is appended.


```jldocstests

julia> a = ones(4,4);

julia> b = zeros(2,4);

julia> dataupdat(a,b)
4×4 Array{Float64,2}:
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0

```
"""

dataupdat(X::Matrix{T}, Xup::Matrix{T}) where T<:AbstractFloat =
  vcat(X,Xup)[1+size(Xup, 1):end, :]


"""

    momentupdat(M::SymmetricTensor{Float, N}, X::Matrix, Xup::Matrix)

Returns SymmetricTensor{Float, N} updated moment, given original moment, original data and update
of data - dataup

```jldocstests
julia> x = ones(6, 2);

julia> m = moment(x, 3);

julia> y = 2*ones(2,2);

julia> momentupdat(m, x, y)
SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]
[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

```
"""

function momentupdat(M::SymmetricTensor{T, m}, X::Matrix{T}, Xup::Matrix{T}) where {T<:AbstractFloat, m}
  tup = size(Xup,1)
  if tup == 0
    return M
  else
    return M + tup/size(X, 1)*(moment(Xup, m, M.bls) - moment(X[1:tup,:], m, M.bls))
  end
end

"""

  momentupdat(M::Vector{SymmetricTensor{T}}, X::Matrix{T}, Xup::Matrix{T})

  Returns Vector{SymmetricTensor} of updated moments

"""

momentupdat(M::Vector{SymmetricTensor{T}}, X::Matrix{T}, Xup::Matrix{T}) where T <: AbstractFloat =
    [momentupdat(M[i], X, Xup) for i in 1:length(M)]


"""
  momentarray(X::Matrix{Float}, m::Int, b::Int)

Returns an array of Symmetric Tensors of moments given data and maximum moment order
"""

momentarray(X::Matrix{T}, m::Int = 4, b::Int = 4) where T <: AbstractFloat =
    [moment(X, i, b) for i in 1:m]

"""

  moms2cums!(M::Vector{SymmetricTensor})

Changes vector of Symmetric Tensors of moments to vector of Symmetric Tensors of cumulants
```jldocstests

julia> m = momentarray(ones(20,3), 3);

julia> moms2cums!(m)

julia> m[3]

SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.0 0.0; 0.0 0.0]
[0.0 0.0; 0.0 0.0] #NULL; #NULL #NULL]
Nullable{Array{Float64,3}}[[0.0 0.0; 0.0 0.0] [0.0; 0.0]; #NULL [0.0]], 2, 2, 3, false)
```

"""


function moms2cums!(M::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  m = length(M)
  for i in 1:m
    f(σ::Int) = outerprodcum(i, σ, M...; exclpartlen = 0)
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


function cums2moms(cum::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  m = length(cum)
  Mvec = Array{SymmetricTensor{T}}(m)
  for i in 1:m
    f(σ::Int) = outerprodcum(i, σ, cum...; exclpartlen = 0)
    @inbounds Mvec[i] = cum[i]
    prods = pmap(f, [(2:m)...])
    for k in 2:m
      @inbounds Mvec[i] += prods[k-1]
    end
  end
  Mvec
end

"""
  cumulantscache(X::Matrix{T}, m::Int = 4, b::Int = 4)

Returns a a vector of SymmetricTensors of cumulants of order 1, ..., m calculated.
Moments array are saved in tmp/cumdata.jld for further cumulants updates,
parameter b is the block size for the block structure of cumulants.
"""

cumulantscache(X::Matrix{T}, m::Int = 4, b::Int = 4) where T <: AbstractFloat =
  cumulantsupdat(X, zeros(0, size(X,2)), m, b)[1]

"""
  cumulantsupdat(X::Matrix{T}, Xup::Matrix{T}, m::Int = 4, b::Int = 4)

Returns a vector of SymmetricTensors of cumulants of order 1, ..., m calculated
for X, updated by Xup. Loads moments arrays and saves updated moment
arrays in /tmp/cumdata.jld b - the block size for the block structure of cumulants.


```jldoctests
julia> x = ones(10,2);

julia> y = 2*ones(2,2);

julia> c = cumulantsupdat(x, 3, 2);

julia> cumulantsupdat(x, y, 3)

3-element Array{SymmetricTensor{Float64,N} where N,1}:
 SymmetricTensor{Float64,1}(Nullable{Array{Float64,1}}[[1.2, 1.2]], 2, 1, 2, true)
 SymmetricTensor{Float64,2}(Nullable{Array{Float64,2}}[[0.16 0.16; 0.16 0.16]], 2, 1, 2, true)
 SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.096 0.096; 0.096 0.096]
[0.096 0.096; 0.096 0.096]], 2, 1, 2, true)
```
"""

function cumulantsupdat(X::Matrix{T}, Xup::Matrix{T}, m::Int = 4, b::Int = 4) where T <: AbstractFloat
  cpath = "/tmp/cumdata.jld"
  d = try load(cpath) catch end
  X1 = try(d["x"]) catch end
  M = try(d["M"]) catch [1.] end
  Xp = dataupdat(X, Xup)
  if ((length(M) == m) & (X1 == X) & (typeof(M[1]) <: SymmetricTensor))
    println("update")
    Mup = momentupdat(M, X, Xup)
  else
    println("compute")
    Mup = momentarray(Xp, m, b)
  end
  try save(cpath, "M", Mup, "x", Xp) catch end
  moms2cums!(Mup)
  Mup, Xp
end
