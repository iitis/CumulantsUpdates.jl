[![Build Status](https://travis-ci.org/ZKSI/CumulantsUpdates.jl.svg?branch=master)](https://travis-ci.org/ZKSI/CumulantsUpdates.jl)
[![Coverage Status](https://coveralls.io/repos/github/ZKSI/CumulantsUpdates.jl/badge.svg?branch=master)](https://coveralls.io/github/ZKSI/CumulantsUpdates.jl?branch=master)

# CumulantsUpdates.jl

Updates following statistics of `n`-variate data
* `m`'th moment tensor,
* an array of moment tensors of order `1,2,...,m`.

Given `t` realisations of `n`-variate data: `X ∈ ℜᵗˣⁿ`, and the update `Xᵤₚ ∈ ℜᵘˣⁿ`
returns array of updated cumulant tensors of order `1,2,...,m`.

To store symmetric tensors uses a `SymmetricTensors` type, requires [SymmetricTensors.jl](https://github.com/ZKSI/SymmetricTensors.jl). To convert to array, run:

```julia
julia> convert(Array, st::SymmetricTensors{T, m})
```
to convert back, run:

```julia
julia>  convert(SymmetricTensor, a::Array{T,m})
```
Requires [Cumulants.jl](https://github.com/ZKSI/Cumulants.jl).

As of 01/01/2017 [kdomino](https://github.com/kdomino) is the lead maintainer of this package.

## Installation

Within Julia, run

```julia
julia> Pkg.clone("https://github.com/ZKSI/CumulantsUpdates.jl")
```

to install the files.  Julia 0.6 is required.

## Parallel computation

For parallel computation just run
```julia
julia> addprocs(n)
julia> @everywhere using CumulantsUpdates
```

## Functions

### Moment update

To update the moment tensor `M::SymmetricTensor{T, m}` for `m ≤ 1` given  data `X ∈ ℜᵗˣⁿ` and the update `Xᵤₚ ∈ ℜᵘˣⁿ` where `u < t` run

```julia
julia> momentupdat(M::SymmetricTensor{T, m}, X::Matrix{T}, Xᵤₚ::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns a `SymmetricTensor{T, m}` of the moment tensor of updated multivariate data: `X' = vcat(X,Xᵤₚ)[1+u:end, :] ∈ Rᵗˣⁿ`. The output of `momentupdat(M, X, Xᵤₚ) = moment(X', m)`,  if `u ≪ t` `momentupdat()` is much faster than a recalculation.

```julia
julia> x = ones(6, 2);

julia> m = moment(x, 3)
SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[1.0 1.0; 1.0 1.0]

[1.0 1.0; 1.0 1.0]],2,1,2,true)

julia> y = 2*ones(2,2);

julia> momentupdat(m, x, y)
SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

julia> moment(vcat(x,y)[3:end,:],3)
SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

```

Function


```julia
julia> momentupdat(M::Array{SymmetricTensor{T, m}}, X::Matrix{T}, Xᵤₚ::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns an `Array{SymmetricTensor{T, m}}` of moment tensors of order `1, ..., m` of updated multivariate data: `X' = vcat(X,Xᵤₚ)[1+u:end, :] ∈ Rᵗˣⁿ`. The output of `momentupdat(M, X, Xᵤₚ) = arraymoment(X', m)`,  if `u ≪ t` `momentupdat()` is much faster than a recalculation.

### Cumulants update

To update the vector of cumulants `cums = [SymmetricTensor{T, 1}, SymmetricTensor{T, 2}, ...,SymmetricTensor{T, m}]` of order `1, ..., m` given original data `X \in R^{t, n}` and the update `Xup \in R^{tup, n}` for `t > tup` just run:

```julia
julia> cumulantsupdat(cums::Vector{SymmetricTensor{T}}, X::Matrix{T}, Xup::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns a vector `[SymmetricTensor{T, 1}, SymmetricTensor{T, 2}, ...,SymmetricTensor{T, m}]` of cumulant tensors of order `1,2,...,m` of updated multivariate data `Xprim = vcat(X,Xup)[1+tup:end, :] \in R^{t, n}`. If `Cumulants.cumulants(X, m) = cums` then the oputpu of `cumulantsupdat(cums, X, Xup)` will correspond to the output of `Cumulants.cumulants(Xprim, m)`. However if `t` is large and `t >> tup` `cumulantsupdat` is much faster. If the input `cums` is not a sequence of cumulants of order `1,2, ..., m` wrror will be returned. If cumulants in `cums` are not computed for the same data, results are meaningless.

```julia
julia> x = ones(10,2);

julia> y = 2*ones(2,2);

julia> c = cumulants(x, 3);

julia> cumulantsupdat(c, x, y)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N},1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Nullable{Array{Float64,1}}[[1.2,1.2]],2,1,2,true)                                             
 SymmetricTensors.SymmetricTensor{Float64,2}(Nullable{Array{Float64,2}}[[0.16 0.16; 0.16 0.16]],2,1,2,true)                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.096 0.096; 0.096 0.096]
[0.096 0.096; 0.096 0.096]],2,1,2,true)

julia> cumulants(vcat(x,y)[3:end, :], 3)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N},1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Nullable{Array{Float64,1}}[[1.2,1.2]],2,1,2,true)                                             
 SymmetricTensors.SymmetricTensor{Float64,2}(Nullable{Array{Float64,2}}[[0.16 0.16; 0.16 0.16]],2,1,2,true)                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.096 0.096; 0.096 0.096]
[0.096 0.096; 0.096 0.096]],2,1,2,true)

```

### Vector norm

```julia
julia> vecnorm(st::SymmetricTensor{T, m}; k::Union{Float64, Int}) where {T<:AbstractFloat, m}
```

Returns a vector norm of the `SymmetricTensors` type if `k != 0` The output of `vecnorm(st, k)` corresponds to the output of `vecnorn(convert(Array, st),k)`. However
`vecnorm(st::SymmetricTensor...)` uses the block structure implemented in `SymmetricTensors` and decreases the computer memory requirement.

```julia
julia> te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154];

julia> st = convert(SymmetricTensor, (reshape(te, (2,2,2))));

julia> vecnorm(st)
0.5273572868359742

julia> vecnorm(st, 2.5)
0.4468668679541424

julia> vecnorm(st, 1)
1.339089
```

### Tensor norms of cumulants

To compute tensor norms of cumulants of `X \in R^{t, n}` run:

```julia
julia>  cumnorms(X::Matrix{T}, m::Int = 4, norm::Bool = true, k::Union{Float64, Int}=2, b::Int = 3, cache::Bool = true) where T <: AbstractFloat
```

```julia
julia> X = ones(10,3);

julia> cumnorms(X, 3, false)
3-element Array{Float64,1}:
 1.73205
 0.0
 0.0

julia> cumnorms(X, 3, false, 1)
3-element Array{Float64,1}:
 3.0
 0.0
 0.0
```

Returns an array of Floats of `k` norms of `1, ..., m` cumulant tensors of `X` i.e. `||C_m|| = (\sum_{c \in C_m(X)} c^k)^(1/k)`,
* if `cache` cumulants are saved in a `/tmp/cumdata.jld` directory;
* if `norm = true` tensor norms for `m > 2` are normalised by `||C_2||^(m/2)` i.e. `h_m = ||C_m||/(||C_2||^(m/2))`,
* the parameter `b` is a block size, in the block structure of cumulants - see `SymmetricTensors` docs.

### Tensor norms of cumulants after the update

Given data, `X \in R^{t,n}`, the update `Xup \in R^{tup,n}` we can compute tensor norms for `Xprim = vcat(X,Xup)[1+tup:end, :] \in R^{t, n}` by running:

```julia
julia> cumupdatnorms(X::Matrix{T}, Xup::Matrix{T}, m::Int = 4, norm::Bool = true, k::Union{Float64, Int}=2, b::Int = 3, cache::Bool = true) where T <: AbstractFloat
```

```julia
julia> X = ones(10,3);

julia> Xup = zeros(5,3);

julia> cumupdatnorms(X, Xup, 3, false)[1]
3-element Array{Float64,1}:
 0.866025
 0.75
 0.0

julia> cumupdatnorms(X, Xup, 3, false)[2]
10×3 Array{Float64,2}:
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

```

Returns an array of Floats of `k` norms of updated `1, ..., m` cumulant tensors, and matrix of updated data: `cumupdatnorms(X, Xup)[2] = Xprim`
Parameters `m`, `norm`, `k`, `b` are as in `cumnorms`. If `cache` updated cumulants are saved in a `/tmp/cumdata.jld` directory. The output of `cumupdatnorms(X, Xup, m, norm, k)[1]` is the same as the output of `cumnorms(X, m, norm, k)`, however if cumulants of `X` were calculated and saved in `/tmp/cumdata.jld` `cumupdatnorms` will use fast `cumulantsupdat`.


# Performance analysis

To analyse the computational time of moments and cumulants updates vs `Cumulants.moment` or `Cumulants.cumulants` recalculation, we supply the executable script `comptimes.jl`.
The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's maksimum order, by default `m = 4`,
* `-n (vararg Int)`: numbers of marginal variables, by default `m = 15, 20, 25, 30`,
* `-t (Int)`: number of realisations of random variable, by default `t = 500000`,
* `-u (vararg Int)`: number of realisations of update, by default `u = 20000, 30000, 40000, 50000`,
* `-b (Int)`: blocks size, by default `b = 3`,
* `-p (Int)`: numbers of processes, by default `p = 1`.

To analyse the computational time of cumulants updates for different block sizes `1 < b =< Int(sqrt(n))`, we supply the executable script `comptimesblocks.jl`.
The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (Int)`: numbers of marginal variables, by default `m = 48`,
* `-u (vararg Int)`: number of realisations of the update, by default `u = 20000, 40000`.
* `-p (Int)`: numbers of processes, by default `p = 1`.

To plot computional times run executable `res/plotcomptimes.jl` on chosen `*.jld` file.


# Citing this work

Krzysztof Domino, Piotr Gawron, *On updates of high order cumulant tensors*, [arXiv:1701.06446](https://arxiv.org/abs/1701.06446)

This project was partially financed by the National Science Centre, Poland – project number 2014/15/B/ST6/05204.
