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

### Data update

To update data `X ∈ ℜᵗˣⁿ` by the update `Xᵤₚ ∈ ℜᵘˣⁿ` in the observation
window of size `t` and get `ℜᵗˣⁿ ∋ X' = vcat(X,Xᵤₚ)[1+u:end, :]` run:

```julia
julia> dataupdat(X::Matrix{T}, Xup::Matrix{T}) where T<:AbstractFloat
```
the condition `u < t` must be fulfilled.

```julia
julia> a = ones(4,4)
4×4 Array{Float64,2}:
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0

julia> b = zeros(2,4)
2×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0

julia> dataupdat(a,b)
4×4 Array{Float64,2}:
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```

### Moment update

To update the moment tensor `M::SymmetricTensor{T, m}` for  data `X ∈ ℜᵗˣⁿ`, given the update `Xᵤₚ ∈ ℜᵘˣⁿ` where `u < t` run

```julia
julia> momentupdat(M::SymmetricTensor{T, m}, X::Matrix{T}, Xᵤₚ::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns a `SymmetricTensor{T, m}` of the moment tensor of updated multivariate data:
`ℜᵗˣⁿ ∈ X' = dataupdat(X,Xᵤₚ)`, i.e.: `M = moment(X, m)`, `momentupdat(M, X, Xᵤₚ) = moment(X', m)`. If `u ≪ t` `momentupdat()` is much faster than a recalculation.

```julia
julia> x = ones(6, 2);

julia> julia> m = moment(x, 3)
SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[1.0 1.0; 1.0 1.0]

[1.0 1.0; 1.0 1.0]], 2, 1, 2, true)

julia> y = 2*ones(2,2);

julia> momentupdat(m, x, y)
SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]], 2, 1, 2, true)

julia> moment(vcat(x,y)[3:end,:],3)
SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]], 2, 1, 2, true)
```

Function `momentarray(X, m)` can be used to compute an array of moments of order `1, ..., m`
of data `X ∈ ℜᵗˣⁿ`

```julia
julia> momentarray(X::Matrix{T}, m::Int, b::Int) where {T<:AbstractFloat, m}
```
`b` - is a `SymmetricTensor` parameter, a block size.

One can update an array of moments by calling:

```julia
julia> momentupdat(M::Array{SymmetricTensor{T, m}}, X::Matrix{T}, Xᵤₚ::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns an `Array{SymmetricTensor{T, m}}` of moment tensors of order `1, ..., m` of updated multivariate data: `ℜᵗˣⁿ ∋ X' = dataupdat(X,Xᵤₚ)`, i.e. `Mₐᵣ = momentarray(X, m)`, `momentupdat(Mₐᵣ, X, Xᵤₚ) = momentarray(X', m)`.  

### Cumulants update

Presented functions are design for sequent update of multivariate cumulant tensors.
Hence it can be applied in a data streaming scheme. Suppose one has data `X ∈ ℜᵗˣⁿ`
and subsequent coming updates `Xᵤₚ ∈ ℜᵘˣⁿ` such that `u ≪ t`. Suppose one want to compute cumulant tensors in an observation window of size `t` each time the update comes.
One can first compute cumulants of order `1, ..., m` and cache moments of order `1, ..., m`
for further fast cumulants updates.

```julia
julia> cumulantscache(X::Matrix{T}, m::Int = 4, b::Int = 4) where {T<:AbstractFloat, m}
```

```julia
julia> X = ones(10,3);

julia> cumulantscache(X, 2,2)[1]

SymmetricTensors.SymmetricTensor{Float64,1}(Union{Array{Float64,1}, Void}[[1.0, 1.0], [1.0]], 2, 2, 3, false)                     
```

The function returns a vector `[SymmetricTensor{T, 1}, ...,SymmetricTensor{T, m}]` of cumulant tensors of order `1,...,m` of `X` and caches the array of moments of `X` in `tmp/cumdata.jld` for further fast cumulants updates.


If the subsequent update `Xᵤₚ ∈ ℜᵘˣⁿ` of `X ∈ ℜᵗˣⁿ` is recorded, to compute
cumulants vector `c = [C₁(X'), ..., Cₘ(X')]` of updated `ℜᵗˣⁿ ∋ X'ₖ = dataupdat(Xₖ, Xᵤₚ_ₖ)` run:

```julia
julia> c, X' = cumulantsupdat(X::Matrix{T}, Xᵤₚ::Matrix{T}, m::Int = 4, b::Int = 4) where {T<:AbstractFloat, m}
```

The function caches the array of updated moments in `tmp/cumdata.jld`, for next update cycle. Now set `X = X'` and wait for another update. If `u ≪ t`  `cumulantsupdat()` is much faster that a cumulants recalculation using `Cumulants.jl`. If `cumulantscache()` is not used first run of `cumulantsupdat()` would compute whole moments tensors in a classical
way and cache updated moments for another fast round.

```julia
julia> x = ones(10,2);

julia> y = 2*ones(2,2);

julia> cumulantscache(x, 3, 2)

julia> c, X' = cumulantsupdat(x, y, 3, 2)

julia> c

3-element Array{SymmetricTensors.SymmetricTensor{Float64,N} where N,1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Union{Array{Float64,1}, Void}[[1.2, 1.2]], 2, 1, 2, true)                                            
 SymmetricTensors.SymmetricTensor{Float64,2}(Union{Array{Float64,2}, Void}[[0.16 0.16; 0.16 0.16]], 2, 1, 2, true)                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[0.096 0.096; 0.096 0.096]
[0.096 0.096; 0.096 0.096]], 2, 1, 2, true)

julia> X'
10×2 Array{Float64,2}:
 1.0  1.0
 1.0  1.0
 1.0  1.0
 1.0  1.0
 1.0  1.0
 1.0  1.0
 1.0  1.0
 1.0  1.0
 2.0  2.0
 2.0  2.0
```

Record that `cumulantsupdat(x, y, m, b) = cumulants(dataupdat(x, y), m, b)`, but
if `u ≪ t` `cumulantsupdat()` is much faster. The update scheme where `X ∈ ℜᵗˣⁿ`
original data, `U` - array of updates with elements `Xᵤₚ_ᵢ ∈ ℜᵘⁱ ˣ ⁿ`

```julia
julia> Xp = X;

julia> c = cumulantscache(X, m, b)

julia>f(c)

julia> for Xᵤₚ in U
        c, Xp =cumulantsupdat(Xp, Xᵤₚ, m, b)
        f(c)
       end

```
Here `f(c)` is some function that extracts information form cumulants, for presentation
purpose we can use `f(c) = println(map(vecnorm, c))`

### Vector norm

```julia
julia> vecnorm(st::SymmetricTensor{T, m}, p::Union{Float64, Int}) where {T<:AbstractFloat, m}
```

Returns a `p`-norm of the tensor stored as `SymmetricTensors`, supported for `k ≠ 0`. The output of `vecnorm(st, p) = vecnorn(convert(Array, st),p)`. However
`vecnorm(st::SymmetricTensor, p)` uses the block structure implemented in `SymmetricTensors`, hence is faster and decreases the computer memory requirement.

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
### Convert cumulants to moments and moments to cumulants

Given `M` a vector of moments of order `1, ..., m` to change it to a vector
of cumulants of order `1, ..., m` using

```julia
julia> function moms2cums!(M::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
```
One can convert a vector of cumulants `c` to a vector of moments by running

```julia
julia> function cums2moms(c::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
```

```julia

julia> m = momentarray(ones(20,3), 3, 2)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N} where N,1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Union{Array{Float64,1}, Void}[[1.0, 1.0], [1.0]], 2, 2, 3, false)                                                                                                                           
 SymmetricTensors.SymmetricTensor{Float64,2}(Union{Array{Float64,2}, Void}[[1.0 1.0; 1.0 1.0] [1.0; 1.0]; nothing [1.0]], 2, 2, 3, false)                                                                                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[1.0 1.0; 1.0 1.0]
[1.0 1.0; 1.0 1.0] nothing; nothing nothing]
Union{Array{Float64,3}, Void}[[1.0 1.0; 1.0 1.0] [1.0; 1.0]; nothing [1.0]], 2, 2, 3, false)


julia> moms2cums!(m)

julia> m
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N} where N,1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Union{Array{Float64,1}, Void}[[1.0, 1.0], [1.0]], 2, 2, 3, false)                                                                                                                       
 SymmetricTensors.SymmetricTensor{Float64,2}(Union{Array{Float64,2}, Void}[[0.0 0.0; 0.0 0.0] [0.0; 0.0]; #undef [0.0]], 2, 2, 3, false)                                                                                             
 SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[0.0 0.0; 0.0 0.0]
[0.0 0.0; 0.0 0.0] #undef; #undef #undef]
Union{Array{Float64,3}, Void}[[0.0 0.0; 0.0 0.0] [0.0; 0.0]; #undef [0.0]], 2, 2, 3, false)


julia>  cums2moms(m)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N} where N,1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Union{Array{Float64,1}, Void}[[1.0, 1.0], [1.0]], 2, 2, 3, false)                                                                                                                       
 SymmetricTensors.SymmetricTensor{Float64,2}(Union{Array{Float64,2}, Void}[[1.0 1.0; 1.0 1.0] [1.0; 1.0]; #undef [1.0]], 2, 2, 3, false)                                                                                             
 SymmetricTensors.SymmetricTensor{Float64,3}(Union{Array{Float64,3}, Void}[[1.0 1.0; 1.0 1.0]
[1.0 1.0; 1.0 1.0] #undef; #undef #undef]
Union{Array{Float64,3}, Void}[[1.0 1.0; 1.0 1.0] [1.0; 1.0]; #undef [1.0]], 2, 2, 3, false)


```
# Performance tests

To analyse the computational time of cumulants updates vs `Cumulants.jl` recalculation, we supply the executable script `comptimes.jl`. The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's maximum order, by default `m = 4`,
* `-n (vararg Int)`: numbers of marginal variables, by default `n = 40`,
* `-t (Int)`: number of realisations of random variable, by default `t = 500000`,
* `-u (vararg Int)`: number of realisations of update, by default `u = 10000, 15000, 20000`,
* `-b (Int)`: blocks size, by default `b = 4`,
* `-p (Int)`: numbers of processes, by default `p = 3`.

To analyse the computational time of cumulants updates for different block sizes `1 < b ≤ Int(√n)+2`, we supply the executable script `comptimesblocks.jl`.
The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (Int)`: numbers of marginal variables, by default `m = 48`,
* `-u (vararg Int)`: number of realisations of the update, by default `u = 10000, 20000`.
* `-p (Int)`: numbers of processes, by default `p = 3`.

To analyse the computational time of cumulants updates for different number of workers, we supply the executable script `comptimesprocs.jl`.
The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (Int)`: numbers of marginal variables, by default `m = 48`,
* `-u (vararg Int)`: number of realisations of the update, by default `u = 10000, 20000`,
* `-b (Int)`: blocks size, by default `b = 4`,
* `-p (Int)`: maximal numbers of processes, by default `p = 6`.

To plot computational times run executable `res/plotcomptimes.jl` on chosen `*.jld` file.


# Citing this work

Krzysztof Domino, Piotr Gawron, *On updates of high order cumulant tensors*, [arXiv:1701.06446](https://arxiv.org/abs/1701.06446)

This project was partially financed by the National Science Centre, Poland – project number 2014/15/B/ST6/05204.
