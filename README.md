[![Build Status](https://travis-ci.org/ZKSI/Cumupdates.jl.svg?branch=master)](https://travis-ci.org/ZKSI/Cumupdates.jl)


# Cumupdates.jl

Updates following statistics of n-variate data
* `m`'th moment tensor,
* a sequence of cumulant tensors of order `1,2,...,m`,
for `m >= 1`.

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
julia> Pkg.clone("https://github.com/ZKSI/Cumupdates.jl")
```

to install the files.  Julia 0.6 is required.

## Parallel computation

For parallel computation just run
```julia
julia> addprocs(n)
julia> @everywhere using Cumupdates
```

## Functions

### Moment update

To update the moment tensor `SymmetricTensor{T, m}` for `m <= 1` given original data `X \in R^{t, n}` and the update `Xup \in R^{tup, n}` and `t > tup` just run

```julia
julia> momentupdat(M::SymmetricTensor{T, m}, X::Matrix{T}, Xup::Matrix{T}) where {T<:AbstractFloat, m}
```

Returns a `SymmetricTensor{T, m}` of the moment tensor of updated multivariate data: `Xprim = vcat(X,Xup)[1+tup:end, :] \in R^{t, n}`. The output of `momentupdat(M, X, Xup)` corresponds to the output of `Cumulants.moment(Xprim, m)`, where `typeof(M) = SymmetricTensor{T, m}`. However if `tup << t` `momentupdat()` is much faster.

```julia
julia> x = ones(6, 2);

julia> m = moment(x, 3)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[1.0 1.0; 1.0 1.0]

[1.0 1.0; 1.0 1.0]],2,1,2,true)

julia> y = 2*ones(2,2);

julia> momentupdat(m, x, y)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

julia> moment(vcat(x,y)[3:end,:],3)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

```

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
* `-n (vararg Int)`: numbers of marginal variables, by default `m = 30`,
* `-t (Int)`: number of realisations of random variable, by default `t = 2500000`,
* `-u (vararg Int)`: number of realisations of update, by default `u = 25000, 30000, 35000, 40000`,
* `-b (Int)`: blocks size, by default `b = 3`,
* `-p (Int)`: numbers of processes, by default `p = 1`.

To analyse the computational time of cumulants updates for different block sizes `1 < b =< Int(sqrt(n))`, we supply the executable script `comptimesblocks.jl`.
The script saves computational times to the `res/*.jld` file. The scripts accept following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (Int)`: numbers of marginal variables, by default `m = 48`,
* `-u (vararg Int)`: number of realisations of the update, by default `u = 20000, 40000`.
* `-p (Int)`: numbers of processes, by default `p = 1`.

To plot computional times run executable `res/plotcomptimes.jl` on chosen `*.jld` file.

# Statistical analysis on generated data

Script `getstats.jl` returns `stats/*.jld` file of statistics computed for data generated from Gaussian copula and updated by data generated from t-Student copula, given following parameters:
* `-m (Int)`: maximum cumulant's order for statistical analysis, by default `m = 4`,
* `-n (Vararg Int)`: numbers of marginal variables, by default `m = [20, 24]`,
* `-t (Int)`: number of data records, by default `t = 500000`,
* `-u (Int)`: size of the update, by default `u = 50000`,
* `-d (Int)`: number of degree of freedom for t-Student copula, by default `d = 10`.

To plot statistics run executable `stats/plotstats.jl` on chosen `*.jld` file.

# Citing this work


Krzysztof Domino, Piotr Gawron, *On updates of high order cumulant tensors*, [arXiv:1701.06446](https://arxiv.org/abs/1701.06446)

