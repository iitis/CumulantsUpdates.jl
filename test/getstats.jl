#!/usr/bin/env julia

using Cumupdates
using Distributions
import Cumupdates: cormatgen, cumulants, tcopulagmarg, gcopulatmarg, tdistdat, normdist
using JLD
using ArgParse

srand(41)
"""
  getstats(n::Int, t::Int, wsize::Int, mu::Int, m::Int)

Returns statistics for randmly generfated data with gaussian marginals and different copulas
"""


function getstats(f::Function, fup::Function, t::Int, nv::Vector{Int}, u::Int, mu::Int, m::Int)
  k = div(t, u)
  tup = [u*i/t for i in 0:k]
  cormat = cormatgen(nv[1]);
  x = f(cormat, t, mu)
  cn = cumnorms(x, m, true, 1)
  skmax = maximum([skewness(x[:,p]) for p in 1:nv[1]])
  kumax = maximum([kurtosis(x[:,p]) for p in 1:nv[1]])
  skmin = minimum([skewness(x[:,p]) for p in 1:nv[1]])
  kumin = minimum([kurtosis(x[:,p]) for p in 1:nv[1]])
  println(nv[1])
  for i in 1:k
    xup = fup(cormat, u, mu)
    x = vcat(x, xup)[(size(xup, 1)+1):end,:]
    skmax = vcat(skmax, maximum([skewness(x[:,p]) for p in 1:nv[1]]))
    kumax = vcat(kumax, maximum([kurtosis(x[:,p]) for p in 1:nv[1]]))
    skmin = vcat(skmin, minimum([skewness(x[:,p]) for p in 1:nv[1]]))
    kumin = vcat(kumin, minimum([kurtosis(x[:,p]) for p in 1:nv[1]]))
    cn = hcat(cn, cumupdatnorms(xup, true, 1))
    println(tup[i])
  end
  for n in nv[2:end]
    cormat = cormatgen(n);
    x = f(cormat, t, mu)
    cn = hcat(cn, cumnorms(x, m, true, 1))
    println(n)
    for i in 1:k
      xup = fup(cormat, u, mu)
      cn = hcat(cn, cumupdatnorms(xup, true, 1))
      println(tup[i])
    end
  end
  cn', tup, skmax, skmin, kumax, kumin
end


function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
      "--order", "-m"
        help = "m, the order of cumulant, ndims of cumulant's tensor"
        default = 4
        arg_type = Int
      "--nvar", "-n"
        nargs = '*'
        default = [20, 24]
        help = "n, numbers of marginal variables"
        arg_type = Int
      "--dats", "-t"
        help = "t, numbers of data records"
        default = 200000
        arg_type = Int
      "--updates", "-u"
        help = "u, size of the update"
        default = 50000
        arg_type = Int
      "--mu", "-d"
        help = "number of degree of freedom for t-Student copula"
        default = 10
        arg_type = Int
    end
  parsed_args = parse_args(s)
  m = parsed_args["order"]
  n = parsed_args["nvar"]
  t = parsed_args["dats"]
  tup = parsed_args["updates"]
  mu = parsed_args["mu"]
  stsdict = Dict{String, Any}()
  f = [normdist, gcopulatmarg, normdist]
  fup = [tcopulagmarg, tdistdat, tdistdat]
  for i in 1:3
    st, y, skmax, skmin, kumax, kumin = getstats(f[i], fup[i], t, n, tup, mu, m)
    str = "stats/stats"*string(i)*replace(string(t)*"_"*string(mu)*string(n)*".jld", "[", "_")
    str = replace(str, "]", "")
    push!(stsdict, "st" => st)
    push!(stsdict, "skmax" => skmax)
    push!(stsdict, "skmin" => skmin)
    push!(stsdict, "kumax" => kumax)
    push!(stsdict, "kumin" => kumin)
    push!(stsdict, "y" => y)
    push!(stsdict, "mu" => mu)
    push!(stsdict, "n" => n)
    save(str, stsdict)
  end
end

main(ARGS)
