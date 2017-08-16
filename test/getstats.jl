#!/usr/bin/env julia

using Cumupdates
using Distributions
import Cumupdates: gendata, cormatgen, cumulants
using JLD
using ArgParse

srand(41)
"""
  getstats(n::Int, t::Int, wsize::Int, mu::Int, m::Int)

Returns statistics for randmly generfated data with gaussian marginals and different copulas
"""


function getstats(t::Int = 200000, n::Int = 20, wsize::Int = 10000, mu::Int = 10,
                                                                    m::Int = 4)
  cormat = cormatgen(n);
  x = transpose(rand(MvNormal(cormat),t));
  k = div(t, wsize)
  cn = cumnorms(x, m, true, 1)
  tup = [wsize*i/t for i in 0:k]
  skmax = maximum([skewness(x[:,p]) for p in 1:n])
  kumax = maximum([kurtosis(x[:,p]) for p in 1:n])
  skmin = minimum([skewness(x[:,p]) for p in 1:n])
  kumin = minimum([kurtosis(x[:,p]) for p in 1:n])
  for i in 1:k
    xup = gendata(cormat, wsize, mu)
    x = vcat(x, xup)[(size(xup, 1)+1):end,:]
    skmax = hcat(skmax, maximum([skewness(x[:,p]) for p in 1:n]))
    kumax = hcat(kumax, maximum([kurtosis(x[:,p]) for p in 1:n]))
    skmin = hcat(skmin, minimum([skewness(x[:,p]) for p in 1:n]))
    kumin = hcat(kumin, minimum([kurtosis(x[:,p]) for p in 1:n]))
    cn = hcat(cn, cumupdatnorms(xup, true, 1))
    println(tup[i])
  end
  cn', tup, skmax', skmin', kumax', kumin'
end


function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
      "--order", "-m"
        help = "m, the order of cumulant, ndims of cumulant's tensor"
        default = 4
        arg_type = Int
      "--nvar", "-n"
        default = 20
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
  st, y, skmax, skmin, kumax, kumin = getstats(t, n, tup, mu, m)
  str = "stats/stats"*string(mu)*string(n)*".jld"
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

main(ARGS)
