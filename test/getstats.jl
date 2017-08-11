#!/usr/bin/env julia

using Cumupdates
using Cumulants
using Distributions
import Cumupdates: gendata, cormatgen, getstats
using JLD
using ArgParse



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
        help = "number of degree of freedom for student copula"
        default = 8
        arg_type = Int
    end
  parsed_args = parse_args(s)
  m = parsed_args["order"]
  n = parsed_args["nvar"]
  t = parsed_args["dats"]
  tup = parsed_args["updates"]
  mu = parsed_args["mu"]
  stsdict = Dict{String, Any}()
  st, sk, k, stats, y = getstats(t, n, tup, mu, m)
  str = "stats/stats"*string(mu)*string(n)*".jld"
  push!(stsdict, "st" => st)
  push!(stsdict, "sk" => sk)
  push!(stsdict, "k" => k)
  push!(stsdict, "stats" => stats)
  push!(stsdict, "y" => y)
  push!(stsdict, "mu" => mu)
  push!(stsdict, "n" => n)
  save(str, stsdict)
end

main(ARGS)
