#!/usr/bin/env julia

using Cumulants
using SymmetricTensors
using Cumupdates
using JLD
using ArgParse


function comptime(c, data::Matrix{Float64}, datup::Matrix{Float64})
  t = time_ns()
  cumulantsupdat(c, data, datup)
  Float64(time_ns()-t)/1.0e9
end


function savect(tup::Vector{Int}, n::Int, m::Int)
  maxb = round(Int, sqrt(n))
  comptimes = zeros(maxb-1, length(tup))
  println("max block size = ", maxb)
  data = randn(maximum(tup)+10, n)
  cumulantsupdat(cumulants(data[1:10, 1:10], 4), data[1:10, 1:10], data[1:5, 1:10])
  for b in 2:maxb
    c = cumulants(data, m, b)
    for k in 1:length(tup)
      datup = randn(tup[k], n)
      comptimes[b-1, k] = comptime(c, data, datup)
      println("tup = ", tup[k])
      println("bloks size = ", b)
    end
  end
  filename = replace("res/$(m)_$(tup)_$(n)_nblocks.jld", "[", "")
  filename = replace(filename, "]", "")
  compt = Dict{String, Any}("cumulants"=> comptimes)
  push!(compt, "t" => tup)
  push!(compt, "n" => n)
  push!(compt, "m" => m)
  push!(compt, "x" => "block size")
  push!(compt, "block size" => [collect(2:maxb)...])
  push!(compt, "functions" => [["cumulants"]])
  save(filename, compt)
end


function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
      "--order", "-m"
        help = "m, the order of cumulant, ndims of cumulant's tensor"
        default = 4
        arg_type = Int
      "--nvar", "-n"
        default = 48
        help = "n, numbers of marginal variables"
        arg_type = Int
      "--tup", "-u"
        help = "u, numbers of data updates"
        nargs = '*'
        default = [20000, 40000]
        arg_type = Int
    end
  parsed_args = parse_args(s)
  m = parsed_args["order"]
  n = parsed_args["nvar"]
  tup = parsed_args["tup"]
  savect(tup::Vector{Int}, n::Int, m::Int)
end

main(ARGS)
