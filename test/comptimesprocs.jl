#!/usr/bin/env julia

using Cumulants
using SymmetricTensors
using CumulantsUpdates
using JLD
using ArgParse


function comptime(c, data::Matrix{Float64}, datup::Matrix{Float64})
  t = time_ns()
  updat(c, data, datup)
  Float64(time_ns()-t)/1.0e9
end

precomp(m::Int, data::Matrix{Float64}) =
  updat(momentarray(data[1:10, 1:10], m), data[1:10, 1:10], data[1:5, 1:10])

function savect(u::Vector{Int}, p::Int, n::Int, m::Int, b::Int)
  comptimes = zeros(p, length(u))
  data = randn(maximum(u)+10, n)
  precomp(m, data)
  M = momentarray(data, m, b)
  for i in 1:p
    rmprocs(procs()[2:end])
    addprocs(i)
    eval(Expr(:toplevel, :(@everywhere using CumulantsUpdates)))
    println("number of workers = ", nworkers())
    for k in 1:length(u)
      datup = randn(u[k], n)
      comptimes[i, k] = comptime(M, data, datup)
      println("u = ", u[k])
    end
  end
  filename = replace("res/$(m)_$(u)_$(n)_$(b)_nprocs.jld", "[", "")
  filename = replace(filename, "]", "")
  filename = replace(filename, " ", "")
  compt = Dict{String, Any}("cumulants"=> comptimes)
  push!(compt, "t" => u)
  push!(compt, "n" => n)
  push!(compt, "m" => m)
  push!(compt, "x" => "number of procs")
  push!(compt, "number of procs" => [collect(1:p)...])
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
        default = [10000, 20000]
        arg_type = Int
      "--blocksize", "-b"
        help = "the size of blocks of the block structure"
        default = 4
        arg_type = Int
      "--nprocs", "-p"
        help = "number of processes"
        default = 6
        arg_type = Int
    end
  parsed_args = parse_args(s)
  m = parsed_args["order"]
  n = parsed_args["nvar"]
  u = parsed_args["tup"]
  p = parsed_args["nprocs"]
  b = parsed_args["blocksize"]
  savect(u, p, n, m, b)
end

main(ARGS)
