#!/usr/bin/env julia

using JLD
using ArgParse
using Cumupdates
using SymmetricTensors
import Cumupdates: cumulants, moment

"""
  compspeedups(fcalc::Function, fup::Function, m::Int, n::Int, t::Int, tup::Vector{Int})

Returns Vector, a computional speedup of m'th statics update of multivariate data
"""

function compspeedups(m::Int, n::Vector{Int}, t::Int, tup::Vector{Int}, f::Function,
                                                                        fup::Function,
                                                                        b::Int)
  compt = zeros(length(tup), length(n))
  updt = copy(compt)
  X = randn(15, 10)
  M = f(X[1:10,:], m, b)
  fup(M, X[1:10,:], X[10:15,:])
  for i in 1:length(n)
    println("n = ", n[i])
    println("compute")
    X = randn(t, n[i])
    t1 = Float64(time_ns())
    M = f(X, m, b)
    calctime = Float64(time_ns())-t1
    for j in 1:length(tup)
      println("update ", fup)
      println("tup = ", tup[j])
      Xup = rand(tup[j], n[i])
      t2 = Float64(time_ns())
      cup = fup(M, X, Xup)
      updt[j,i] = Float64(time_ns()) - t2
      compt[j,i] = calctime
    end
  end
  compt, updt
end

function savecomptime(m::Int, n::Vector{Int}, t::Int, tup::Vector{Int}, b::Int, p::Int)
  filename = replace("res/$(m)_$(t)_$(n)_$(tup)_$(p).jld", "[", "")
  filename = replace(filename, "]", "")
  filename = replace(filename, " ", "")
  compt = Dict{String, Any}()
  momtime = compspeedups(m, n, t, tup, moment, momentupdat, b)
  cumtime = compspeedups(m, n, t, tup, cumulants, cumulantsupdat, b)
  push!(compt, "moment" => momtime[1])
  push!(compt, "moment update" => momtime[2])
  push!(compt, "cumulants" => cumtime[1])
  push!(compt, "cumulants updat" => cumtime[2])
  push!(compt, "tm" => t./(2*tup))
  push!(compt, "t" => t)
  push!(compt, "n" => n)
  push!(compt, "m" => m)
  push!(compt, "tup" => tup)
  push!(compt, "x" => "tup")
  push!(compt, "functions" => [["moment", "moment update"], ["cumulants", "cumulants updat"]])
  save(filename, compt)
end

"""
  main(args)

Returns plots of the speedup of updates. Takes optional arguments from bash
"""
function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
      "--order", "-m"
        help = "m, the order of cumulant, ndims of cumulant's tensor"
        default = 4
        arg_type = Int
        "--blocksize", "-b"
        help = "the size of blocks of the block structure"
        default = 3
        arg_type = Int
      "--nvar", "-n"
        nargs = '*'
        default = [30]
        help = "n, numbers of marginal variables"
        arg_type = Int
      "--dats", "-t"
        help = "t, numbers of data records"
        default = 300000
        arg_type = Int
      "--updates", "-u"
        help = "u, size of the update"
        nargs = '*'
        default = [30000, 40000, 50000]
        arg_type = Int
      "--nprocs", "-p"
        help = "number of processes"
        default = 1
        arg_type = Int
    end
  parsed_args = parse_args(s)
  m = parsed_args["order"]
  n = parsed_args["nvar"]
  t = parsed_args["dats"]
  tup = parsed_args["updates"]
  b = parsed_args["blocksize"]
  p = parsed_args["nprocs"]
  if p > 1
    addprocs(p)
    eval(Expr(:toplevel, :(@everywhere using Cumupdates)))
  end
  println("number of workers = ", nworkers())
  savecomptime(m, n, t, tup, b, p)
end

main(ARGS)
