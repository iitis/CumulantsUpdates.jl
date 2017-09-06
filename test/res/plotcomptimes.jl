#!/usr/bin/env julia

using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
using JLD
using ArgParse

o = 2
function singleplot(filename::String, name::String, compare::String = "")
    d = load(filename*".jld")
  if compare == ""
    comptimes = d[name]
    ylab = "computional time [s]"
  else
    comptimes = (d[name]./d[compare])
    ylab = "speedup"
  end
  x = d["x"]
  t = d["t"]
  n = d["n"]
  m = d["m"]
  mpl.rc("font", family="serif", size = 7)
  fig, ax = subplots(figsize = (2.5, 2.))
  col = ["red", "blue", "green", "yellow", "orange", "cyan"]
  marker = [":s", ":o", ":v", ":<", ":>", ":d"]
  for i in 1:size(comptimes, 2)
    if contains(filename, "nblocks")
      str = replace(@sprintf("%.0e", t[i]), "0", "")
      tt = "\$"*replace(str, "e+", "\\times 10^{")*"}\$"
      ax[:plot](d[x], comptimes[:,i], marker[i], label= "\$ t_{up} \$ = $tt", color = col[i], markersize=2.5, linewidth = 1)
    else
      tt =  n[i]
      ax[:plot](d[x][o:end], comptimes[o:end,i], marker[i], label= "n = $tt", color = col[i], markersize=2.5, linewidth = 1)
      i == size(comptimes, 2)? ax[:plot](d[x][o:end], d["tm"][o:end], ":<", label= "theoretical", color = "black", markersize=2.5, linewidth = 1):()
    end
  end
  subplots_adjust(bottom = 0.16, top=0.92, left = 0.12, right = 0.92)
  fx = matplotlib[:ticker][:ScalarFormatter]()
  fx[:set_powerlimits]((-1, 4))
  ax[:xaxis][:set_major_formatter](fx)
  if contains(filename, "nblocks")
    subplots_adjust(left = 0.15)
    PyPlot.ylabel(ylab, labelpad = 0.6)
    PyPlot.xlabel(x, labelpad = 0.6)
    ax[:legend](fontsize = 5, loc = 2, ncol = 1)
  else
    PyPlot.ylabel(ylab, labelpad = -1.0)
    PyPlot.xlabel("\$ t_{up} \$", labelpad = 0.6)
    ax[:legend](fontsize = 6, loc = 1, ncol = 1)
  end
  if maximum(comptimes) > 10
    f = matplotlib[:ticker][:ScalarFormatter]()
    f[:set_powerlimits]((-3, 2))
  else
    f = matplotlib[:ticker][:FormatStrFormatter]("%.1f")
  end
  ax[:yaxis][:set_major_formatter](f)
  fig[:savefig](name*filename*".eps")
end


"""
  pltspeedup(comptimes::Array{Float}, m::Int, n::Vector{Int}, T::Vector{Int}, label::String)

Returns a figure in .eps format of the computional speedup of cumulants function

"""

function pltspeedup(filename::String)
  d = load(filename)
  filename = replace(filename, ".jld", "")
  for f in d["functions"]
    singleplot(filename::String, f...)
  end
end


function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
    "file"
    help = "the file name"
    arg_type = String
  end
  parsed_args = parse_args(s)
  pltspeedup(parsed_args["file"])
end

main(ARGS)
