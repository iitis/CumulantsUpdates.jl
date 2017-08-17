#!/usr/bin/env julia

using JLD
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
using ArgParse



function plotstats(st, y, s, ylab, fname, legloc, of::Float64 = -2.0)
  mpl.rc("font", family="serif", size = 5)
  fig, ax = subplots(figsize = (2.5, 2.))
  col = ["red", "blue", "green", "black", "orange", "cyan"]
  marker = [":s", ":o", ":v", ":<", ":>", ":d"]
  for i in 1:size(st, 2)
    ax[:plot](y, st[:,i], marker[i], label = s[i], color = col[i], markersize=2.5, linewidth = 1)
  end
  PyPlot.ylabel(ylab, labelpad = of)
  PyPlot.xlabel("\$ t_{up} / t \$", labelpad = -2.0)
  ax[:legend](fontsize = 5, loc = legloc, ncol = 1)
  subplots_adjust(bottom = 0.12,top=0.92)
  f = matplotlib[:ticker][:ScalarFormatter]()
  f[:set_powerlimits]((-1, 2))
  ax[:yaxis][:set_major_formatter](f)
  fig[:savefig](fname)
end

function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
    "file"
    help = "the file name"
    arg_type = String
  end
  parsed_args = parse_args(s)
  d = load(parsed_args["file"])
  str = replace(parsed_args["file"], ".jld", "")
  l = length(d["y"])
  c12 = d["st"][1:l,1:2]
  c34 = d["st"][1:l,3:4]
  leg12 = ["m = 1 n = $(d["n"][1])" "m = 2 n = $(d["n"][1])"]
  leg34 = ["m = 3 n = $(d["n"][1])" "m = 4 n = $(d["n"][1])"]
  for j in 2:div(length(d["st"][:,1]), l)
    c12 = hcat(c12, d["st"][(l*(j-1)+1):j*l,1:2])
    c34 = hcat(c34, d["st"][(l*(j-1)+1):j*l,3:4])
    leg12 = hcat(leg12, ["m = 1 n = $(d["n"][j])" "m = 2 n = $(d["n"][j])"])
    leg34 = hcat(leg34, ["m = 3 n = $(d["n"][j])" "m = 4 n = $(d["n"][j])"])
  end
  plotstats(c12, d["y"], leg12, "\$ ||C_{m}|| \$", "c1c2"*str*".eps", 5, 1.)
  plotstats(c34, d["y"], leg34, "\$ ||C_{m}||/||C_{2}||^{m/2}\$", "c3c4"*str*".eps", 2)
  plotstats(hcat(d["skmax"], d["skmin"], d["kumax"], d["kumin"]), d["y"], ["max assym.", "min assym.", "max kurt.", "min kurt."], " ", "1d"*str*".eps", 5)
end

main(ARGS)
