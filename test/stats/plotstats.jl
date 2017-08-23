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


function pltslvnorm(str::String)
  d = load(str)
  d1 = load(replace(str, "25", "30"))
  d2 = load(replace(str, "25", "35"))
  d3 = load(replace(str, "25", "40"))
  mpl.rc("font", family="serif", size = 5)
  fig, ax = subplots(figsize = (2.5, 2.))
  ax[:plot](d["nc2"], d["mVinv(sl)"], "s", color = "r", label = "n = $(d["n"])", markersize=2.5)
  ax[:plot](d1["nc2"], d1["mVinv(sl)"], "o", color = "b", label = "n = $(d1["n"])", markersize=2.5)
  ax[:plot](d2["nc2"], d2["mVinv(sl)"], "<", color = "g", label = "n = $(d2["n"])", markersize=2.5)
  ax[:plot](d3["nc2"], d3["mVinv(sl)"], ">", color = "gray", label = "n = $(d3["n"])", markersize=2.5)
  ax[:plot](d["nc2"], [fill(0.25, length(d["nc2"]))...], ":", color = "black", label = "proposed model", linewidth = 0.5)
  ax[:legend](fontsize = 5, loc = 2, ncol = 1)
  PyPlot.xlabel("\$||\\Sigma||\$", labelpad = 0.2)
  f = matplotlib[:ticker][:FormatStrFormatter]("%.2f")
  ax[:yaxis][:set_major_formatter](f)
  PyPlot.ylabel("slope", labelpad = -2.0)
  subplots_adjust(bottom = 0.12,top=0.92)
  str = replace(str, ".jld", ".eps")
  fig[:savefig](replace(str, "25", ""))
end


function pltinvslvmu(str::String)
  d = load(str)
  mpl.rc("font", family="serif", size = 5)
  fig, ax = subplots(figsize = (2.5, 2.))
  col = ["red", "blue", "green", "black", "orange", "cyan"]
  marker = ["s", "o", "v", "<", ">", "d"]
  for i in 1: length(d["nc2"])
    a,b = linreg(d["sl"][i,:], d["mu"])
    ax[:plot](d["sl"][i,:], d["mu"], marker[i], label = "\$||\\Sigma||\$ =$(round(d["nc2"][i], 1))", color = col[i], markersize=2.5)
    ax[:plot](d["sl"][i,:], [a+b*j for j in d["sl"][i,:]], color = col[i], ":", linewidth = 0.25)
  end
  ax[:legend](fontsize = 5, loc = 4, ncol = 1)
  PyPlot.xlabel("1/a", labelpad = -1.8)
  PyPlot.ylabel("\$\\mu\$", labelpad = -1.0)
  f = matplotlib[:ticker][:FormatStrFormatter]("%.0f")
  ax[:yaxis][:set_major_formatter](f)
  subplots_adjust(bottom = 0.12,top=0.92)
  str = replace(str, ".jld", ".eps")
  fig[:savefig](str)
end

function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
    "file"
    help = "the file name"
    arg_type = String
  end
  parsed_args = parse_args(s)
  file = parsed_args["file"]
  if contains(file, "nvmu")
    if contains(file, "nvmu25")
      pltslvnorm(file)
    end
    pltinvslvmu(file)
  else
    d = load(file)
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
    plotstats(c34, d["y"], leg34, "\$ h_{m}(\\mathbf{X})\$", "c3c4"*str*".eps", 2, -1.3)
    plotstats(hcat(d["skmax"], d["skmin"], d["kumax"], d["kumin"]), d["y"], ["max assym.", "min assym.", "max kurt.", "min kurt."], " ", "1d"*str*".eps", 5)
  end
end

main(ARGS)
