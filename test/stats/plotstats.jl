#!/usr/bin/env julia

using JLD
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
using ArgParse



function plotstats(st, y, s, ylab, fname, legloc)
  mpl.rc("font", family="serif", size = 7)
  fig, ax = subplots(figsize = (2.5, 2.))
  col = ["red", "blue", "green", "black", "orange", "cyan"]
  marker = [":s", ":o", ":v", ":<", ":>", ":d"]
  for i in 1:size(st, 2)
    ax[:plot](y, st[:,i], marker[i], label = s[i], color = col[i], markersize=2.5, linewidth = 1)
  end
  PyPlot.ylabel(ylab, labelpad = -2.3)
  PyPlot.xlabel("\$ t_{up} / t \$", labelpad = -2)
  ax[:legend](fontsize = 6, loc = legloc, ncol = 1)
  subplots_adjust(bottom = 0.12,top=0.92)
  f = matplotlib[:ticker][:ScalarFormatter]()
  f[:set_powerlimits]((-2, 3))
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
  str = string(d["mu"])*string(d["n"])
  plotstats(d["st"][:,1:2], d["y"], ["m = 1" "m = 2"], "\$ ||C_{m}|| \$", "c1c2"*str*".eps", 5)
  plotstats(d["st"][:,3:4], d["y"], ["m = 3" "m = 4"], "\$ ||C_{m}|| \$", "c3c4"*str*".eps", 2)
  plotstats(hcat(d["skmax"], d["skmin"], d["kumax"], d["kumin"]), d["y"], ["max assymetry", "min assymetry", "max kurtosis", "min kurtosis"], " ", "1dstats"*str*".eps", 1)
end

main(ARGS)

function moin(args)
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
  st, sk, k, stats, y = getstats(t, n, tup, mu, m)
  str = "stats/"*string(mu)*string(n)
  plotstats(st, y, ["m = 1" "m = 2" "m = 3" "m = 4"], "\$ ||C_{m}|| \$", str*"normcums.eps", 2)
  plotstats(hcat(sk, k), y, ["1d assymetry", "1d curtosis"], " ", str*"1dstats.eps", 1)
  plotstats(stats, y, ["mean", "var", "assymetry", "curtosis"], " ", str*"chstats.eps", 6)
end
