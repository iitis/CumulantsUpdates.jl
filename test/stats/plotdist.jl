#!/usr/bin/env julia

using Distributions
using Cumupdates
import Cumupdates: cormatgen, tcopulagmarg
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
srand(42)

t = 200000
cormat = cormatgen(10)
x = transpose(rand(MvNormal(cormat),t));
x1 = tcopulagmarg(cormat, t, 10);

y = x[:,1]
y5 = x[:,5]
y1_3 = vcat(x[1:div(2*t, 3),1], x1[(div(2*t, 3)+1):end,1])
y1_35 = vcat(x[1:div(2*t, 3),5], x1[(div(2*t, 3)+1):end,5])
y2_3 = vcat(x[1:div(t, 3),1], x1[(div(t, 3)+1):end,1])
y2_35 = vcat(x[1:div(t, 3),5], x1[(div(t, 3)+1):end,5])
y1 = x1[:,1]
y15 = x1[:,5]

function plotscatter(y, y1, file)
  fig, ax = subplots(figsize = (2.5, 2.))
  ax[:plot](y, y1, "o", color = "blue", markersize=1, alpha=0.1)
  fig[:savefig](file)
end


function plothist(y, file)
  fig, ax = subplots(figsize = (2.5, 2.))
  ax[:hist](y[:,1], 20)
  f = matplotlib[:ticker][:ScalarFormatter]()
  f[:set_powerlimits]((-1, 2))
  ax[:yaxis][:set_major_formatter](f)
  fig[:savefig](file)
end

function main()
  plotscatter(y, y5, "tup0.png")
  plotscatter(y1_3, y1_35, "tup13.png")
  plotscatter(y2_3, y2_35, "tup23.png")
  plotscatter(y1, y15, "tup1.png")
  plothist(y, "h0.eps")
  plothist(y1_3, "h13.eps")
  plothist(y2_3, "h23.eps")
end

main()
