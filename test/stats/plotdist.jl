#!/usr/bin/env julia

using Distributions
using Cumupdates
import Cumupdates: gendata, cormatgen
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
srand(42)

t = 2000000
cormat = cormatgen(10)
x = transpose(rand(MvNormal(cormat),t));
x1 = gendata(cormat, t, 10);

y = x[:,1]
y1 = x[:,5]
yy = vcat(x[1:div(t, 2),1], x1[(div(t, 2)+1):end,1])
yy1 = vcat(x[1:div(t, 2),10], x1[(div(t, 2)+1):end,5])
yyy = x1[:,1]
yyy1 = x1[:,5]

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
  plotscatter(y, y1, "pics/tup0.png")
  plotscatter(yy, yy1, "pics/tup12.png")
  plotscatter(yyy, yyy1, "pics/tup1.png")
  plothist(y, "pics/h0.eps")
  plothist(yy, "pics/h12.eps")
  plothist(yyy, "pics/h1.eps")
end

main()
