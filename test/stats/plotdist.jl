#!/usr/bin/env julia

using Distributions
using CumulantsUpdates
import CumulantsUpdates: cormatgen, tcopulagmarg
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
mpl.use("Agg")
using PyPlot
srand(42)

t = 10000
cormat = cormatgen(10)
x = transpose(rand(MvNormal(cormat),t));
x1 = tcopulagmarg(cormat, t, 4);

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
  mpl.rc("font", family="serif", size = 7)
  subplots_adjust(bottom = 0.14, top=0.92, left = 0.10, right = 0.92)
  PyPlot.ylabel("\$ X_5 \$", labelpad = -7.)
  PyPlot.xlabel("\$ X_1 \$", labelpad = 0.5)
  ax[:plot](y, y1, "o", color = "blue", markersize=1., alpha = 0.1)
  fig[:savefig](file)
end


function main()
  plotscatter(y, y5, "tup0.png")
  plotscatter(y, y5, "tup0.pdf")
  plotscatter(y1_3, y1_35, "tup13.pdf")
  plotscatter(y2_3, y2_35, "tup23.pdf")
  plotscatter(y1, y15, "tup1.pdf")
end

main()
