import Pkg
Pkg.activate(".")
Pkg.instantiate()
using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, ComplexMixtures, PDBTools

# Will use moving averages for more pretty graphs
ma(data) = movingaverage(data,10).x

# Plot defaults
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1.5,
    framestyle=:box,
    label=nothing,
    grid=false,
    dpi=300,
    legend=:topleft,
    minorticks=false)

plot()

proteins = [ "p53",
             "p63",
             "p73"
            ]

# Pallate of colors :tab10
c = ["#1f77b4", "#9467bd", "#8c564b"]

for i in 1:length(proteins)

    R = ComplexMixtures.load("./../analyses/hydration/$(proteins[i])/results-water.json") 
    plot!(xaxis = ("r/Å",0:2.5:7.5),ylabel=L"\mathrm{MDDF}",
          )
    plot!(R.d,(ma(R.mddf)),
          color=c[i],
          xlim=[1,4.0],
          ylim=[0.8,1.9],
          label = proteins[i],
          legend=:topright,
          )

  end

plot!(size=(600,400))

#savefig("./../figures/fig_MDDF.svg")
