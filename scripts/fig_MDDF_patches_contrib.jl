import Pkg
Pkg.activate(".")
Pkg.instantiate()
using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, ComplexMixtures, PDBTools

# Will use moving averages for more pretty graphs
ma(data) = movingaverage(data,10).x

# Plot defaults
#margin = 1.5Plots.mm
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

plot(layout = (2,3),
    leftmargin=0.2Plots.Measures.cm,
    rightmargin=0.4Plots.Measures.cm,
    bottommargin=0.4Plots.Measures.cm
    )

# Pallate of colors :tab10
c = ["#1f77b4", "#ff7f0e", "#2ca02c",
     ]

######## p53 ########
patches_p53 = [ 108 116
                128 131
                147 156
                164 167
                179 188
                198 210
              ]

for i in 1:length(patches_p53[:,1])

    # Calcuylate the group contributions
    atoms = readPDB("./../coordinates/p53/p53_initial.pdb")
    protein = select(atoms,"protein")
    water = select(atoms,"resname SOL")
    solute = AtomSelection(protein,nmols=1)
    solvent = AtomSelection(water,natomspermol=3)

    R = ComplexMixtures.load("./../analyses/hydration/p53/results-water.json") 
    start = patches_p53[i,1]
    stop = patches_p53[i,2]
    patches = contributions(R, SoluteGroup(PDBTools.select(atoms,"protein and resnum >= $start and resnum <= $stop")))
    plot!(R.d,(ma(patches)),
          xaxis = ("r/Å",0:2:8.0),ylabel=L"\mathrm{g_{ps} \ (r)}",
          color=c[1],
          xlim=[0,8.0],
          label = "p53",
          legend=:bottomright,
          title = "P$i - residues $(start)-$(stop)",
          titlefontsize = 10,
          sp=i
          )

end

######## p63 ########
patches_p63 = [ 137 145
                157 160
                176 185
                193 196
                208 219
                229 241
              ]

for i in 1:length(patches_p63[:,1])

    # Calcuylate the group contributions
    atoms = readPDB("./../coordinates/p63/p63_initial.pdb")
    protein = select(atoms,"protein")
    water = select(atoms,"resname SOL")
    solute = AtomSelection(protein,nmols=1)
    solvent = AtomSelection(water,natomspermol=3)

    R = ComplexMixtures.load("./../analyses/hydration/p63/results-water.json") 
    start = patches_p63[i,1]
    stop = patches_p63[i,2]
    patches = contributions(R, SoluteGroup(PDBTools.select(atoms,"protein and resnum >= $start and resnum <= $stop")))
    plot!(R.d,(ma(patches)),
          xaxis = ("r/Å",0:2:8.0),ylabel=L"\mathrm{g_{ps} \ (r)}",
          color=c[2],
          xlim=[0,8.0],
          label = "p63",
          legend=:bottomright,
          sp=i
          )

end

######## p73 ########
patches_p73 = [ 126 134
                146 149
                165 174
                182 185
                197 208 
                218 230
              ]

for i in 1:length(patches_p73[:,1])

    # Calcuylate the group contributions for the A system in urea
    atoms = readPDB("./../coordinates/p73/p73_initial.pdb")
    protein = select(atoms,"protein")
    water = select(atoms,"resname SOL")
    solute = AtomSelection(protein,nmols=1)
    solvent = AtomSelection(water,natomspermol=3)

    R = ComplexMixtures.load("./../analyses/hydration/p73/results-water.json") 
    start = patches_p73[i,1]
    stop = patches_p73[i,2]
    patches = contributions(R, SoluteGroup(PDBTools.select(atoms,"protein and resnum >= $start and resnum <= $stop")))
    plot!(R.d,(ma(patches)),
          xaxis = ("r/Å",0:2:8.0),ylabel=L"\mathrm{g_{ps} \ (r)}",
          color=c[3],
          xlim=[0,8.0],
          label = "p73",
          legend=:bottomright,
          sp=i
          )

end

plot!(size=(840*0.9,540*0.9))

#savefig("./../figures/fig_MDDF_patches_contrib_p53_p63_p73.svg")
