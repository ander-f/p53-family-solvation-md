import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles
using LaTeXStrings
using Plots
using EasyFit
using ColorSchemes
using ComplexMixtures
using PDBTools
using Formatting

PROTEIN1 = "p53"
PROTEIN2 = "p63"
PROTEIN3 = "p73"

plot(layout = (2,1),
    leftmargin=3Plots.Measures.cm,
    bottommargin=3Plots.Measures.cm)

# Will use moving averages for more pretty graphs
ma(data) = movingaverage(data, 10).x

s = 10
# Plot defaults
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    xtickfont=font(s, plot_font),
    ytickfont=font(s, plot_font),
    legendfont=font(s, plot_font),
    guidefont=font(s, plot_font),
    linewidth=0.1,
    framestyle=:box,
    label=nothing,
    grid=false,
    dpi=300,
    xguidefontsize=s,
    xtickfontsize=s,
    yguidefontsize=s,
    ytickfontsize=s,
    minorticks=true)

###### PDB file of the system simulated - p53 WT
pdb = readPDB("./../coordinates/$(PROTEIN1)/$(PROTEIN1)_initial.pdb")
# Inform which is the solute
protein = select(pdb, "protein")
solute = AtomSelection(protein, nmols=1)
# Inform which is the solvent
water = select(pdb, "resname SOL")
solvent = AtomSelection(water, natomspermol=3)
# Obtain pretty labels for the residues in the x-axis
residues = collect(eachresidue(protein))
labels_WT = PDBTools.oneletter.(resname.(residues)) .* format.(resnum.(residues))
# Save the labels for later
open("./../labels_$(PROTEIN1).txt", "w") do io
    for label in labels_WT
        println(io, label)
    end
end

# Load results of a ComplexMixtures run
R = load("./../analyses/hydration/$(PROTEIN1)/results-water.json")
# columns equal to the number of residues
rescontrib = zeros(length(R.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib[:, ires] .= contributions(R, SoluteGroup(residue))
end
# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, R.d)
idmax = findfirst(d -> d > 4.2, R.d)
#sequence = getseq(protein)
rescontrib_p53 = rescontrib[:,1:198]

###### PDB file of the system simulated - p63
pdb = readPDB("./../coordinates/$(PROTEIN2)/$(PROTEIN2)_initial.pdb")
# Inform which is the solute
protein = select(pdb, "protein")
solute = AtomSelection(protein, nmols=1)
# Inform which is the solvent
water = select(pdb, "resname SOL")
solvent = AtomSelection(water, natomspermol=3)
# Obtain pretty labels for the residues in the x-axis
residues = collect(eachresidue(protein))
labels_p63 = PDBTools.oneletter.(resname.(residues)) .* format.(resnum.(residues))
# Save the labels for later
open("./../labels_$(PROTEIN2).txt", "w") do io
    for label in labels_p63
        println(io, label)
    end
end

# Load results of a ComplexMixtures run
R = load("./../analyses/hydration/$(PROTEIN2)/results-water.json")
# columns equal to the number of residues
rescontrib = zeros(length(R.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib[:, ires] .= contributions(R, SoluteGroup(residue))
end
# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, R.d)
idmax = findfirst(d -> d > 4.2, R.d)
#sequence = getseq(protein)
rescontrib_p63 = rescontrib[:,[1:91;94:200]]

###### PDB file of the system simulated - p73
pdb = readPDB("./../coordinates/$(PROTEIN3)/$(PROTEIN3)_initial.pdb")
# Inform which is the solute
protein = select(pdb, "protein")
solute = AtomSelection(protein, nmols=1)
# Inform which is the solvent
water = select(pdb, "resname SOL")
solvent = AtomSelection(water, natomspermol=3)
# Obtain pretty labels for the residues in the x-axis
residues = collect(eachresidue(protein))
labels_p73 = PDBTools.oneletter.(resname.(residues)) .* format.(resnum.(residues))
# Save the labels for later
open("./../labels_$(PROTEIN3).txt", "w") do io
    for label in labels_p73
        println(io, label)
    end
end

# Load results of a ComplexMixtures run
R = load("./../analyses/hydration/$(PROTEIN3)/results-water.json")
# columns equal to the number of residues
rescontrib = zeros(length(R.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib[:, ires] .= contributions(R, SoluteGroup(residue))
end
# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, R.d)
idmax = findfirst(d -> d > 4.2, R.d)
#sequence = getseq(protein)
rescontrib_p73 = rescontrib[:,[1:91;94:200]]

sp=1
# p53 - p63
rescontrib = rescontrib_p53 - rescontrib_p63
clims = [  minimum(rescontrib), maximum(rescontrib) ]
println(clims)
# We will plot only the range ARGS[2], for clarity
irange_arg = "1:198"
start, stop = parse.(Int, split(irange_arg, ":"))
irange = start:stop
# Use irange as UnitRange{Int64}
println(irange)
contourf!(irange, R.d[idmin:idmax], rescontrib[idmin:idmax, irange],
#    clim=(round.(-maximum(rescontrib); digits=2), (round.(maximum(rescontrib); digits=2))),
    clims=(-0.0301, 0.0301),
    color=cgrad(:RdBu, rev=true), # high values in red, small values in blue
    linewidth=0.1, linecolor=:black,
    colorbar=:true, levels=10,
    xlabel=L"\mathrm{Residue}", ylabel=L"\mathrm{r/\AA}",
    xticks=(1:5:198, labels_WT[1:5:198] ), xrotation=90,
    #xticks=(irange, labels_WT[1:198]), xrotation=90,
    xtickfont=font(8, plot_font),
    #size=(700, 400),
    subplot=sp
)

sp=2
# p53 - p73
rescontrib = rescontrib_p53 - rescontrib_p73
clims = [  minimum(rescontrib), maximum(rescontrib) ]
println(clims)
# We will plot only the range ARGS[2], for clarity
irange_arg = "1:198"
start, stop = parse.(Int, split(irange_arg, ":"))
irange = start:stop
# Use irange as UnitRange{Int64}
println(irange)
contourf!(irange, R.d[idmin:idmax], rescontrib[idmin:idmax, irange],
    #clims=(-0.01871,0.01871),
    clims=(-0.0301,0.0301), # considering sp=1
    color=cgrad(:RdBu, rev=true), # high values in red, small values in blue
    linewidth=0.1, linecolor=:black,
    colorbar=:true, levels=10,
    xlabel=L"\mathrm{Residue}", ylabel=L"\mathrm{r/\AA}",
    xticks=(1:5:198, labels_WT[1:5:198] ), xrotation=90,
    #xticks=([17, 21, 55, 74, 87, 91, 93, 108, 111], labels_WT[[17, 21, 55, 74, 87, 91, 93, 108, 111]] ), xrotation=90,
    xtickfont=font(8, plot_font),
    #size=(700, 400),
    subplot=sp
)

plot!(size=(1200*0.9, 700*0.9))

#savefig("./../figures/fig_cm_density_p63_p73_complete_DIFFERENCE.svg")
