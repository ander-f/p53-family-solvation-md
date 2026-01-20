using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, ComplexMixtures, PDBTools, CSV, DataFrames

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
    ylim=(-0.0005,0.008),
    minorticks=true)

base_dir="./../analyses/zscore"

#plot((layout=@layout [ a [ b ; c ] ]),
plot(layout = (4,1),
    leftmargin=0.8Plots.Measures.cm,
    rightmargin=0.4Plots.Measures.cm,
    bottommargin=0.4Plots.Measures.cm
    )

# Pallate of colors :tab10
c = ["#1f77b4", "#ff7f0e", "#2ca02c",
     ]

# data files
df_global = CSV.read("$base_dir/per_residue_global_hydration.csv", DataFrame)
df_p53_blocks = CSV.read("$base_dir/per_residue_p53_hydration.csv", DataFrame)
df_p63_blocks = CSV.read("$base_dir/per_residue_p63_hydration.csv", DataFrame)
df_p73_blocks = CSV.read("$base_dir/per_residue_p73_hydration.csv", DataFrame)

# Subplot 1: Global hydration contribution per residue
res_labels_p53 = df_global.p53_res .* string.(df_global.p53_resnum)
p53_hydration = df_global.p53_18A

res_labels_p63 = df_global.p63_res .* string.(df_global.p63_resnum.+39)
p63_hydration = df_global.p63_18A

res_labels_p73 = df_global.p73_res .* string.(df_global.p73_resnum)
p73_hydration = df_global.p73_18A

sp=1
plot!(p53_hydration,
    label="p53 (global)",
    color=c[1],
    xlabel="p53 Sequence",
    ylabel="Hydration Score (1.8 Å)",
    xticks=(1:5:length(res_labels_p53), string.(res_labels_p53[1:5:end])),
    xrotation=90,
    subplot=sp)
plot!(df_p53_blocks.mean_18A,
    ribbon=(df_p53_blocks.std_18A,),
    label="Mean of blocks ± SD",
    color=c[1],
    ls=:dash,
    alpha=0.5, 
    subplot=sp)

###############################
sp=2
plot!(p63_hydration,
    label="p63 (global)",
    color=c[2],
    xlabel="p63 Sequence",
    ylabel="Hydration Score (1.8 Å)",
    xticks=(1:5:length(res_labels_p63), string.(res_labels_p63[1:5:end])),
    xrotation=90,
    subplot=sp)
plot!(df_p63_blocks.mean_18A,
    ribbon=(df_p63_blocks.std_18A,),
    label="Mean of blocks ± SD",
    color=c[2],
    ls=:dash,
    alpha=0.5,
    subplot=sp)

###############################
sp=3
plot!(p73_hydration,
    label="p73 (global)",
    color=c[3],
    xlabel="p73 Sequence",
    ylabel="Hydration Score (1.8 Å)",
    xticks=(1:5:length(res_labels_p73), string.(res_labels_p73[1:5:end])),
    xrotation=90,
    subplot=sp)
plot!(df_p73_blocks.mean_18A,
    ribbon=(df_p73_blocks.std_18A,),
    label="Mean of blocks ± SD",
    color=c[3],
    ls=:dash,
    alpha=0.5,
    subplot=sp)

sp=4
plot!(p53_hydration,
    label="p53 (global)",
    color=c[1],
    xlabel="Sequence Position",
    ylabel="Hydration Score (1.8 Å)",
    xticks=(1:5:length(res_labels_p53)),
    xrotation=90,
    subplot=sp)

plot!(p63_hydration,
    label="p63 (global)",
    color=c[2],
    subplot=sp)

plot!(p73_hydration,
    label="p73 (global)",
    color=c[3],
    subplot=sp)

sp=1
annotate!(196,0.007,L"\mathrm{A}",14, subplot=sp)
sp=2
annotate!(196,0.007,L"\mathrm{B}",14, subplot=sp)
sp=3
annotate!(196,0.007,L"\mathrm{C}",14, subplot=sp)
sp=4
annotate!(196,0.007,L"\mathrm{D}",14, subplot=sp)

plot!(size=(840,1260))

#savefig("$base_dir/figures/fig_trapz_compare_18A.pdf")
