using MolSimToolkit
using Plots
using DelimitedFiles

base_dir="./../analyses/convergence"

### p53 system ###
p53_residues = [110, 114]  # R110, L114
for i in 1:length(p53_residues)

    data_rep1 = "$base_dir/coord_number_residue_p53_rep1_$(p53_residues[i]).dat"
    coord_number_rep1 = vec(readdlm(data_rep1))
    data_rep2 = "$base_dir/coord_number_residue_p53_rep2_$(p53_residues[i]).dat"
    coord_number_rep2 = vec(readdlm(data_rep2))

    mean_coord_number = (coord_number_rep1 .+ coord_number_rep2) ./ 2
    b = block_average(mean_coord_number)

    plot(b)

#    savefig("$base_dir/figures/coord_number_convergence_residue_p53_$(p53_residues[i]).svg")

end


### p63 system ###
p63_residues = [139, 143]  # D139, Q143
for i in 1:length(p63_residues)

    data_rep1 = "$base_dir/coord_number_residue_p63_rep1_$(p63_residues[i]).dat"
    coord_number_rep1 = vec(readdlm(data_rep1))
    data_rep2 = "$base_dir/coord_number_residue_p63_rep2_$(p63_residues[i]).dat"
    coord_number_rep2 = vec(readdlm(data_rep2))

    mean_coord_number = (coord_number_rep1 .+ coord_number_rep2) ./ 2
    b = block_average(mean_coord_number)

    plot(b)

#    savefig("$base_dir/figures/coord_number_convergence_residue_p63_$(p63_residues[i]).svg")

end

### p73 system ###
p73_residues = [128, 132]  # E128, Q132
for i in 1:length(p73_residues)

    data_rep1 = "$base_dir/coord_number_residue_p73_rep1_$(p73_residues[i]).dat"
    coord_number_rep1 = vec(readdlm(data_rep1))
    data_rep2 = "$base_dir/coord_number_residue_p73_rep2_$(p73_residues[i]).dat"
    coord_number_rep2 = vec(readdlm(data_rep2))

    mean_coord_number = (coord_number_rep1 .+ coord_number_rep2) ./ 2
    b = block_average(mean_coord_number)

    plot(b)

#    savefig("$base_dir/figures/coord_number_convergence_residue_p73_$(p73_residues[i]).svg")

end
