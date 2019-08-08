import collections
import bioUtilities.files as files
import bioUtilities.seq as seq
import bioUtilities.sim as sim
from bioUtilities.useful_sets.codons import stop_codons
from pybiomart import Dataset
import numpy as np


def position_density(sequences, motif_set):
    counts = collections.Counter()
    for sequence in sequences:
        hit_positions = seq.motif_hit_positions(sequence, motif_set)
        for position in hit_positions:
            counts[position] += 1
    densities = {i: np.divide(counts[i], len(sequences)) for i in sorted(counts)}
    return densities

length = 50


coding_exons_file = "source_data/human/coding_exons.fa"
introns_file = "source_data/human/introns.fa"

coding_exons = files.read_fasta(coding_exons_file)
coding_exons = {name.split("(")[0]: coding_exons.sequences[i] for i, name in enumerate(coding_exons.ids) if len(coding_exons.sequences[i]) >= 200}

introns = files.read_fasta(introns_file)
introns = {name.split("(")[0].split("-")[0]: introns.sequences[i] for i, name in enumerate(introns.ids) if len(introns.sequences[i]) >= 200}

exon_list = []
intron_list = []

for i, id in enumerate(coding_exons):
    # if i < 10:
    if i:
        if id in introns:
            exon_list.append(coding_exons[id][-length:])
            intron_list.append(introns[id][:length])


exon_density = seq.calc_motif_density(exon_list, stop_codons)
intron_density = seq.calc_motif_density(intron_list, stop_codons)

exon_position_densities = position_density(exon_list, stop_codons)
intron_position_densities = position_density(intron_list, stop_codons)


output_file = "results/temp/stop_codon_position_densities.csv"
with open(output_file, "w") as outfile:
    outfile.write("position,density\n")
    [outfile.write("{0},{1}\n".format(-length+i, exon_position_densities[i])) for i in exon_position_densities]
    [outfile.write("{0},{1}\n".format(1+i, intron_position_densities[i])) for i in intron_position_densities]
