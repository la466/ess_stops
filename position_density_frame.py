import collections
import bioUtilities.files as files
import bioUtilities.seq as seq
import bioUtilities.sim as sim
from bioUtilities.useful_sets.codons import stop_codons
from pybiomart import Dataset
import numpy as np



# def chunk_motif_hits


def position_density_frame(sequences, motif_set):

    counts = collections.defaultdict(lambda: collections.Counter())
    for id in sequences:
        if id:
            query_seq = sequences[id][0]
            start_frame = sequences[id][1]

            hits = seq.motif_hit_positions(query_seq, motif_set, flatten = False)
            for hit in hits:
                hit_start = hit[0]
                hit_frame = (hit_start + start_frame) % 3
                for pos in hit:
                    counts[hit_frame][pos] += 1

    frame_densities = collections.defaultdict(lambda: collections.defaultdict())

    for frame in counts:
        for pos in counts[frame]:
            frame_densities[frame][pos] = np.divide(counts[frame][pos], len(sequences))

    return frame_densities


length = 50


coding_exons_file = "source_data/human/coding_exons.fa"
query_file = "source_data/constitutive_decoys.fa"
introns_file = "source_data/human/introns.fa"
transcripts_file = "source_data/human/cds_transcripts.fa"

coding_exons = files.read_fasta(coding_exons_file)
coding_exons = {name.split("(")[0]: coding_exons.sequences[i] for i, name in enumerate(coding_exons.ids) if len(coding_exons.sequences[i]) >= 50}

query_ids = [i.split("|")[0] for i in files.read_fasta(query_file).ids]
coding_exons = {i: coding_exons[i] for i in coding_exons if i in query_ids}

introns = files.read_fasta(introns_file)
introns = {name.split("(")[0].split("-")[0]: introns.sequences[i] for i, name in enumerate(introns.ids) if len(introns.sequences[i]) >= 50}

transcripts = files.read_fasta(transcripts_file)
transcripts = {name: transcripts.sequences[i] for i, name in enumerate(transcripts.ids)}

exon_list = []
intron_list = []



def get_sequences(tn = None):

    exon_sequences = {}
    intron_sequences = {}

    for id in coding_exons:
        if id in introns:
            transcript_id = id.split(".")[0]
            exon = coding_exons[id]
            intron = introns[id]
            cds = transcripts[transcript_id]

            if tn and not len(intron) % 3 == 0:
                pass
            else:

                exon_start_index = cds.index(exon)
                exon_query_start = exon_start_index + len(exon) - length
                exon_query_frame = exon_query_start % 3
                exon_sequences[id] = [exon[-length:], exon_query_frame]


                intron_start_index = exon_start_index + len(exon)
                intron_start_frame = intron_start_index % 3
                intron_sequences[id] = [intron[:length], intron_start_frame]

    return exon_sequences, intron_sequences


exon_sequences, intron_sequences = get_sequences(tn = True)

exon_densities = position_density_frame(exon_sequences, stop_codons)
intron_densities = position_density_frame(intron_sequences, stop_codons)

output_file = "results/temp/constitutive_decoy_3n_frame_densities.csv"
with open(output_file, "w") as outfile:
    outfile.write("position,frame,density\n")
    for frame in [0, 1, 2]:
        for position in range(length):
            outfile.write("{0},{1},{2}\n".format(-length + position, frame, exon_densities[frame][position] if frame in exon_densities and position in exon_densities[frame] else 0))
        for position in range(length):
            outfile.write("{0},{1},{2}\n".format(1 + position, frame, intron_densities[frame][position] if frame in intron_densities and position in intron_densities[frame] else 0))
