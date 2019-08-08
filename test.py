import collections
import bioUtilities.files as files
from maxentpy import maxent  # use normal version of maxent
from maxentpy.maxent import load_matrix5, load_matrix3
import re
import numpy as np


exons_file = "source_data/human/coding_exons.fasta"
introns_file = "source_data/human/introns.fasta"
transcripts_file = "source_data/human/cds_transcripts.fasta"

output_file = "splice_sites.fa"
decoy_file = "decoys.fa"

ess_file = "source_data/motif_sets/ess_fas_hex3.txt"


def get_sequences(length):

    exon_names, exon_seqs = files.read_fasta(exons_file)
    exons = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    [exons[name.split('(')[0].split('.')[0]][int(name.split('(')[0].split('.')[1])].append(exon_seqs[i]) for i, name in enumerate(exon_names) if len(exon_seqs[i]) > length]
    exons = {id: {exon_id: exons[id][exon_id][0] for exon_id in exons[id]} for id in exons}

    intron_names, intron_seqs = files.read_fasta(introns_file)
    introns = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    [introns[name.split('(')[0].split('.')[0]][int(name.split('(')[0].split('.')[1].split('-')[0])].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split('.')[0] in exons]
    introns = {id: {intron_id: introns[id][intron_id][0] for intron_id in introns[id]} for id in introns}




    #
    with open(output_file, 'w') as outfile:
        for id in exons:
            for exon_id in exons[id]:
                if id in introns:
                    if exon_id in introns[id]:
                        outfile.write(">{0}.{1}\n{2}{3}\n".format(id, exon_id, exons[id][exon_id][-length:].lower(), introns[id][exon_id][:length]))


    matrix5 = load_matrix5()
    entries = files.read_fasta(output_file)
    entries = {id: entries.sequences[i] for i, id in enumerate(entries.ids)}


    decoys = []

    for id in entries:
        seq = entries[id]

        splice_site = int(len(seq)/2)

        splice_site_seq = seq[splice_site-3:splice_site+6]
        real_splice_site_max_ent = maxent.score5(splice_site_seq, matrix=matrix5)


        kept = False
        for i in range(1, len(seq) - splice_site - 5):
            if not kept:
                query = seq[splice_site + i - 3:splice_site + i + 6]
                query = "{0}{1}".format(query[:3].lower(), query[3:])
                max_ent_score = maxent.score5(query, matrix=matrix5)
                if max_ent_score >= real_splice_site_max_ent:
                    print(id, real_splice_site_max_ent, i, query, max_ent_score)
                    decoys.append(id)
                    kept = True


    with open(decoy_file, "w") as outfile:
        [outfile.write(">{0}\n{1}\n".format(id, entries[id])) for id in decoys]


def calc_hits(seq, motif_set):
    motif_search = re.compile("(?=({0}))".format("|".join(motif_set)))
    hits = []
    matches = re.finditer(motif_search, seq)
    [hits.extend(list(range(hit.span()[0], hit.span()[0] + len(hit.group(1))))) for hit in matches]
    hits = sorted(list(set(hits)))
    return hits



def calc_densities(seq_list, motif_set):
    motif_search = re.compile("(?=({0}))".format("|".join(motif_set)))
    hit_count = collections.Counter()
    for id in seq_list:
        seq = seq_list[id]
        seq = seq.upper()
        hits = []
        matches = re.finditer(motif_search, seq)
        [hits.extend(list(range(hit.span()[0], hit.span()[0] + len(hit.group(1))))) for hit in matches]
        hits = sorted(list(set(hits)))
        for i in hits:
            hit_count[i-length] += 1
    densities = {i: np.divide(hit_count[i], len(seq_list)) for i in hit_count}
    return densities


length = 50
# get_sequences(length)

motifs = [i[0] for i in files.read_many_fields(ess_file, "\t")]

decoys = files.read_fasta(decoy_file)
decoys = {id: decoys.sequences[i] for i, id in enumerate(decoys.ids)}

all = files.read_fasta(output_file)
non_decoys = {id: all.sequences[i] for i, id in enumerate(all.ids) if id not in decoys}


cds_entries = files.read_fasta(transcripts_file)
cds_entries = {id: cds_entries.sequences[i] for i, id in enumerate(cds_entries.ids)}

seq = "ATCAGCAGTCAG"
query = "GCA"
index = seq.index(query)
print((index + len(query)) % 3)



def in_frame_densities(seq_list, frame):
    hit_count = collections.Counter()
    for i, id in enumerate(seq_list):
        # if i < 10:
        cds = cds_entries[id.split(".")[0]]
        splice_seq = seq_list[id]
        exon_seq = splice_seq[:length].upper()
        intron_seq = splice_seq[length:].upper()

        index = cds.index(exon_seq)
        intron_frame = (index + len(exon_seq)) % 3

        hits = calc_hits(intron_seq, ["TAA", "TAG", "TGA"])
        for hit in hits:
            if (hit+intron_frame) % 3 == frame:
                hit_count[hit] += 1
                hit_count[hit+1] += 1
                hit_count[hit+2] += 1
    hit_density = {i: np.divide(hit_count[i], len(seq_list)) if i in hit_count else 0 for i in range(length)}

    return hit_density

decoys_0 = in_frame_densities(decoys, 0)
decoys_1 = in_frame_densities(decoys, 1)
decoys_2 = in_frame_densities(decoys, 2)
# non_decoys_in_frame = in_frame_densities(non_decoys)
with open("decoy_stop_frames.csv", "w") as outfile:
    outfile.write("pos,frame_0,frame_1,frame_2\n")
    [outfile.write("{0},{1},{2},{3}\n".format(i, decoys_0[i], decoys_1[i], decoys_2[i])) for i in sorted(decoys_0)]



stop_motifs = [i for i in motifs if len(re.findall("(TAA|TAG|TGA)", i)) > 0]
non_stop_motifs = [i for i in motifs if i not in stop_motifs]

# stop_decoy_densities = calc_densities(decoys, stop_motifs)
# non_stop_decoy_densities = calc_densities(decoys, non_stop_motifs)
#
#
# with open("stop_non_stop_decoys.csv", "w") as outfile:
#     outfile.write("pos,stop,non_stop\n")
#     [outfile.write("{0},{1},{2}\n".format(i, stop_decoy_densities[i], non_stop_decoy_densities[i])) for i in sorted(stop_decoy_densities)]
#


# decoy_densities = calc_densities(decoys, motifs)
# non_decoy_densities = calc_densities(non_decoys, motifs)

# decoy_stops = calc_densities(decoys, ["TAA", "TAG", "TGA"])
# non_decoy_stops = calc_densities(non_decoys, ["TAA", "TAG", "TGA"])

# density_file = "densities.csv"
# with open(density_file, 'w') as outfile:
#     outfile.write("pos,decoy,non_decoy\n")
#     [outfile.write("{0},{1},{2}\n".format(i, decoy_densities[i], non_decoy_densities[i])) for i in sorted(decoy_densities)]
#
# stop_file = "stops.csv"
# with open(stop_file, 'w') as outfile:
#     outfile.write("pos,decoy,non_decoy\n")
#     [outfile.write("{0},{1},{2}\n".format(i, decoy_stops[i], non_decoy_stops[i])) for i in sorted(decoy_stops)]
