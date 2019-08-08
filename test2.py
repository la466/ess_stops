import collections
import bioUtilities.files as files
import bioUtilities.seq as seq
from maxentpy import maxent  # use normal version of maxent
from maxentpy.maxent import load_matrix5, load_matrix3
import re
import numpy as np

ess_file = "source_data/motif_sets/ess_fas_hex3.txt"

introns = files.read_fasta("source_data/human/introns.fasta")
introns = {id.split("-")[0]: introns.sequences[i] for i, id in enumerate(introns.ids)}




decoys = files.read_fasta("decoys.fa")
decoys = {id: decoys.sequences[i] for i, id in enumerate(decoys.ids)}

all = files.read_fasta("splice_sites.fa")
non_decoys = {id: all.sequences[i] for i, id in enumerate(all.ids) if id not in decoys}

cds = files.read_fasta("source_data/human/cds_transcripts.fasta")
cds = {id: cds.sequences[i] for i, id in enumerate(cds.ids)}


stops = ["TAA", "TAG", "TGA"]

decoy_introns = [decoys[id][50:] for id in decoys]
non_decoy_introns = [non_decoys[id][50:] for id in non_decoys]

decoy_density = seq.calc_motif_density(decoy_introns, stops)
non_decoy_density = seq.calc_motif_density(non_decoy_introns, stops)
# print(decoy_density)
# print(non_decoy_density)

motif_search = re.compile("(?=({0}))".format("|".join(stops)))

decoy_hits = 0
decoy_nts = 0

frame_hits = collections.Counter()
nts = 0

for i, id in enumerate(non_decoys):
    if i == 3:
    # if i:
        print(id)
        transcript_id = id.split(".")[0]
        cds_seq = cds[transcript_id]
        exon = non_decoys[id][:50].upper()
        intron = non_decoys[id][50:]
        full_intron = introns[id]

        if len(full_intron) % 3 == 0:


            decoy_nts += len(intron)
            index = cds_seq.index(exon)
            exon_end_index = index + len(exon)
            intron_start_frame = exon_end_index % 3

            print(intron_start_frame)

            nts += len(intron)
            for i in range(3):
                print(i, (i+intron_start_frame) % 3)
                frame_hits[(i+intron_start_frame) % 3] += 3*len([i for i in re.findall('.{3}', intron[i:]) if i in stops])

print(frame_hits)
print(nts)

            # print(intron)
            # print(full_intron)
            #
            #
            #
            # matches = re.finditer(motif_search, intron)
            # hits = []
            # [hits.extend(list(range(hit.span()[0], hit.span()[0] + len(hit.group(1))))) for hit in matches]
            # hits = sorted(list(set(hits)))
            # # print(intron_start_frame, hits)
            # for i in range(0, len(hits), 3):
            #     if (hits[i] + intron_start_frame) % 3 == 0:
            #         decoy_hits += 3

print(np.divide(decoy_hits, decoy_nts))
