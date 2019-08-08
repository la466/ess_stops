import collections

from bioUtilities.files import fasta_from_bed, read_many_fields
from bioUtilities.bed import convert_chr_name

input_file = "source_data/alternative_5_splice_site_exons_chr.bed"
entries_file = "source_data/alternative_5_splice_site_exons.bed"
output_fasta = "results/alternative_5_splice_site_exons.fa"
genome_fasta = "../source_data/Genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# convert_chr_name(input_file, entries_file)

entries = read_many_fields(input_file)


entry_list = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
for entry in entries[1:]:
    entry_list[entry[0]][int(entry[2])].append(entry)

output_ids = collections.Counter()

with open(entries_file, "w") as outfile:
    outfile.write("#chr\tstart\tend\tid\t.\tstrand\t{0}\n".format("\t".join(entries[0][4:-1])))
    for chr in sorted(entry_list):
        print(chr)
        for start in sorted(entry_list[chr]):
            for entry in  entry_list[chr][start]:
                entry_chr = entry[0].strip("chr")
                strand = entry[1]
                start = entry[2]
                end = entry[3]
                id = entry[-1]
                info = entry[4:-1]

                output_ids[id] += 1

                output = [entry_chr, start, end, "{0}.{1}".format(id, output_ids[id]), ".", strand] + info
                outfile.write("{0}\n".format("\t".join(output)))


fasta_from_bed(entries_file, genome_fasta, output_fasta, names = True)
