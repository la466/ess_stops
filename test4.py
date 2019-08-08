import collections
import bioUtilities.files as files
import bioUtilities.seq as seq
import re

gtf = "../source_data/Genomes/hg38/Homo_sapiens.GRCh38.94.gtf"

entries = [i for i in files.read_many_fields(gtf, "\t") if not i[0].startswith('#')]


genes = collections.defaultdict(lambda: collections.defaultdict(lambda: []))


for i in entries[:5000]:
    type = i[2]

    print(i)

    if type == "exon":
        info = i[-1]
        try:
            gene_id = re.findall('gene_id "(.*?)"', info)[0]
            transcript_id = re.findall('transcript_id "(.*?)"', info)[0]
            biotype = re.findall('transcript_biotype "(.*?)"', info)[0]
            exon_number = int(re.findall('exon_number "(.*?)"', info)[0])

            if biotype == "protein_coding":
                # print(gene_id, biotype)
                # print(info)

                genes[gene_id][transcript_id].append(exon_number)
        except:
            pass


for gene_id in genes:
    print(gene_id)
    for transcript_id in genes[gene_id]:
        print(transcript_id, sorted(genes[gene_id][transcript_id]))
    print("\n")
