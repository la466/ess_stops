import ftp_ops
from bioUtilities.files import read_many_fields
from bioUtilities.dir import create_directory

link = "ftp://ftp.ensembl.org/pub/release-96/bamcov/homo_sapiens/genebuild/"

output_dir = "human_bodymap_bams"
create_directory(output_dir)

files = [i[0] for i in read_many_fields("filelist.txt", "\t")]

host = "ftp.ensembl.org"
user = None
password = None

needed_dir = "/pub/release-96/bamcov/homo_sapiens/genebuild/"

ftp = ftp_ops.ftp_connect(host, user, password, directory = needed_dir)

for file in files:
    ftp = ftp_ops.ftp_retrieve(ftp, host, user, password, needed_dir, file, destination = output_dir)
