import subprocess

args = ["pip", "install", "git+https://github.com/la466/bioUtilities", "--upgrade", "-t", "."]
process = subprocess.Popen(args, shell = False)
stdout, stderr = process.communicate()
