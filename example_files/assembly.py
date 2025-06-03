#!/bin/python3.13
import os, subprocess, sys, argparse, signal,csv, time
from datetime import timedelta

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
        
signal.signal(signal.SIGHUP, signal.SIG_IGN) #ignore bash quit
parser = argparse.ArgumentParser(
                    prog='Assembly quality',
                    description='Assess the quality of an assembly',
                    epilog='Made by Guilhem Zeitoun, by and from IMGT team')
parser.add_argument('-m', '--minimap', required=True)      # option that takes a value
parser.add_argument('-l', '--locusfile', required=True)      # option that takes a value
parser.add_argument('-s', '--species', required=True)      # option that takes a value
parser.add_argument('-g', '--genelist', required=False)      # option that takes a value
parser.add_argument('--totalread', required=False, action='store_true')
parser.add_argument('-q', '--fastqc', required=False, action='store_true')
arg = parser.parse_args()
print("Répertoire de travail actuel :", os.getcwd())
if not os.path.exists("assembly.fasta"):
    print("Assembly is not found. Exiting")
    sys.exit(1)
starttime = time.monotonic()
reads=[]
bamreads=[]
for file in os.listdir(os.getcwd()):
    if file.startswith("SRR") and file.endswith(".bam"):
        bamreads.append(file)
if len(bamreads)>0:
    print("BAM files were found, converting")
    for file in bamreads:
        for path in execute(f"samtools fastq -@ 32 \"{file}\" > \"{file}.fastq.gz\" && rm {file}"):
            print(path, end="",flush=True)
for file in os.listdir(os.getcwd()):
    if file.startswith("SRR") and file.endswith(".gz") or file=="full.fastq.gz":
        reads.append(file)
if len(reads)==0:
    print("No reads found, have you downloaded them? Exiting")
    sys.exit(1)
elif reads[0]=="full.fastq.gz":
    pass
elif len(reads)==1:
    for path in execute("mv " + reads[0] + " full.fastq.gz"):
        print(path, end="",flush=True)
else:
    for path in execute("cat SRR*.gz > full.fastq.gz && rm SRR*.gz"):
        print(path, end="",flush=True)
print("Reads found, continuing")
if arg.fastqc:
	print("Analysis in progress")
	for path in execute("fastqc -t 28 full.fastq.gz"):
            print(path,end="",flush=True)
	print("Analysis finished, exiting")
	sys.exit(0)
else:
	print("No analysis done, continuing")
if arg.minimap != "map-pb" and arg.minimap != "map-ont" and arg.minimap != "lr:hq" and arg.minimap != "map-hifi":
    print("You can have minimap as map-pb, map-ont, map-hifi or lr:hq. Choose carefully and retry")
    sys.exit(1)
print("Alignment starting (minimap)")
if os.path.exists("full_cs.bam") or os.path.exists("full_cs.sam"):
    print("Full cs bam exists, skipping")
else:
    for path in execute("minimap2 -ax " + arg.minimap + " -t 14 --cs --eqx assembly.fasta full.fastq.gz > full_cs.sam"):
       print(path, end="",flush=True)
print("\nMinimap success")
print("\nSamtools starting")
if os.path.exists("full_cs.bam"):
    print("Full cs bam exists, skipping")
else:
    for path in execute("samtools sort -@ 28 full_cs.sam -o full_cs.bam"):
      print(path, end="",flush=True)
    for path in execute("samtools index -c -@ 28 full_cs.bam"):
      print(path, end="",flush=True)
    for path in execute("rm full_cs.sam"):
      print(path, end="",flush=True)
print("Samtools exiting")
print("Starting pileups...")
locus=[]
with open(arg.locusfile, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in spamreader:
        if len(row)!=5:
            print(f"Row does not contain 5 tab-separated fields, found " + ",".join(row))
        gene=row[0].replace(" ","_")
        haplo="pri" if row[1].strip().lower() =="primary" else "alt"
        contig=row[2]
        start=int(row[3])
        end=int(row[4])
        if start>end:
            end,start=start,end
        locus.append(f"{contig}:{start}-{end}")
        for path in execute("samtools mpileup -Q 0 -q 0 --ff SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -r " + contig + ":" + str(start) + "-" + str(end) + " -aa -f assembly.fasta full_cs.bam > " + gene + "_" + haplo + "_pileup.txt"):
            print(path, end="",flush=True)
print("End pileups")
print("Create compact BAM with only interested loci")
for path in execute("samtools view -b -h full_cs.bam \"" + "\" \"".join(locus) + "\" > locus_" + locus + "_" + haplo + ".bam && samtools index -c locus.bam"):
    print(path, end="",flush=True)
print("\nStart IMGT_assembly")
totalread=""
if arg.totalread:
	totalread="--totalread"
genelist=""
if arg.genelist:
        genelist="-g \"" + arg.genelist + "\""
paths=f"IMGT_StatAssembly -f full_cs.bam {totalread} {genelist} -s {arg.species} --locuspos {arg.locusfile} --outdir results/"
for path in execute(paths):
    print(path,end="",flush=True)
print("End IMGT assembly")
endtime=time.monotonic()
time_elapsed = timedelta(seconds=endtime-starttime)
print(f"Done in {time_elapsed}")
