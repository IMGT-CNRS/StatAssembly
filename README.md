# IMGT/StatAssembly
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.15016809.svg)](https://zenodo.org/doi/10.5281/zenodo.15016810)
![GitHub License](https://img.shields.io/github/license/DorianCoding/IMGT_StatAssembly)

IMGT/StatAssembly uses BAM file to assess the quality of the assembly, including order of genes and validation of alleles in IG and TR loci. 
<p align="middle">
<img src = "images/logo_software.png" width ="200" />
</p>

## Source
It is a script written in Rust, compiled as a optimized binary.
<p align="middle">
<img src = "images/rust.png" width ="50" />
</p>
The script was made by <a href="//www.imgt.org">IMGT team</a>.
<p align="middle">
<img src= "images/logo_imgt.png" width ="150" />
</p>

## Script input files and data
* The BAM file from analysis and its index, the presence of a cigar with `=`/`X` (match; substitution), a MD tag or a cs tag is recommended. *Some analysis won't be available without*.
* A TSV file with the following information, separated by a tabular:
```tsv
Locus Haplotype contig  start end
```
Locus must be one of the following:
* IGH
* IGK
* IGL
* TRA[^1]
* TRB
* TRG

[^1]: TRD is inside TRA locus and so both loci are analyzed together.

Haplotype must be one of the following:
* Primary or pri or p (case insensitive)
* Alternate or alt or a (case insensitive)

The rest is ***case sensitive***. You can only have one alternate per primary (the line just after the primary) and as many primary as you want. Primary and Alternate are compared in graph.

Contig, start and end should match SAM regions (1-based position). If start is greater than end, the locus would be considered reverse.
Example in test files.
* A CSV file containing gene position on the chromosome (optional). ***Header must be preserved***:
```csv
"gene","chromosome","strand","start","end"
"IGHA1","NC_060938.1","1","99976277","99980553"
"IGHA2","NC_060938.1","1","99837189","99841426"
"IGHD","NC_060938.1","1","100108615","100117138"
```
Strand can be 0 or 1 (reverse), + or - (reverse), plus or minus (reverse).
Chromosome, start and end should be 1-based position. Start should be less than end.
Example in test files.

### Generation of a BAM file

To generate the BAM file used in the analysis, you can follow those steps. Keep in mind minimap2 requires at least 32 GB of memory.

* Download the assembly of T2T-CHM13v2.0 from [NCBI website](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and name it CHM13v2.0.fasta.
* Download HiFi reads of T2T-CHMv2.0 from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/?term=SRX789768*+CHM13) and name it reads.fastq.gz (compress it if needed).
* Install dependancies if not existing
```console
# apt install minimap2 samtools
```
* Launch minimap from bash terminal:
```bash
minimap2 -ax map-hifi -t 32 --eqx --cs CHM13v2.0.fasta reads.fastq.gz > reads.sam
samtools sort -@ 32 -o result.bam result.sam
samtools index -c result.bam
```
The output of the commands should look like:
```console
[M::mm_idx_gen::42.277*1.51] collected minimizers
[M::mm_idx_gen::47.081*2.06] sorted minimizers
[M::main::47.081*2.06] loaded/built the index for 24 target sequence(s)
[M::mm_mapopt_update::50.206*1.99] mid_occ = 177
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 24
[M::mm_idx_stat::52.608*1.95] distinct minimizers: 215124360 (92.06% are singletons); average occurrences: 1.455; average spacing: 9.958; total length: 3117275501
[M::worker_pipeline::64.187*2.84] mapped 7641 sequences
[M::main] Version: 2.27-r1193
[M::main] CMD: minimap2 -ax map-hifi -t 32 --eqx --cs genome.fasta reads.fastq.gz
[M::main] Real time: 64.475 sec; CPU: 182.440 sec; Peak RSS: 12.281 GB
[bam_sort_core] merging from 0 files and 32 in-memory blocks...
```
The BAM file as a result is different from the BAM in example_folder as it is restricted to specific portion of the assembly but contains the same results and can be used as input file.

## How to install

### Binaries

Download the binaries from binaries folder depending on your OS. Then type:
```bash
IMGT_StatAssembly -h
```
to access the help.

### Source code

- [ ] Install rust if not installed.
- [ ] Check Rust version, should be >= 1.85.
- [ ] Do a git clone and then `cargo build --release` to compile the software.


## Execution (Test)

Here is the command to execute with example files from the repo folder on linux 64bits:
```bash
binaries/IMGT_StatAssembly_linux_x64_86 -f example_files/CHM13v2.0.bam -s human -l example_files/CHM13v2.0loc.csv -g example_files/CHM13v2.0geneloc.csv -o results/
```

## Output

The expected output from execution is present in `example_files/results/`.

### Description of generated file in example folder

- *break.txt* lists where breaks are present (default: 3,parameter: breaks). An empty file means no break.
- *mismatchresult.txt* shows two graphs.
    - The first graph shows the PHRED quality score (`rgb(0,0,0)` (black) curve) with the right axis. The rate of mismatches (`rgb(126,87,194)`) and misalign (`rgb(239,83,80)`) is also shown for each position. A misalign is a read that has an indel at this position and a mismatch a read with a substitution.
    - The bottom graph shows the number of mismatch rate for all reads overlapping the position indicated (`rgb(255, 171, 145)`).
- *readresult.png* shows over the locus (position on the chromosome and on the locus displayed) the number of reads based on their quality score, as well as secondary, supplementary and overlapping reads. The number of breaks is displayed as red bars at the bottom panel if existing.
- *positionresult.csv* lists all the information of both graphs. However mismatches and misalign represents a number and not a rate as in the graph.
- If gene list given:
    - *allele_confidence.csv*: List all suspicious (shown as ! in Excel and `rgb(239,83,80)` on charts) and warning positions (shown as ~ in Excel and `rgb(255, 183, 77)` on charts).
    - A folder containing a graph for each gene, with number of total reads for each position (total reads), reads without indels (sequence match) and sequence match. The number of reads on the entire region with 100% match are displayed with the `rgb(0,0,0)` (black) curve.
    - *geneanalysis.csv*: List all genes, their chromosome, strand, start and end. It displays the average read coverage (how many times larger the reads are compared to the length of the given region), the number of reads on this region. Then for each position, the number of reads in total with the number of reads with identical sequence (=), ones with substitutions (X) and ones with indels (ID). Readsfull counts the number of reads spanning the entire region, whereas reads100 and reads100m shows respectively the number of reads matching without indels or with perfect match the full region. Coveragex shows how much position are covered by at least x reads (default: 10, parameter: coverage).

### Results analysis

For a better overview of IMGT rules based on this result, check [IMGT assembly quality rules](https://imgt.org/IMGTScientificChart/Assemblies/IMGTassemblyquality.php).

## How to cite

If you use IMGT/StatAssembly in your work, please cite the version you used, for example:

> Institut de Génétique Humaine. (2025). IMGT StatAssembly (v0.1.7). Zenodo. https://doi.org/10.5281/zenodo.15234695

## Memory consumption
The script uses hundreds of Mo up to some Gb for a several Mo locus. Some Gb of memory should be reserved depending on the BAM file.
![Memory consumption](/images/memory.png)