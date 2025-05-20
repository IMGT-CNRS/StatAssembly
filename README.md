# IMGT/StatAssembly
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.15016809.svg)](https://zenodo.org/doi/10.5281/zenodo.15016809)
![GitLab Release](https://img.shields.io/gitlab/v/release/imgt-igh%2Fstatassembly?gitlab_url=https%3A%2F%2Fsrc.koda.cnrs.fr)
![GitLab License](https://img.shields.io/gitlab/license/imgt-igh%2Fstatassembly?gitlab_url=https%3A%2F%2Fsrc.koda.cnrs.fr%2F)


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

> [!TIP]
> If your BAM file does not contain equal CIGAR format or a CS/MD tag, you can recalculate this tag without relaunching the analysis completely if you have the bam file and the assembly like:
> ```console
> samtools calmd -b -@ 28 full.bam assembly.fasta > align.bam
> ```
> Then you can use the new BAM file to have full results in the software.

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

The rest is ***case sensitive***. You can only have one alternate per primary (the line just after the primary) and as many primary as you want. Primary and Alternate are compared and shown together in graphs.

Contig name has to match reference ID, start and end should match SAM regions (1-based position). If start is greater than end, the locus would be considered reverse.
Example in test files.
* A CSV file containing gene position on the chromosome (optional). ***Header must be preserved***:
```csv
"gene","chromosome","strand","start","end"
"IGHA1","NC_060938.1","minus","99976277","99980553"
"IGHA2","NC_060938.1","minus","99837189","99841426"
"IGHD","NC_060938.1","minus","100108615","100117138"
```
Strand can be 0 or 1 (reverse), + or - (reverse), plus or minus (reverse). Strand is related to the chromosome.
Chromosome, start and end should be 1-based position. *Start should be less than end*.
Example in test files.

### Generation of a BAM file

To generate the BAM file used in the analysis, you can follow those steps. Keep in mind minimap2 requires at least 32 GB of memory.

* Download the assembly of T2T-CHM13v2.0 from [NCBI website](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and name it CHM13v2.0.fasta.
* Download HiFi reads of T2T-CHMv2.0 from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/?term=SRX789768*+CHM13) and name it reads.fastq.gz (compress it if needed).
* Install dependancies if not existing
```console
# apt install minimap2 samtools
```
* Launch minimap from bash terminal and create the BAM file and its index:
```bash
minimap2 -ax map-hifi -t 32 --eqx --cs CHM13v2.0.fasta reads.fastq.gz > reads.sam
samtools sort -@ 32 -o reads.bam reads.sam
samtools index -c reads.bam
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

Download the binaries from binaries folder or releases depending on your OS and put it in your path. Then type:
```bash
IMGT_StatAssembly -h
```
to access the help.

### Source code

- [ ] Install [rust](https://www.rust-lang.org/fr/learn/get-started) if not installed.
- [ ] Check Rust version `rustc -V`, should be >= 1.85.
- [ ] Do a `git clone` and then `cargo build --release` to compile the software.


## Execution  (Test)

Here is the command to execute with example files from the repository folder on linux 64bits:
```bash
binaries/IMGT_StatAssembly_linux_x64_86 -f example_files/CHM13v2.0.bam -s human -l example_files/CHM13v2.0loc.csv -g example_files/CHM13v2.0geneloc.csv -o results/
```

## Output

The expected output from execution is present in `example_files/results/`.

### Description of generated file in example folder

- *break.txt* lists where breaks are present (default: 3, parameter: breaks).
- *mismatchresult.txt* shows two graphs.
    - The first graph shows the PHRED quality score (`rgb(0,0,0)` (black) curve) with the legend on the right axis. The rate of mismatches (`rgb(126,87,194)`) and misalign (`rgb(239,83,80)`) is also shown for each position with the legend on the left axis. A misalign is a read that has an indel at this position and a mismatch a read with a substitution.
    - The bottom graph shows the number of mismatch rate for all reads which alignment cover the position indicated (`rgb(255, 171, 145)`).
- *readresult.png* shows over the locus (position on the chromosome and on the locus displayed) the number of reads based on their quality score, as well as secondary, supplementary and overlapping reads. The number of breaks is displayed as red bars at the bottom panel if existing.
- *positionresult.csv* lists all the information of both graphs. However mismatches and misalign represents a number and not a rate as in the graph. The rate could be recalculated by dividing with the sum of reads in the column map60,map1 and map0.
- If gene list is provided:
    - *allele_confidence.csv*: List all suspicious (shown as ! in Excel and `rgb(239,83,80)` on charts) and warning positions (shown as ~ in Excel and `rgb(255, 183, 77)` on charts). By default:
        - Warning positions (`rgb(255, 183, 77)`) are positions where less than x reads (parameter: minreads default 10) are present and/or the rate of reads matching the base compared to the number of reads present at this position is above the suspicious position rate and below the treeshold (parameter: percentwarning default 0.8).
        - Suspicious positions (`rgb(239,83,80)`) are positions where the rate of reads matching the base compared to the number of reads present at this position is less than the treeshold (parameter: percentalerting default 0.6).
    - A folder containing a graph for each gene, with number of total reads for each position (total reads), reads without indels (sequence match) and sequence match. The number of reads that covers the entire region with 100% match are displayed with the `rgb(0,0,0)` (black) curve.
    - *geneanalysis.csv*: List all genes, their chromosome, strand, start and end. It displays the average read coverage (how many times larger the reads are compared to the length of the given region), the number of reads on this region. Then for each position, the number of reads in total with the number of reads with identical sequence (=), ones with substitutions (X) and ones with indels (ID). Readsfull column counts the number of reads spanning the entire region, whereas reads100 and reads100m shows respectively the number of reads matching without indels or with perfect match the full region. Coveragex shows how much position are covered by at least x reads (default: 10, parameter: coverage).

### Results analysis

For a better overview of IMGT rules based on this result, check [IMGT assembly quality rules](https://imgt.org/IMGTScientificChart/Assemblies/IMGTassemblyquality.php).

## How to cite

If you use IMGT/StatAssembly in your work, please cite the version you used, for example:

> Institut de Génétique Humaine. (2025). IMGT/StatAssembly (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.15396812

## License

IMGT/StatAssembly: Copyright Guilhem Zeitoun (IMGT), 2025, licensed under the EUPL (European Union Public Licence) v1.2.
The IMGT logo and the software logo remain the property of IMGT and all rights are reserved.

The [Rust crab](https://www.rustacean.net/) is under [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/).


## Memory consumption
The script uses hundreds of Mo up to some Gb for a several Mo locus. Some Gb of memory should be reserved depending on the BAM file.
![Memory consumption](/images/memory.png)