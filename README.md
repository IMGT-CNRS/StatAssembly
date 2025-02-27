# IMGT/StatAssembly
IMGT/StatAssembly uses BAM file to assess the quality of the assembly, including order of genes and validation of alleles in IG and TR loci. It is a script written in Rust, compiled as a optimized binary.
<p align="middle">
<img src = "images/rust.png" width ="50" />
</p>
The script was made by <a href="//www.imgt.org">IMGT team</a>.
<p align="middle">
<img src= "images/logo_imgt.png" width ="150" />
</p>

## Script input files and data
* The BAM file from analysis and its index, the presence of a cigar with `=`/`X` (match; substitution), a MD tag or a cs tag is recommended. *Some analysis won't be available without*.
* A CSV file with the following information, separated by a tabular:
```
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
* Primary
* Alternate

It is ***case sensitive***. You can only have one alternate per primary (the line just after the primary) and as many primary as you want. Primary and Alternate are compared in graph.

Contig, start and end should match SAM regions.
Example in test files.

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

## Execution

Here is the command to execute with example files from the repo folder on linux 64bits:
```bash
binaries/IMGT_StatAssembly_linux_x64_x86 -f example_files/CHM13v2.0.bam -s human -l example_files/CHM13v2.0loc.csv -g example_files/CHM13v2.0geneloc.csv -o results/
```

## Output



### Memory consumption
The script uses hundreds of Mo up to some Gb for a several Mo locus. Some Gb of memory should be reserved depending on the BAM file.
![Memory consumption](/images/memory.png)