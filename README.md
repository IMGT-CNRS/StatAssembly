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
* The BAM file from analysis, the presence of a cigar with `=`/`X` (match; substitution) is compulsory.
* A CSV file with the following information, separated by a tabular:
```
Locus Haplotype contig  start end
```
Locus must be IGH,IGK,IGL,TRA,TRB or TRG. Haplotype must be primary or alternate (case sensitive). Contig, start and end should match SAM regions.
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
- [ ] Check Rust version, should be >= 1.81.
- [ ] Do a git clone and then `cargo build --release` to compile the software.

## Output

### Memory consumption
The script uses hundreds of Mo up to some Gb for a several Mo locus. Some Gb of memory should be reserved depending on the BAM file.
![Memory consumption](/images/memory.png)