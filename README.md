# bamdam-extras

If you're looking for the bamdam software, see [here](https://github.com/bdesanctis/bamdam#).

This is an informal collection of standalone ancient DNA scripts that I wrote for myself and collaborators. It is not extensively tested, not formally supported, and may break. If you find it useful, please cite the [the bamdam paper](https://doi.org/10.1186/s13059-025-03879-x) because everything here is just repurposed [bamdam](https://github.com/bdesanctis/bamdam#) code.

To run, clone this github repo, make executable, and run the python scripts as standalone functions. For example:
```
git clone https://github.com/bdesanctis/bamdam-extras.git
cd bamdam-extras
chmod +x *.py
./compute-nolca.py --h
```

Requires python v 3.8 or higher. Feel free to contribute or make (straightforward) requests. You can reach me at bddesanctis@gmail.com. 

All functions require bams to have MD tags. If your bam doesn't have those, you can first add them with samtools:
```
samtools calmd -b in.bam reference.fa > out.bam
```
Also bams must be read-sorted.

## Example tasks

Q: I want statistics about my bam like you'd get from bamdam compute, but I don't have an LCA file (e.g. it's all one taxon, or I want to know stats about each reference separately).   
A: Use compute-nolca. If you want separate stats for each ref ID, use --per-contig.

Q: I want a damage plot from a bam, but I don't have an LCA file.   
A: Use compute-nolca first, and then feed the output subs file into plotdamage-nolca. This will essentially replicate mapDamage behavior without requiring a reference fasta or a single reference. 

Q: I want a damage plot from a bam, I don't have an LCA file, and I want a separate line for each reference ID. For example: I have reads mapped back against a metagenome-assembled genome I made, and I want a separate damage line for each contig (like pyDamage).    
A: Use compute-nolca with --per-contig, then feed the output subs file into plotdamage-nolca with --per-contig. 

Q: I want to annotate a bam with PMD scores, like PMDtools, but I don't have an LCA file.   
A: Use the pmd function. PMDtools, as of the time of writing, has a serious bug in single-stranded mode, whereas this function does not. 

Q: I want a damage plot, but I want to plot all the substitution types, not lump them in "other".   
A: Use plotallsubs on any subs file. You can specify --tax or not. 

Q: I hate merging, read-sorting, and stripping unnecessary headers from my bams.   
A: Chenxi Zhou wrote a lovely function to do this faster and with less RAM than samtools merge/sort, and it also strips your bam of unused headers on the way. You can find it implemented in Richard Durbin's onebam as "onebam bamsort". https://github.com/richarddurbin/onebam

## Usage

### compute-nolca
  ```
usage: compute-nolca.py [-h] --in_bam IN_BAM --out_tsv OUT_TSV --out_subs OUT_SUBS --stranded
                        {ss,ds} [--k K] [--mode {1,2,3}] [--per-contig] [--mincount MINCOUNT]
                        [--show_progress]

Compute bamdam statistics without an LCA file

optional arguments:
  -h, --help           show this help message and exit
  --in_bam IN_BAM      Path to the BAM file (required)
  --out_tsv OUT_TSV    Path to the output tsv file (required)
  --out_subs OUT_SUBS  Path to the output subs file (required)
  --stranded {ss,ds}   Either ss for single stranded or ds for double stranded (required)
  --k K                Value of k for per-node counts of unique k-mers and duplicity
                       (default: 29)
  --mode {1,2,3}       Mode to calculate stats. 1: use best alignment (recommended), 2:
                       average over reads, 3: average over alignments (default: 1)
  --per-contig         Treat each reference ID as a separate node (default: not set)
  --mincount MINCOUNT  Minimum read count to keep a contig; only used with --per-contig
                       (default: 100)
  --show_progress      Print a progress bar (default: not set)
  ```
  
The output tsv here is standard bamdam format so you can use it with bamdam krona and combine, but some of the columns are obviously irrelevant now. See normal bamdam compute docs for futher info. 

The --per-contig flag will trigger --mincount 100 by default, but you can edit it if you want. This is so you don't get insane damage lines from contigs with very few reads in your downstream plotdamage-nolca with --per-contig. This flag will also throw an error if you try to change the mode, because that doesn't make any sense.

### plotdamage-nolca
  ```
usage: plotdamage-nolca.py [-h] --in_subs IN_SUBS [--tax TAX] [--outplot OUTPLOT]
                           [--ymax YMAX] [--per-contig]

Plot damage from a subs file (designed for compute-nolca output)

optional arguments:
  -h, --help         show this help message and exit
  --in_subs IN_SUBS  Input subs file from compute-nolca or bamdam compute (required)
  --tax TAX          Taxonomic node ID or name (optional; defaults to first row in file)
  --outplot OUTPLOT  Filename for output plot, ending in .png or .pdf (default:
                     damage_plot.png)
  --ymax YMAX        Maximum for y axis (optional, default: auto)
  --per-contig       Plot all contigs/nodes from the subs file as separate lines
  ```
  
Acts similarly to bamdam plotdamage, but takes in only one subs file, and plots a C>T, G>A and Other line for every "tax" (or contig, in the MAG use case discussed above).  Specifying a tax ID is optional, and incompatible with the --per-contig flag.

### pmd
```
usage: pmd.py [-h] --in_bam IN_BAM --out_bam OUT_BAM --stranded {ss,ds} [--show_progress]

Annotate a BAM file with PMD scores

optional arguments:
  -h, --help          show this help message and exit
  --in_bam IN_BAM     Input BAM file (required). Must have MD tags.
  --out_bam OUT_BAM   Output BAM file with PMD scores in DS:Z tag (required)
  --stranded {ss,ds}  Library prep: ss for single-stranded, ds for double-stranded (required)
  --show_progress     Show progress bar (requires tqdm)
```
  
Like PMDtools but without a bug in single-stranded mode and written with bamdam code. Note that there are some parameters involved in the PMD calculation, which are set by default to the same ones used in the PMDtools software. You can edit the source code easily if you want to fiddle with them; they are defined at the top with some info about what they are. 

### plotallsubs

```
usage: plotallsubs.py [-h] --in_subs IN_SUBS --tax TAX [--outplot OUTPLOT] [--ymax YMAX]

Plot all substitution types from a bamdam subs file

optional arguments:
  -h, --help         show this help message and exit
  --in_subs IN_SUBS  Input subs file from bamdam compute (required)
  --tax TAX          Taxonomic node ID (required)
  --outplot OUTPLOT  Filename for output plot, ending in .png or .pdf (default:
                     allsubs_plot.png)
  --ymax YMAX        Maximum for y axis (optional, default: auto)
```

Basically the same as bamdam plotdamage, but plots a line for every type of substitution, not just the C>T, G>A and Other in the normal bamdam plotdamage command. Can only take in one subs file because this would be a totally overwhelming number of lines with multiple samples.










