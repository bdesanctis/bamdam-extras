#!/usr/bin/env python3

"""
pmd.py - Annotate a BAM file with PMD scores

Part of bamdam-extras. Please cite the bamdam paper if you use this.
https://github.com/bdesanctis/bamdam

PMD scores are from Skoglund et al. 2014 (doi:10.1073/pnas.1318934111).
Note: PMDtools has a bug in single-stranded mode; this implementation is correct.
Details of PMDtools bug: in lines 868 and 900 of the main python script, it multiplies the likelihood by the double-stranded models, and then there is an if statement that multiplies by the single-stranded models resulting in totally incorrect pmd scores for single-stranded mode. the reimplementation here should be correct. You can of course just fix the original too.

Usage:
    python pmd.py --in_bam FILE.bam --out_bam FILE.pmd.bam --stranded ds

The output BAM will have PMD scores in the DS:Z tag.
"""

import sys
import argparse
import math

try:
    import pysam
    pysam_imported = True
except ImportError:
    pysam_imported = False


# PMD score parameters (from PMDtools) -- see the original paper if you want to understand/fiddle with these. They are tuned for the neanderthal implementation
# P = probability of damage at position 1
# C = background substitution rate
# pi = probability of sequencing error causing apparent damage
# Dn = null model damage rate (PMDtools code uses 0.001, not 0 as stated in paper, which I think is reasonable)
PMD_P = 0.3
PMD_C = 0.01
PMD_PI = 0.001
PMD_DN = 0.001


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate a BAM file with PMD scores"
    )
    parser.add_argument(
        "--in_bam",
        required=True,
        help="Input BAM file (required). Must have MD tags."
    )
    parser.add_argument(
        "--out_bam",
        required=True,
        help="Output BAM file with PMD scores in DS:Z tag (required)"
    )
    parser.add_argument(
        "--stranded",
        required=True,
        choices=["ss", "ds"],
        help="Library prep: ss for single-stranded, ds for double-stranded (required)"
    )
    parser.add_argument(
        "--show_progress",
        action="store_true",
        help="Show progress bar (requires tqdm)"
    )
    return parser.parse_args()


def parse_cigar(cigar):
    """Parse CIGAR string into list of operations."""
    import re
    return re.findall(r'(\d+)([MIDNSHP=X])', cigar)


def parse_md(md):
    """Parse MD tag into list of elements."""
    import re
    return re.findall(r'(\d+|\^[A-Za-z]+|[A-Za-z])', md)


def rev_complement(seq):
    """Return reverse complement of a sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'N': 'N', '-': '-', 's': 's'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def get_mismatches(seq, cigar, md):
    """
    Parse read sequence, CIGAR, and MD tag to reconstruct alignment.
    Returns [mismatch_list, ref_seq, read_seq].
    """
    cig = parse_cigar(cigar)
    md_list = parse_md(md)
    
    ref_seq = ''
    read_seq = ''
    query_pos = 0
    read_pos = 0
    
    for bases_str, cat in cig:
        bases = int(bases_str)
        
        if cat == 'H':
            continue
        
        if cat == 'S':
            md_extra = "softclip" + str(bases)
            if (bases_str + cat) == cig[-1][0] + cig[-1][1]:
                md_list.append(md_extra)
            elif (bases_str + cat) == cig[0][0] + cig[0][1]:
                md_list.insert(0, md_extra)
            
            read_seq += 's' * bases
            ref_seq += 's' * bases
            query_pos += bases
            read_pos += bases
            
        elif cat in ['M', '=', 'X']:
            ref_seq += seq[query_pos:query_pos + bases]
            read_seq += seq[read_pos:read_pos + bases]
            query_pos += bases
            read_pos += bases
            
        elif cat in ['D', 'N']:
            ref_seq += 'N' * bases
            read_seq += '-' * bases
            
        elif cat == 'I':
            read_seq += seq[read_pos:read_pos + bases]
            ref_seq += '-' * bases
            query_pos += bases
            read_pos += bases
    
    rec_pos = 0
    mismatch_list = []
    
    for x in md_list:
        if x.startswith('softclip'):
            num = int(x[len("softclip"):])
            rec_pos += num
        elif x.startswith('^'):
            num_chars_after_hat = len(x) - 1
            rec_pos += num_chars_after_hat
        else:
            if x.isdigit():
                rec_pos += int(x)
            else:
                refhere = x
                readhere = read_seq[rec_pos] if rec_pos < len(read_seq) else 'N'
                char_list = list(ref_seq)
                if rec_pos < len(char_list):
                    char_list[rec_pos] = x
                    ref_seq = "".join(char_list)
                
                read_up_to_here = read_seq[0:rec_pos]
                read_pos_actual = len(read_up_to_here) - read_up_to_here.count("-")
                mismatch_list.append([refhere, readhere, read_pos_actual + 1])
                rec_pos += 1
    
    return [mismatch_list, ref_seq, read_seq]


def calculate_pmd(read, stranded):
    """
    Calculate PMD score for a pysam read object.
    
    This is a reimplementation that fixes the PMDtools bug in single-stranded mode.
    """
    seq = read.query_sequence
    cigar = read.cigarstring
    
    try:
        md = read.get_tag('MD')
    except KeyError:
        return None  # No MD tag
    
    rawphred = list(read.query_qualities) if read.query_qualities else [30] * len(seq)
    flagsum = read.flag
    
    P = PMD_P
    C = PMD_C
    pi = PMD_PI
    Dn = PMD_DN
    
    # Check if reverse strand
    backwards = (flagsum & 16) != 0
    
    mmsc, ref_seq, read_seq = get_mismatches(seq, cigar, md)
    
    # Adjust phred for indels
    phred = []
    phred_idx = 0
    for base in read_seq:
        if base in ['-', 'N']:
            phred.append(0)
        else:
            if phred_idx < len(rawphred):
                phred.append(rawphred[phred_idx])
                phred_idx += 1
            else:
                phred.append(30)
    
    if backwards:
        ref_seq = rev_complement(ref_seq)
        read_seq = rev_complement(read_seq)
        phred = phred[::-1]
    
    refseqlist = list(ref_seq)
    readseqlist = list(read_seq)
    readlength = len(readseqlist)
    
    pmd_lik = 1.0
    null_lik = 1.0
    pos = 0
    
    if stranded == "ss":
        for b in range(readlength):
            if readseqlist[b] == "-":
                continue
            
            if refseqlist[b] == "C" and readseqlist[b] in ["T", "C"]:
                epsilon = (1/3) * 10**(-phred[b] / 10)
                z = pos + 1
                y = readlength - pos
                Dz = ((1-P)**(z-1)) * P + C
                Dy = ((1-P)**(y-1)) * P + C
                
                if readseqlist[b] == "T":
                    pmd_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dz)*(1-Dy) + 
                                    (1-pi)*epsilon*Dz*(1-Dy) + 
                                    (1-pi)*epsilon*Dy*(1-Dz) + 
                                    pi*epsilon*(1-Dz)*(1-Dy))
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn)*(1-Dn) + 
                                     (1-pi)*epsilon*Dn*(1-Dn) + 
                                     (1-pi)*epsilon*Dn*(1-Dn) + 
                                     pi*epsilon*(1-Dn)*(1-Dn))
                else:  # C
                    pmd_lik *= ((1-pi)*(1-epsilon)*(1-Dz)*(1-Dy) + 
                                (1-pi)*epsilon*Dz*(1-Dy) + 
                                (1-pi)*epsilon*Dy*(1-Dz) + 
                                pi*epsilon*(1-Dz)*(1-Dy))
                    null_lik *= ((1-pi)*(1-epsilon)*(1-Dn)*(1-Dn) + 
                                 (1-pi)*epsilon*Dn*(1-Dn) + 
                                 (1-pi)*epsilon*Dn*(1-Dn) + 
                                 pi*epsilon*(1-Dn)*(1-Dn))
            pos += 1
    
    elif stranded == "ds":
        for b in range(readlength):
            if readseqlist[b] == "-":
                continue
            
            if refseqlist[b] == "C" and readseqlist[b] in ["T", "C"]:
                epsilon = (1/3) * 10**(-phred[b] / 10)
                z = pos + 1
                Dz = ((1-P)**(z-1)) * P + C
                
                if readseqlist[b] == "T":
                    pmd_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dz) + 
                                    (1-pi)*epsilon*Dz + 
                                    pi*epsilon*(1-Dz))
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn) + 
                                     (1-pi)*epsilon*Dn + 
                                     pi*epsilon*(1-Dn))
                else:  # C
                    pmd_lik *= ((1-pi)*(1-epsilon)*(1-Dz) + 
                                (1-pi)*epsilon*Dz + 
                                pi*epsilon*(1-Dz))
                    null_lik *= ((1-pi)*(1-epsilon)*(1-Dn) + 
                                 (1-pi)*epsilon*Dn + 
                                 pi*epsilon*(1-Dn))
            
            if refseqlist[b] == "G" and readseqlist[b] in ["A", "G"]:
                epsilon = (1/3) * 10**(-phred[b] / 10)
                z = readlength - pos
                Dz = ((1-P)**(z-1)) * P + C
                
                if readseqlist[b] == "A":
                    pmd_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dz) + 
                                    (1-pi)*epsilon*Dz + 
                                    pi*epsilon*(1-Dz))
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn) + 
                                     (1-pi)*epsilon*Dn + 
                                     pi*epsilon*(1-Dn))
                else:  # G
                    pmd_lik *= ((1-pi)*(1-epsilon)*(1-Dz) + 
                                (1-pi)*epsilon*Dz + 
                                pi*epsilon*(1-Dz))
                    null_lik *= ((1-pi)*(1-epsilon)*(1-Dn) + 
                                 (1-pi)*epsilon*Dn + 
                                 pi*epsilon*(1-Dn))
            
            pos += 1
    
    if pmd_lik == 0 or null_lik == 0:
        return 0.0
    else:
        return math.log(pmd_lik / null_lik)


def main():
    args = parse_args()
    
    if not pysam_imported:
        print("Error: pysam is required. Install with: pip install pysam")
        sys.exit(1)
    
    # Check for tqdm if progress requested
    tqdm = None
    if args.show_progress:
        try:
            from tqdm import tqdm as tqdm_import
            tqdm = tqdm_import
        except ImportError:
            print("Warning: tqdm not available, progress bar disabled.")
    
    # Count alignments for progress bar (only if requested, before anything else)
    total_reads = None
    if tqdm is not None:
        print("Counting alignments...")
        bamfile_count = pysam.AlignmentFile(args.in_bam, "rb", require_index=False)
        total_reads = sum(1 for _ in bamfile_count)
        bamfile_count.close()
    
    # Open input and output
    print("Reading header...")
    bamfile = pysam.AlignmentFile(args.in_bam, "rb", require_index=False)
    outfile = pysam.AlignmentFile(args.out_bam, "wb", header=bamfile.header)
    
    # Check first read for MD tag
    first_read = None
    for read in bamfile:
        first_read = read
        break
    
    if first_read is None:
        print("Error: BAM file appears to be empty.")
        bamfile.close()
        outfile.close()
        sys.exit(1)
    
    try:
        first_read.get_tag('MD')
    except KeyError:
        print("Error: BAM file does not have MD tags.")
        print("You can add them with: samtools calmd -b in.bam reference.fa > out.bam")
        bamfile.close()
        outfile.close()
        sys.exit(1)
    
    # Now start progress bar (after counting and header reading)
    progress_bar = None
    if tqdm is not None:
        print("Processing reads...")
        progress_bar = tqdm(total=total_reads, unit=' reads')
    
    # Process first read
    pmd_score = calculate_pmd(first_read, args.stranded)
    if pmd_score is not None:
        first_read.set_tag('DS', round(pmd_score, 3), 'f')
    outfile.write(first_read)
    if progress_bar:
        progress_bar.update(1)
    
    # Process remaining reads
    reads_processed = 1
    for read in bamfile:
        pmd_score = calculate_pmd(read, args.stranded)
        if pmd_score is not None:
            read.set_tag('DS', round(pmd_score, 3), 'f')
        outfile.write(read)
        reads_processed += 1
        
        if progress_bar:
            progress_bar.update(1)
    
    if progress_bar:
        progress_bar.close()
    
    bamfile.close()
    outfile.close()
    
    print(f"Processed {reads_processed} reads. Output written to {args.out_bam}")


if __name__ == "__main__":
    main()
