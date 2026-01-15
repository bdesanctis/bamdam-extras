#!/usr/bin/env python3

"""
compute-nolca.py - Compute bamdam statistics without an LCA file

Part of bamdam-extras. Please cite the bamdam paper if you use this.
https://github.com/bdesanctis/bamdam

This replicates bamdam compute functionality but without requiring an LCA file.
All reads are treated as belonging to a single "all" node, or with --per-contig,
each reference ID in the BAM becomes its own node.

Usage:
    python compute-nolca.py --in_bam FILE.bam --out_tsv FILE.tsv --out_subs FILE.subs.txt --stranded ds
    python compute-nolca.py --in_bam FILE.bam --out_tsv FILE.tsv --out_subs FILE.subs.txt --stranded ds --per-contig
"""

import sys
import re
import csv
import math
import argparse
import os
import random
from functools import lru_cache

try:
    import pysam
    pysam_imported = True
except ImportError:
    pysam_imported = False

try:
    import hyperloglog
    hyperloglog_imported = True
except ImportError:
    hyperloglog_imported = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute bamdam statistics without an LCA file"
    )
    parser.add_argument(
        "--in_bam",
        type=str,
        required=True,
        help="Path to the BAM file (required)"
    )
    parser.add_argument(
        "--out_tsv",
        type=str,
        required=True,
        help="Path to the output tsv file (required)"
    )
    parser.add_argument(
        "--out_subs",
        type=str,
        required=True,
        help="Path to the output subs file (required)"
    )
    parser.add_argument(
        "--stranded",
        type=str,
        required=True,
        choices=["ss", "ds"],
        help="Either ss for single stranded or ds for double stranded (required)"
    )
    parser.add_argument(
        "--k",
        type=int,
        default=29,
        help="Value of k for per-node counts of unique k-mers and duplicity (default: 29)"
    )
    parser.add_argument(
        "--mode",
        type=int,
        default=1,
        choices=[1, 2, 3],
        help="Mode to calculate stats. 1: use best alignment (recommended), 2: average over reads, 3: average over alignments (default: 1)"
    )
    parser.add_argument(
        "--per-contig",
        action="store_true",
        help="Treat each reference ID as a separate node (default: not set)"
    )
    parser.add_argument(
        "--mincount",
        type=int,
        default=100,
        help="Minimum read count to keep a contig; only used with --per-contig (default: 100)"
    )
    parser.add_argument(
        "--show_progress",
        action="store_true",
        help="Print a progress bar (default: not set)"
    )
    parser.add_argument(
        "--upto",
        help=argparse.SUPPRESS  # Hidden argument to catch and error
    )
    return parser.parse_args()


# ============== Utility functions (from bamdam utils.py) ==============

@lru_cache(maxsize=1024)
def parse_cigar_cached(cigar):
    return re.findall(r"\d+\D", cigar)

@lru_cache(maxsize=1024)
def parse_md_cached(md):
    return re.compile(r"\d+|\^[A-Za-z]+|[A-Za-z]").findall(md)

@lru_cache(maxsize=128)
def is_reverse_strand(flagsum):
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4
    return len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == "1"

@lru_cache(maxsize=128)
def rev_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-', 's': 's'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

@lru_cache(maxsize=128)
def get_rep_kmer(seq):
    return min(seq, rev_complement(seq))

def calculate_dust(seq):
    """Calculate DUST score for sequence complexity (0-100, higher = less complex)"""
    readlength = len(seq)
    if readlength < 3:
        return 0
    
    w = 64
    k = 3
    firstwindowend = min(readlength, w)
    l = firstwindowend - 2
    maxpossibledust = l * (l - 1) / 2
    
    kmer_counts = {}
    for i in range(firstwindowend - k + 1):
        kmer = seq[i:i + k]
        if not all(base in {'A', 'C', 'T', 'G'} for base in kmer):
            return -1
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    
    currentdust = sum((count * (count - 1)) / 2 for count in kmer_counts.values())
    
    if firstwindowend == readlength:
        return currentdust * (100 / maxpossibledust)
    
    maxdust = currentdust
    for i in range(1, readlength - w + 1):
        oldkmer = seq[(i-1):(i+2)]
        newkmer = seq[(i+w-3):(i+w)]
        if not all(base in {'A', 'C', 'T', 'G'} for base in newkmer):
            return -1
        kmer_counts[oldkmer] -= 1
        if kmer_counts[oldkmer] == 0:
            del kmer_counts[oldkmer]
        kmer_counts[newkmer] = kmer_counts.get(newkmer, 0) + 1
        currentdust = sum((count * (count - 1)) / 2 for count in kmer_counts.values())
        if currentdust > maxdust:
            maxdust = currentdust
    
    return maxdust * 100 / maxpossibledust


# ============== Alignment utilities (from bamdam alignment_utils.py) ==============

def get_mismatches(seq, cigar, md):
    """Parse read, CIGAR, and MD to determine mismatches and positions."""
    cig = parse_cigar_cached(cigar)
    md_list = list(parse_md_cached(md))  # make mutable copy
    
    ref_seq = ''
    read_seq = ''
    query_pos = 0
    read_pos = 0
    
    for x in cig:
        cat = x[-1]
        if cat == 'H':
            continue
        
        bases = int(x[:-1])
        
        if cat == 'S':
            md_extra = "softclip" + str(bases)
            if x == cig[-1]:
                md_list.append(md_extra)
            elif x == cig[0]:
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
            num_chars = len(x) - 1
            rec_pos += num_chars
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


def mismatch_table(read, cigar, md, flagsum):
    """Wrapper for get_mismatches that handles strand orientation and position mirroring."""
    COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
    
    mms, refseq, readseq = get_mismatches(read, cigar, md)
    readlength = len(readseq)
    
    non_dash_indices = [i for i, char in enumerate(readseq) if char != "-"]
    readseq_no_read_dashes = "".join(readseq[i] for i in non_dash_indices)
    refseq_no_read_dashes = "".join(refseq[i] for i in non_dash_indices)
    
    matchs = []
    for i, (read_char, ref_char) in enumerate(zip(readseq_no_read_dashes, refseq_no_read_dashes)):
        if read_char == ref_char and read_char in "ACTG":
            pos = non_dash_indices[i] + 1
            matchs.append([ref_char, read_char, pos])
    
    backwards = is_reverse_strand(flagsum)
    
    if backwards:
        mmsc = []
        matchsc = []
        for entry in mms:
            new_entry = [
                COMPLEMENT.get(entry[0], "N"),
                COMPLEMENT.get(entry[1], "N"),
                readlength - entry[2] + 1,
            ]
            mmsc.append(new_entry)
        for entry in matchs:
            new_entry = [
                COMPLEMENT.get(entry[0], "N"),
                COMPLEMENT.get(entry[1], "N"),
                readlength - entry[2] + 1,
            ]
            matchsc.append(new_entry)
    else:
        mmsc = mms
        matchsc = matchs
    
    # Mirror positions around middle of read
    for entry in mmsc:
        pos = entry[2]
        if pos > readlength / 2:
            entry[2] = -(readlength - pos + 1)
    
    for entry in matchsc:
        pos = entry[2]
        if pos > readlength / 2:
            entry[2] = -(readlength - pos + 1)
    
    return mmsc, matchsc, refseq


# ============== K-mer utilities ==============

def get_hll_info(seq, k):
    """Get representative k-mers for HyperLogLog counting."""
    rep_kmers = []
    total_kmers = 0
    if len(seq) > k:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if all(base in {'A', 'C', 'T', 'G'} for base in kmer):
                rep_kmers.append(get_rep_kmer(kmer))
                total_kmers += 1
    return rep_kmers, total_kmers


# ============== Main compute functions ==============

def gather_subs_and_kmers_per_contig(bamfile_path, kn, stranded, show_progress=False):
    """
    Gather substitution and k-mer metrics per contig (reference ID).
    Each read contributes to exactly one contig (its best alignment by AS score).
    Assumes BAM is read-sorted (all alignments for a read are consecutive).
    """
    print("\nGathering substitution and k-mer metrics per contig...")
    
    node_data = {}
    ref_to_node_id = {}
    next_node_id = 1
    
    bamfile = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
    
    # Check for MD tags on first alignment
    first_alignment = None
    for aln in bamfile:
        first_alignment = aln
        break
    
    if first_alignment is None:
        print("Error: BAM file appears to be empty.")
        sys.exit(1)
    
    try:
        first_alignment.get_tag('MD')
    except KeyError:
        print("Error: BAM file does not have MD tags.")
        print("You can add them with: samtools calmd -b in.bam reference.fa > out.bam")
        sys.exit(1)
    
    # Check for AS tags (needed for mode 1)
    try:
        first_alignment.get_tag('AS')
    except KeyError:
        print("Error: BAM file does not have AS (alignment score) tags.")
        print("These are needed to pick the best alignment per read in --per-contig mode.")
        sys.exit(1)
    
    # Check for PMD scores
    are_pmds_in_bam = True
    try:
        first_alignment.get_tag('DS')
    except KeyError:
        are_pmds_in_bam = False
    
    # Reopen to start from beginning
    bamfile.close()
    bamfile = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
    
    # Progress bar setup
    progress_bar = None
    if show_progress:
        try:
            from tqdm import tqdm
            bamfile_count = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
            total_alignments = sum(1 for _ in bamfile_count)
            bamfile_count.close()
            print("Computing...")
            progress_bar = tqdm(total=total_alignments, unit=' alignments')
        except ImportError:
            print("Warning: tqdm not available, progress bar disabled.")
    
    # Initialize tracking variables
    oldreadname = ""
    best_alignments = []
    best_alignment_score = -1000
    has_multi_alignments = False
    readswithNs = 0
    num_alignments_this_read = 0
    
    for alignment in bamfile:
        readname = alignment.query_name
        
        # When we see a new read name, process the previous read
        if readname != oldreadname and oldreadname != "":
            # Check for multi-alignment reads
            if num_alignments_this_read > 1:
                has_multi_alignments = True
            
            # Pick best alignment (randomly if tied)
            if len(best_alignments) > 1:
                best_alignment = random.choice(best_alignments)
            elif len(best_alignments) == 1:
                best_alignment = best_alignments[0]
            else:
                # No valid alignments for this read (all unmapped?), skip
                oldreadname = readname
                best_alignments = []
                best_alignment_score = -1000
                num_alignments_this_read = 0
                if progress_bar:
                    progress_bar.update(1)
                continue
            
            # Get the reference this read maps to
            ref_name = best_alignment.reference_name
            
            # Assign node ID
            if ref_name not in ref_to_node_id:
                ref_to_node_id[ref_name] = next_node_id
                next_node_id += 1
            
            # Initialize node if needed
            if ref_name not in node_data:
                if are_pmds_in_bam:
                    node_data[ref_name] = {
                        'total_reads': 0, 'pmdsover2': 0, 'pmdsover4': 0,
                        'meanlength': 0, 'total_alignments': 0, 'ani': 0,
                        'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                        'tax_path': ref_name, 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                        'totalkmers': 0, 'unaggregatedreads': 0
                    }
                else:
                    node_data[ref_name] = {
                        'total_reads': 0, 'meanlength': 0, 'total_alignments': 0,
                        'ani': 0, 'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                        'tax_path': ref_name, 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                        'totalkmers': 0, 'unaggregatedreads': 0
                    }
            
            nd = node_data[ref_name]
            
            # Process the best alignment
            seq = best_alignment.query_sequence
            readlength = len(seq)
            cigar = best_alignment.cigarstring
            md = best_alignment.get_tag('MD')
            nm = best_alignment.get_tag('NM')
            flagsum = best_alignment.flag
            
            # Get substitutions
            subs, matches, refseq = mismatch_table(seq, cigar, md, flagsum)
            allsubs = subs + matches
            for sub in allsubs:
                key = tuple(sub)
                nd['subs'][key] = nd['subs'].get(key, 0) + 1
            
            # ANI
            ani_for_read = (readlength - nm) / readlength
            nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
            
            # Ref GC
            refgc = (refseq.count('C') + refseq.count('G')) / len(refseq) if len(refseq) > 0 else 0
            nd['avgrefgc'] = (refgc + nd['avgrefgc'] * nd['total_reads']) / (nd['total_reads'] + 1)
            
            # Mean length
            nd['meanlength'] = ((nd['meanlength'] * nd['total_reads']) + readlength) / (nd['total_reads'] + 1)
            
            # DUST
            dust = calculate_dust(seq)
            if dust != -1:
                nd['avgdust'] = ((nd['avgdust'] * nd['total_reads']) + dust) / (nd['total_reads'] + 1)
            else:
                readswithNs += 1
            
            # Read GC
            gc_content = (seq.count('C') + seq.count('G')) / readlength
            nd['avgreadgc'] = ((nd['avgreadgc'] * nd['total_reads']) + gc_content) / (nd['total_reads'] + 1)
            
            # K-mers
            rep_kmers, total_kmers = get_hll_info(seq, kn)
            for kmer in rep_kmers:
                nd['hll'].add(kmer)
            nd['totalkmers'] += total_kmers
            
            # PMD scores
            if are_pmds_in_bam:
                try:
                    pmd = float(best_alignment.get_tag('DS'))
                except KeyError:
                    pmd = 0
                if pmd > 2:
                    nd['pmdsover2'] += 1
                if pmd > 4:
                    nd['pmdsover4'] += 1
            
            nd['total_reads'] += 1
            nd['total_alignments'] += 1
            nd['unaggregatedreads'] += 1
            
            # Reset for next read
            best_alignments = []
            best_alignment_score = -1000
            num_alignments_this_read = 0
        
        # Process current alignment - track best for this read
        if alignment.reference_name is not None:  # Skip unmapped
            alignment_score = alignment.get_tag('AS')
            
            if not best_alignments:
                best_alignments = [alignment]
                best_alignment_score = alignment_score
            elif alignment_score > best_alignment_score:
                best_alignments = [alignment]
                best_alignment_score = alignment_score
            elif alignment_score == best_alignment_score:
                best_alignments.append(alignment)
        
        num_alignments_this_read += 1
        oldreadname = readname
        
        if progress_bar:
            progress_bar.update(1)
    
    # Process the last read
    if oldreadname != "" and best_alignments:
        if num_alignments_this_read > 1:
            has_multi_alignments = True
        
        if len(best_alignments) > 1:
            best_alignment = random.choice(best_alignments)
        else:
            best_alignment = best_alignments[0]
        
        ref_name = best_alignment.reference_name
        
        if ref_name not in ref_to_node_id:
            ref_to_node_id[ref_name] = next_node_id
            next_node_id += 1
        
        if ref_name not in node_data:
            if are_pmds_in_bam:
                node_data[ref_name] = {
                    'total_reads': 0, 'pmdsover2': 0, 'pmdsover4': 0,
                    'meanlength': 0, 'total_alignments': 0, 'ani': 0,
                    'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                    'tax_path': ref_name, 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                    'totalkmers': 0, 'unaggregatedreads': 0
                }
            else:
                node_data[ref_name] = {
                    'total_reads': 0, 'meanlength': 0, 'total_alignments': 0,
                    'ani': 0, 'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                    'tax_path': ref_name, 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                    'totalkmers': 0, 'unaggregatedreads': 0
                }
        
        nd = node_data[ref_name]
        seq = best_alignment.query_sequence
        readlength = len(seq)
        cigar = best_alignment.cigarstring
        md = best_alignment.get_tag('MD')
        nm = best_alignment.get_tag('NM')
        flagsum = best_alignment.flag
        
        subs, matches, refseq = mismatch_table(seq, cigar, md, flagsum)
        allsubs = subs + matches
        for sub in allsubs:
            key = tuple(sub)
            nd['subs'][key] = nd['subs'].get(key, 0) + 1
        
        ani_for_read = (readlength - nm) / readlength
        nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
        
        refgc = (refseq.count('C') + refseq.count('G')) / len(refseq) if len(refseq) > 0 else 0
        nd['avgrefgc'] = (refgc + nd['avgrefgc'] * nd['total_reads']) / (nd['total_reads'] + 1)
        
        nd['meanlength'] = ((nd['meanlength'] * nd['total_reads']) + readlength) / (nd['total_reads'] + 1)
        
        dust = calculate_dust(seq)
        if dust != -1:
            nd['avgdust'] = ((nd['avgdust'] * nd['total_reads']) + dust) / (nd['total_reads'] + 1)
        else:
            readswithNs += 1
        
        gc_content = (seq.count('C') + seq.count('G')) / readlength
        nd['avgreadgc'] = ((nd['avgreadgc'] * nd['total_reads']) + gc_content) / (nd['total_reads'] + 1)
        
        rep_kmers, total_kmers = get_hll_info(seq, kn)
        for kmer in rep_kmers:
            nd['hll'].add(kmer)
        nd['totalkmers'] += total_kmers
        
        if are_pmds_in_bam:
            try:
                pmd = float(best_alignment.get_tag('DS'))
            except KeyError:
                pmd = 0
            if pmd > 2:
                nd['pmdsover2'] += 1
            if pmd > 4:
                nd['pmdsover4'] += 1
        
        nd['total_reads'] += 1
        nd['total_alignments'] += 1
        nd['unaggregatedreads'] += 1
    
    if progress_bar:
        progress_bar.close()
    
    bamfile.close()
    
    # Warn about multi-alignment reads
   # if has_multi_alignments:
# eh, whatever.    

    if readswithNs > 0:
        print(f"\nSkipped DUST/k-mer calculation for {readswithNs} reads with non-ACGT characters.")
    
    print(f"\nGathered data for {len(node_data)} contig(s). Writing output files...")
    
    return node_data, are_pmds_in_bam, ref_to_node_id


def gather_subs_and_kmers(bamfile_path, kn, mode, stranded, per_contig, show_progress=False):
    """Main function to gather substitution and k-mer metrics."""
    print("\nGathering substitution and k-mer metrics...")
    
    node_data = {}
    ref_to_node_id = {}  # maps ref names to numeric IDs for --per-contig
    next_node_id = 1
    
    bamfile = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
    
    # Check for MD tags on first alignment
    first_alignment = None
    for aln in bamfile:
        first_alignment = aln
        break
    
    if first_alignment is None:
        print("Error: BAM file appears to be empty.")
        sys.exit(1)
    
    try:
        first_alignment.get_tag('MD')
    except KeyError:
        print("Error: BAM file does not have MD tags.")
        print("You can add them with: samtools calmd -b in.bam reference.fa > out.bam")
        sys.exit(1)
    
    # Check for PMD scores
    are_pmds_in_bam = True
    try:
        first_alignment.get_tag('DS')
    except KeyError:
        are_pmds_in_bam = False
    
    # Reopen to start from beginning
    bamfile.close()
    bamfile = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
    
    # Progress bar setup
    progress_bar = None
    if show_progress:
        try:
            from tqdm import tqdm
            bamfile_count = pysam.AlignmentFile(bamfile_path, "rb", require_index=False)
            total_alignments = sum(1 for _ in bamfile_count)
            bamfile_count.close()
            print("Computing...")
            progress_bar = tqdm(total=total_alignments, unit=' alignments')
        except ImportError:
            print("Warning: tqdm not available, progress bar disabled.")
    
    # Initialize tracking variables
    oldreadname = ""
    oldmd = ""
    oldcigar = ""
    oldflagsum = ""
    num_alignments = 0
    currentsubdict = {}
    nms = 0
    pmdsover2 = 0
    pmdsover4 = 0
    readswithNs = 0
    best_alignments = []
    best_alignment_score = -1000
    refgc = 0
    current_ref_name = None  # For --per-contig mode
    
    alignment_count = 0
    
    for alignment in bamfile:
        readname = alignment.query_name
        
        # When we encounter a new read, process the previous read
        if readname != oldreadname and oldreadname != "":
            # Determine which node(s) to update
            if per_contig:
                # Use the reference from the best alignment (mode 1) or first alignment
                if mode == 1 and best_alignments:
                    ref_name = best_alignments[0].reference_name
                else:
                    ref_name = current_ref_name
                
                if ref_name not in ref_to_node_id:
                    ref_to_node_id[ref_name] = next_node_id
                    next_node_id += 1
                node_name = ref_name
            else:
                node_name = "all"
            
            # Process based on mode
            if mode == 1:
                if len(best_alignments) > 1:
                    best_alignment = random.choice(best_alignments)
                elif len(best_alignments) == 1:
                    best_alignment = best_alignments[0]
                else:
                    print("Error: No best alignment found for read.")
                    sys.exit(1)
                
                seq = best_alignment.query_sequence
                readlength = len(seq)
                cigar = best_alignment.cigarstring
                md = best_alignment.get_tag('MD')
                nms = best_alignment.get_tag('NM')
                
                if are_pmds_in_bam:
                    try:
                        pmd = float(best_alignment.get_tag('DS'))
                    except KeyError:
                        pmd = 0
                    if pmd > 2:
                        pmdsover2 = 1
                    if pmd > 4:
                        pmdsover4 = 1
                
                flagsum = best_alignment.flag
                subs, matches, reference_seq = mismatch_table(seq, cigar, md, flagsum)
                refgc = (reference_seq.count('C') + reference_seq.count('G')) / len(reference_seq) if len(reference_seq) > 0 else 0
                
                allsubs = subs + matches
                for sub in allsubs:
                    key = tuple(sub)
                    currentsubdict[key] = currentsubdict.get(key, 0) + 1
            else:
                # Mode 2 or 3: seq was already set during alignment iteration
                seq = oldseq
                readlength = len(seq)
            
            # Calculate DUST and k-mer info
            dust = calculate_dust(seq)
            rep_kmers, total_kmers = get_hll_info(seq, kn)
            
            # Initialize node if needed
            if node_name not in node_data:
                if are_pmds_in_bam:
                    node_data[node_name] = {
                        'total_reads': 0, 'pmdsover2': 0, 'pmdsover4': 0,
                        'meanlength': 0, 'total_alignments': 0, 'ani': 0,
                        'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                        'tax_path': "", 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                        'totalkmers': 0, 'unaggregatedreads': 0
                    }
                else:
                    node_data[node_name] = {
                        'total_reads': 0, 'meanlength': 0, 'total_alignments': 0,
                        'ani': 0, 'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                        'tax_path': "", 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                        'totalkmers': 0, 'unaggregatedreads': 0
                    }
            
            nd = node_data[node_name]
            
            # Update node data
            nd['meanlength'] = ((nd['meanlength'] * nd['total_reads']) + readlength) / (nd['total_reads'] + 1)
            
            if dust != -1:
                nd['avgdust'] = ((nd['avgdust'] * nd['total_reads']) + dust) / (nd['total_reads'] + 1)
            else:
                readswithNs += 1
            
            # ANI calculation
            if mode == 1:
                ani_for_read = (readlength - nms) / readlength
                nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
            elif mode == 2:
                ani_for_read = (readlength - (nms / num_alignments)) / readlength
                nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
            elif mode == 3:
                ani_for_read = (readlength * num_alignments - nms) / readlength
                nd['ani'] = (ani_for_read + nd['ani'] * nd['total_alignments']) / (nd['total_alignments'] + num_alignments)
            
            # GC content
            gc_content = (seq.count('C') + seq.count('G')) / readlength
            nd['avgreadgc'] = ((nd['avgreadgc'] * nd['total_reads']) + gc_content) / (nd['total_reads'] + 1)
            
            if mode == 1:
                nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_reads']) + refgc) / (nd['total_reads'] + 1)
            elif mode == 2:
                nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_reads']) + (refgc / num_alignments)) / (nd['total_reads'] + 1)
            elif mode == 3:
                nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_alignments']) + refgc) / (nd['total_alignments'] + num_alignments)
            
            nd['total_alignments'] += num_alignments
            
            # PMD scores
            if are_pmds_in_bam:
                if mode == 2:
                    nd['pmdsover2'] += pmdsover2 / num_alignments
                    nd['pmdsover4'] += pmdsover4 / num_alignments
                else:
                    nd['pmdsover2'] += pmdsover2
                    nd['pmdsover4'] += pmdsover4
            
            # K-mer updates
            for kmer in rep_kmers:
                nd['hll'].add(kmer)
            nd['totalkmers'] += total_kmers
            
            # Substitution updates
            if currentsubdict:
                for sub, count in currentsubdict.items():
                    if mode == 1 or mode == 3:
                        nd['subs'][sub] = nd['subs'].get(sub, 0) + count
                    elif mode == 2:
                        nd['subs'][sub] = nd['subs'].get(sub, 0) + count / num_alignments
            
            # Set tax_path
            if nd['tax_path'] == "":
                nd['tax_path'] = node_name
            
            nd['total_reads'] += 1
            nd['unaggregatedreads'] += 1
            
            # Reset for next read
            oldreadname = readname
            oldmd = ""
            oldcigar = ""
            oldflagsum = ""
            currentsubdict = {}
            num_alignments = 0
            best_alignments = []
            best_alignment_score = -1000
            nms = 0
            refgc = 0
            pmdsover2 = 0
            pmdsover4 = 0
        
        # Process current alignment
        num_alignments += 1
        current_ref_name = alignment.reference_name
        
        if mode == 1:
            try:
                alignment_score = alignment.get_tag("AS")
            except KeyError:
                print("Error: No AS (alignment score) tags in BAM. Required for mode 1.")
                print("Use mode 2 or 3, or add AS tags to your BAM.")
                sys.exit(1)
            
            if not best_alignments:
                best_alignments = [alignment]
                best_alignment_score = alignment_score
            elif alignment_score > best_alignment_score:
                best_alignments = [alignment]
                best_alignment_score = alignment_score
            elif alignment_score == best_alignment_score:
                best_alignments.append(alignment)
        
        elif mode == 2 or mode == 3:
            seq = alignment.query_sequence
            oldseq = seq  # Save for later
            readlength = len(seq)
            cigar = alignment.cigarstring
            md = alignment.get_tag('MD')
            nms += alignment.get_tag('NM')
            
            if are_pmds_in_bam:
                try:
                    pmd = float(alignment.get_tag('DS'))
                except KeyError:
                    pmd = 0
                if pmd > 2:
                    pmdsover2 += 1
                if pmd > 4:
                    pmdsover4 += 1
            
            flagsum = alignment.flag
            
            if (readname != oldreadname) or (cigar != oldcigar) or (md != oldmd) or (flagsum != oldflagsum):
                subs, matches, refseq = mismatch_table(seq, cigar, md, flagsum)
                oldcigar = cigar
                oldmd = md
                oldflagsum = flagsum
            
            refgc += (refseq.count('C') + refseq.count('G')) / len(refseq) if len(refseq) > 0 else 0
            
            allsubs = subs + matches
            for sub in allsubs:
                key = tuple(sub)
                currentsubdict[key] = currentsubdict.get(key, 0) + 1
        
        if oldreadname == "":
            oldreadname = readname
            if mode == 2 or mode == 3:
                oldseq = alignment.query_sequence
        
        alignment_count += 1
        if progress_bar:
            progress_bar.update(1)
    
    # Process the last read
    if oldreadname != "":
        if per_contig:
            if mode == 1 and best_alignments:
                ref_name = best_alignments[0].reference_name
            else:
                ref_name = current_ref_name
            if ref_name not in ref_to_node_id:
                ref_to_node_id[ref_name] = next_node_id
            node_name = ref_name
        else:
            node_name = "all"
        
        if mode == 1:
            if best_alignments:
                best_alignment = random.choice(best_alignments) if len(best_alignments) > 1 else best_alignments[0]
                seq = best_alignment.query_sequence
                readlength = len(seq)
                cigar = best_alignment.cigarstring
                md = best_alignment.get_tag('MD')
                nms = best_alignment.get_tag('NM')
                
                if are_pmds_in_bam:
                    try:
                        pmd = float(best_alignment.get_tag('DS'))
                    except KeyError:
                        pmd = 0
                    pmdsover2 = 1 if pmd > 2 else 0
                    pmdsover4 = 1 if pmd > 4 else 0
                
                flagsum = best_alignment.flag
                subs, matches, reference_seq = mismatch_table(seq, cigar, md, flagsum)
                refgc = (reference_seq.count('C') + reference_seq.count('G')) / len(reference_seq) if len(reference_seq) > 0 else 0
                
                allsubs = subs + matches
                for sub in allsubs:
                    key = tuple(sub)
                    currentsubdict[key] = currentsubdict.get(key, 0) + 1
        else:
            seq = oldseq
            readlength = len(seq)
        
        dust = calculate_dust(seq)
        rep_kmers, total_kmers = get_hll_info(seq, kn)
        
        if node_name not in node_data:
            if are_pmds_in_bam:
                node_data[node_name] = {
                    'total_reads': 0, 'pmdsover2': 0, 'pmdsover4': 0,
                    'meanlength': 0, 'total_alignments': 0, 'ani': 0,
                    'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                    'tax_path': "", 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                    'totalkmers': 0, 'unaggregatedreads': 0
                }
            else:
                node_data[node_name] = {
                    'total_reads': 0, 'meanlength': 0, 'total_alignments': 0,
                    'ani': 0, 'avgdust': 0, 'avgreadgc': 0, 'avgrefgc': 0,
                    'tax_path': "", 'subs': {}, 'hll': hyperloglog.HyperLogLog(0.01),
                    'totalkmers': 0, 'unaggregatedreads': 0
                }
        
        nd = node_data[node_name]
        nd['meanlength'] = ((nd['meanlength'] * nd['total_reads']) + readlength) / (nd['total_reads'] + 1)
        
        if dust != -1:
            nd['avgdust'] = ((nd['avgdust'] * nd['total_reads']) + dust) / (nd['total_reads'] + 1)
        
        if mode == 1:
            ani_for_read = (readlength - nms) / readlength
            nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
        elif mode == 2:
            ani_for_read = (readlength - (nms / num_alignments)) / readlength if num_alignments > 0 else 0
            nd['ani'] = (ani_for_read + nd['ani'] * nd['total_reads']) / (nd['total_reads'] + 1)
        elif mode == 3:
            ani_for_read = (readlength * num_alignments - nms) / readlength if num_alignments > 0 else 0
            nd['ani'] = (ani_for_read + nd['ani'] * nd['total_alignments']) / (nd['total_alignments'] + num_alignments) if (nd['total_alignments'] + num_alignments) > 0 else 0
        
        gc_content = (seq.count('C') + seq.count('G')) / readlength
        nd['avgreadgc'] = ((nd['avgreadgc'] * nd['total_reads']) + gc_content) / (nd['total_reads'] + 1)
        
        if mode == 1:
            nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_reads']) + refgc) / (nd['total_reads'] + 1)
        elif mode == 2 and num_alignments > 0:
            nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_reads']) + (refgc / num_alignments)) / (nd['total_reads'] + 1)
        elif mode == 3:
            nd['avgrefgc'] = ((nd['avgrefgc'] * nd['total_alignments']) + refgc) / (nd['total_alignments'] + num_alignments) if (nd['total_alignments'] + num_alignments) > 0 else 0
        
        nd['total_alignments'] += num_alignments
        
        if are_pmds_in_bam:
            if mode == 2 and num_alignments > 0:
                nd['pmdsover2'] += pmdsover2 / num_alignments
                nd['pmdsover4'] += pmdsover4 / num_alignments
            else:
                nd['pmdsover2'] += pmdsover2
                nd['pmdsover4'] += pmdsover4
        
        for kmer in rep_kmers:
            nd['hll'].add(kmer)
        nd['totalkmers'] += total_kmers
        
        if currentsubdict:
            for sub, count in currentsubdict.items():
                if mode == 1 or mode == 3:
                    nd['subs'][sub] = nd['subs'].get(sub, 0) + count
                elif mode == 2 and num_alignments > 0:
                    nd['subs'][sub] = nd['subs'].get(sub, 0) + count / num_alignments
        
        if nd['tax_path'] == "":
            nd['tax_path'] = node_name
        
        nd['total_reads'] += 1
        nd['unaggregatedreads'] += 1
    
    if progress_bar:
        progress_bar.close()
    
    bamfile.close()
    
    if readswithNs > 0:
        print(f"\nSkipped DUST/k-mer calculation for {readswithNs} reads with non-ACGT characters.")
    
    print(f"\nGathered data for {len(node_data)} node(s). Writing output files...")
    
    return node_data, are_pmds_in_bam, ref_to_node_id


def format_subs(subs, nreads):
    """Format substitution dictionary for output."""
    formatted_subs = []
    
    for key, value in subs.items():
        parts = [str(key[0]), str(key[1]), str(key[2])]
        pos = int(parts[2])
        if (-15 <= pos <= 15) and (parts[0] in {'A', 'C', 'T', 'G'}) and (parts[1] in {'A', 'C', 'T', 'G'}):
            formatted_key = "".join(parts)
            formatted_value = round(value / nreads, 3)
            formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))
    
    formatted_subs.sort(key=lambda x: (x[0] > 0, x[0]))
    return " ".join(sub[1] for sub in formatted_subs)


def calculate_node_damage(subs, stranded):
    """Calculate damage statistics from substitution dictionary."""
    ctp1 = 0  # C>T at 5' position 1
    ctm1 = 0  # C>T at 3' position -1 (ss)
    gam1 = 0  # G>A at 3' position -1 (ds)
    c_p1 = 0  # total C at 5' position 1
    c_m1 = 0  # total C at 3' position -1 (ss)
    g_m1 = 0  # total G at 3' position -1 (ds)
    
    for key, count in subs.items():
        from_base, to_base, pos = str(key[0]), str(key[1]), int(key[2])
        
        if from_base == 'C' and pos == 1:
            if to_base == 'T':
                ctp1 += count
            c_p1 += count
        
        if stranded == "ss" and from_base == 'C' and pos == -1:
            if to_base == 'T':
                ctm1 += count
            c_m1 += count
        
        if stranded == "ds" and from_base == 'G' and pos == -1:
            if to_base == 'A':
                gam1 += count
            g_m1 += count
    
    dp1 = ctp1 / c_p1 if c_p1 > 0 else 0
    dm1 = (ctm1 / c_m1 if c_m1 > 0 else 0) if stranded == "ss" else (gam1 / g_m1 if g_m1 > 0 else 0)
    
    return dp1, dm1


def write_output(node_data, tsv_path, subs_path, stranded, pmds_in_bam, mode, ref_to_node_id, per_contig):
    """Write TSV and subs output files."""
    statsfile = open(tsv_path, 'w', newline='')
    subsfile = open(subs_path, 'w', newline='')
    
    # Header
    if pmds_in_bam:
        header = ['TaxNodeID', 'TaxName', 'TotalReads', 'Duplicity', 'MeanDust',
                  'Damage+1', 'Damage-1', 'MeanLength', 'ANI', 'AvgReadGC', 'AvgRefGC',
                  'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments',
                  'PMDsover2', 'PMDSover4', 'UnaggregatedReads', 'taxpath']
    else:
        header = ['TaxNodeID', 'TaxName', 'TotalReads', 'Duplicity', 'MeanDust',
                  'Damage+1', 'Damage-1', 'MeanLength', 'ANI', 'AvgReadGC', 'AvgRefGC',
                  'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments',
                  'UnaggregatedReads', 'taxpath']
    
    statsfile.write("\t".join(header) + "\n")
    
    writer = csv.writer(statsfile, delimiter="\t", quotechar='"',
                        quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
    subswriter = csv.writer(subsfile, delimiter="\t", quotechar='"',
                            quoting=csv.QUOTE_NONE, lineterminator="\n")
    
    rows = []
    subsrows = {}
    
    for node_name in node_data:
        nd = node_data[node_name]
        
        # Get node ID
        if per_contig:
            node_id = ref_to_node_id.get(node_name, 1)
        else:
            node_id = 1
        
        # Format subs
        if mode == 1 or mode == 2:
            fsubs = format_subs(nd['subs'], nd['total_reads'])
        else:  # mode 3
            fsubs = format_subs(nd['subs'], nd['total_alignments'])
        
        # Duplicity
        numuniquekmers = len(nd['hll'])
        if numuniquekmers > 0:
            duplicity = nd['totalkmers'] / numuniquekmers
            uniquekmersperread = numuniquekmers / nd['total_reads']
        else:
            duplicity = 0
            uniquekmersperread = 0
        
        if uniquekmersperread < 0:
            uniquekmersperread = 0
        
        # Damage
        dp1, dm1 = calculate_node_damage(nd['subs'], stranded)
        
        # Build row
        if pmds_in_bam:
            if mode == 1 or mode == 2:
                pd2 = nd['pmdsover2'] / nd['total_reads'] if nd['total_reads'] > 0 else 0
                pd4 = nd['pmdsover4'] / nd['total_reads'] if nd['total_reads'] > 0 else 0
            else:
                pd2 = nd['pmdsover2'] / nd['total_alignments'] if nd['total_alignments'] > 0 else 0
                pd4 = nd['pmdsover4'] / nd['total_alignments'] if nd['total_alignments'] > 0 else 0
            
            row = [node_id, node_name, nd['total_reads'], round(duplicity, 2),
                   round(nd['avgdust'], 2), round(dp1, 4), round(dm1, 4),
                   round(nd['meanlength'], 2), round(nd['ani'], 4),
                   round(nd['avgreadgc'], 3), round(nd['avgrefgc'], 3),
                   numuniquekmers, round(uniquekmersperread, 3), nd['total_alignments'],
                   round(pd2, 3), round(pd4, 3), nd['unaggregatedreads'], nd['tax_path']]
        else:
            row = [node_id, node_name, nd['total_reads'], round(duplicity, 2),
                   round(nd['avgdust'], 2), round(dp1, 4), round(dm1, 4),
                   round(nd['meanlength'], 2), round(nd['ani'], 4),
                   round(nd['avgreadgc'], 3), round(nd['avgrefgc'], 3),
                   numuniquekmers, round(uniquekmersperread, 3), nd['total_alignments'],
                   nd['unaggregatedreads'], nd['tax_path']]
        
        rows.append(row)
        subsrows[node_id] = [node_id, node_name, fsubs]
    
    # Sort by read count
    rows.sort(key=lambda x: x[2], reverse=True)
    
    for row in rows:
        writer.writerow(row)
    
    for row in rows:
        subswriter.writerow(subsrows[row[0]])
    
    statsfile.close()
    subsfile.close()
    
    print(f"Wrote {tsv_path} and {subs_path}")


def main():
    args = parse_args()
    
    # Check for --upto flag
    if args.upto is not None:
        print("Error: --upto is not supported in compute-nolca (there is no taxonomy).")
        sys.exit(1)
    
    if not pysam_imported:
        print("Error: pysam is required. Install with: pip install pysam")
        sys.exit(1)
    
    if not hyperloglog_imported:
        print("Error: hyperloglog is required. Install with: pip install hyperloglog")
        sys.exit(1)
    
    if not os.path.exists(args.in_bam):
        print(f"Error: BAM file {args.in_bam} does not exist.")
        sys.exit(1)
    
    # Handle --per-contig mode logic
    if args.per_contig:
        if args.mode != 1:
            print("Error: --per-contig only supports mode 1 (best alignment per read).")
            print("Modes 2 and 3 don't make sense here because reads could contribute")
            print("to multiple contigs' stats. Use mode 1 (the default) or pre-filter")
            print("your BAM to one alignment per read.")
            sys.exit(1)
        
        node_data, pmds_in_bam, ref_to_node_id = gather_subs_and_kmers_per_contig(
            args.in_bam, args.k, args.stranded, args.show_progress
        )
        mode = 1
        
        # Filter by mincount
        if args.mincount > 0:
            original_count = len(node_data)
            node_data = {k: v for k, v in node_data.items() if v['total_reads'] >= args.mincount}
            filtered_count = original_count - len(node_data)
            if filtered_count > 0:
                print(f"Filtered out {filtered_count} contigs with fewer than {args.mincount} reads (--mincount).")
    else:
        mode = args.mode
        node_data, pmds_in_bam, ref_to_node_id = gather_subs_and_kmers(
            args.in_bam, args.k, mode, args.stranded, False, args.show_progress
        )
    
    write_output(node_data, args.out_tsv, args.out_subs, args.stranded,
                 pmds_in_bam, mode, ref_to_node_id, args.per_contig)
    
    print("Done!")


if __name__ == "__main__":
    main()
