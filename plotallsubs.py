#!/usr/bin/env python3

"""
plotallsubs.py - Plot all substitution types from a bamdam subs file

Part of bamdam-extras. Please cite the bamdam paper if you use this.
https://github.com/bdesanctis/bamdam

Usage:
    python plotallsubs.py --in_subs FILE --tax TAXID --outplot OUTPUT.png
"""

import sys
import re
import argparse
import os

try:
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except ImportError:
    matplotlib_imported = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot all substitution types from a bamdam subs file"
    )
    parser.add_argument(
        "--in_subs",
        required=True,
        help="Input subs file from bamdam compute (required)"
    )
    parser.add_argument(
        "--tax",
        required=True,
        help="Taxonomic node ID (required)"
    )
    parser.add_argument(
        "--outplot",
        default="allsubs_plot.png",
        help="Filename for output plot, ending in .png or .pdf (default: allsubs_plot.png)"
    )
    parser.add_argument(
        "--ymax",
        default=0,
        help="Maximum for y axis (optional, default: auto)"
    )
    return parser.parse_args()


def calculate_all_subs(items):
    """
    Parse subs data and calculate frequencies for all 12 substitution types.
    Returns dict of {sub_type: {position: frequency}} for both 5' and 3' ends.
    """
    bases = ['A', 'C', 'G', 'T']
    
    # All 12 substitution types (excluding self-matches)
    sub_types = []
    for from_base in bases:
        for to_base in bases:
            if from_base != to_base:
                sub_types.append(f"{from_base}>{to_base}")
    
    # Initialize counters for 5' (positions 1-15) and 3' (positions -1 to -15)
    # For each position, track counts of each substitution and total of each "from" base
    counts_5prime = {sub: {i: 0 for i in range(1, 16)} for sub in sub_types}
    counts_3prime = {sub: {i: 0 for i in range(1, 16)} for sub in sub_types}
    totals_5prime = {base: {i: 0 for i in range(1, 16)} for base in bases}
    totals_3prime = {base: {i: 0 for i in range(1, 16)} for base in bases}
    
    for item in items:
        mutation, proportion = item.split(":")
        
        # Skip CpG/nonCpG entries if present (from --udg mode)
        if mutation.startswith("CpG") or mutation.startswith("nonCpG"):
            continue
            
        try:
            from_base, to_base, pos = mutation[0], mutation[1], int(mutation[2:])
        except (ValueError, IndexError):
            print(f"Warning: Could not parse mutation '{mutation}', skipping.")
            continue
        
        if from_base not in bases or to_base not in bases:
            continue
            
        count = float(proportion)
        
        # 5' positions (1 to 15)
        if 1 <= pos <= 15:
            totals_5prime[from_base][pos] += count
            if from_base != to_base:
                sub_type = f"{from_base}>{to_base}"
                counts_5prime[sub_type][pos] += count
        
        # 3' positions (-15 to -1)
        elif -15 <= pos <= -1:
            abs_pos = abs(pos)
            totals_3prime[from_base][abs_pos] += count
            if from_base != to_base:
                sub_type = f"{from_base}>{to_base}"
                counts_3prime[sub_type][abs_pos] += count
    
    # Calculate proportions
    prop_5prime = {}
    prop_3prime = {}
    
    for sub in sub_types:
        from_base = sub[0]
        prop_5prime[sub] = {
            i: counts_5prime[sub][i] / totals_5prime[from_base][i] 
            if totals_5prime[from_base][i] > 0 else 0 
            for i in range(1, 16)
        }
        prop_3prime[sub] = {
            i: counts_3prime[sub][i] / totals_3prime[from_base][i] 
            if totals_3prime[from_base][i] > 0 else 0 
            for i in range(1, 16)
        }
    
    return prop_5prime, prop_3prime, sub_types


def make_allsubs_plot(in_subs, tax, plotfile, ymax=0):
    """
    Create a damage-style plot showing all 12 substitution types.
    """
    if not matplotlib_imported:
        print("Error: Cannot find matplotlib library for plotting. Try: pip install matplotlib")
        sys.exit(1)
    
    if not os.path.exists(in_subs):
        print(f"Error: File {in_subs} does not exist.")
        sys.exit(1)
    
    with open(in_subs, 'r') as f:
        file_content = f.readlines()
    
    # Match tax ID at start of line followed by tab
    tax_pattern = f"^{tax}\\t" if tax.isdigit() else tax
    matched_line = [line for line in file_content if re.match(tax_pattern, line)]
    
    if len(matched_line) == 0:
        print(f"Error: {in_subs} does not contain an entry for {tax}.")
        sys.exit(1)
    elif len(matched_line) > 1:
        print(f"Error: More than one matching line found. Please use a numeric tax ID.")
        sys.exit(1)
    
    split_line = matched_line[0].split('\t')
    tax_id, tax_name, data_part = split_line[0], split_line[1], split_line[2]
    data_items = data_part.split()
    
    prop_5prime, prop_3prime, sub_types = calculate_all_subs(data_items)
    
    # Color palette - keep C>T and G>A as in original, add colors for the rest
    # Original: C>T = '#F8766D' (red), G>A = '#56B4E9' (blue)
    color_palette = {
        'C>T': '#F8766D',  # red (original)
        'G>A': '#56B4E9',  # blue (original)
        'C>A': '#E69F00',  # orange
        'C>G': '#009E73',  # green
        'T>A': '#CC79A7',  # pink
        'T>C': '#0072B2',  # dark blue
        'T>G': '#D55E00',  # vermillion
        'A>C': '#999999',  # gray
        'A>G': '#882255',  # purple
        'A>T': '#44AA99',  # teal
        'G>C': '#DDCC77',  # sand
        'G>T': '#117733',  # dark green
    }
    
    positions_5prime = list(range(1, 16))
    positions_3prime = list(range(-15, 0))
    
    # Determine y-axis max
    if ymax == 0 or ymax == "0":
        all_values = []
        for sub in sub_types:
            all_values.extend(prop_5prime[sub].values())
            all_values.extend(prop_3prime[sub].values())
        max_y = min(1.0, max(all_values) * 1.2) if all_values else 0.5
    else:
        max_y = float(ymax)
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 5' plot (left)
    for sub in sub_types:
        values = [prop_5prime[sub][i] for i in positions_5prime]
        ax1.plot(
            [str(p) for p in positions_5prime],
            values,
            label=sub,
            color=color_palette[sub],
            linewidth=2
        )
    
    ax1.set_ylim(0, max_y)
    ax1.set_xlabel('Position', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=14)
    ax1.yaxis.set_ticks_position('left')
    ax1.tick_params(right=False, labelsize=12)
    ax1.tick_params(axis='x', labelsize=10)
    
    # 3' plot (right)
    for sub in sub_types:
        values = [prop_3prime[sub][abs(i)] for i in positions_3prime]
        ax2.plot(
            [str(p) for p in positions_3prime],
            values,
            label=sub,
            color=color_palette[sub],
            linewidth=2
        )
    
    ax2.set_ylim(0, max_y)
    ax2.set_xlabel('Position', fontsize=14)
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_visible(False)  # hide right y-axis, redundant with left
    ax2.tick_params(left=False, labelsize=12)
    ax2.tick_params(axis='x', labelsize=10)
    
    # Legend on the right side, outside the plot
    ax2.legend(
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
        frameon=True
    )
    
    plt.suptitle(f'All Substitutions for {tax_name} (tax ID {tax_id})', fontsize=16)
    plt.tight_layout(rect=[0, 0, 0.85, 0.95])  # leave room for legend
    
    # Save
    file_extension = os.path.splitext(plotfile)[1].lower()
    if file_extension == '.pdf':
        plt.savefig(plotfile, format='pdf', bbox_inches='tight')
    else:
        if file_extension != '.png':
            print("Warning: Invalid plot file suffix. Saving as PNG.")
        plt.savefig(plotfile, bbox_inches='tight')
    
    plt.close()
    print(f"Wrote plot to {plotfile}")


def main():
    args = parse_args()
    make_allsubs_plot(args.in_subs, args.tax, args.outplot, args.ymax)


if __name__ == "__main__":
    main()
