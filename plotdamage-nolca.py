#!/usr/bin/env python3

"""
plotdamage-nolca.py - Plot damage from a subs file without requiring a tax ID

Part of bamdam-extras. Please cite the bamdam paper if you use this.
https://github.com/bdesanctis/bamdam

This is designed to work with subs files from compute-nolca, but also works
with regular bamdam subs files if you specify --tax.

Usage:
    python plotdamage-nolca.py --in_subs FILE.subs.txt --outplot damage.png
    python plotdamage-nolca.py --in_subs FILE.subs.txt --outplot damage.png --tax 12345
    python plotdamage-nolca.py --in_subs FILE.subs.txt --outplot damage.png --per-contig
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
        description="Plot damage from a subs file (designed for compute-nolca output)"
    )
    parser.add_argument(
        "--in_subs",
        required=True,
        help="Input subs file from compute-nolca or bamdam compute (required)"
    )
    parser.add_argument(
        "--tax",
        default=None,
        help="Taxonomic node ID or name (optional; defaults to first row in file)"
    )
    parser.add_argument(
        "--outplot",
        default="damage_plot.png",
        help="Filename for output plot, ending in .png or .pdf (default: damage_plot.png)"
    )
    parser.add_argument(
        "--ymax",
        default=0,
        help="Maximum for y axis (optional, default: auto)"
    )
    parser.add_argument(
        "--per-contig",
        action="store_true",
        help="Plot all contigs/nodes from the subs file as separate lines"
    )
    return parser.parse_args()


def calculate_damage_for_plot(items):
    """Parse subs data and calculate C>T, G>A, and other frequencies."""
    
    ctp_5prime = {i: 0 for i in range(1, 16)}
    gap_5prime = {i: 0 for i in range(1, 16)}
    total_c_5prime = {i: 0 for i in range(1, 16)}
    total_g_5prime = {i: 0 for i in range(1, 16)}
    
    ctp_3prime = {i: 0 for i in range(1, 16)}
    total_c_3prime = {i: 0 for i in range(1, 16)}
    gap_3prime = {i: 0 for i in range(1, 16)}
    total_g_3prime = {i: 0 for i in range(1, 16)}
    
    other_5prime = {i: 0 for i in range(1, 16)}
    other_3prime = {i: 0 for i in range(1, 16)}
    total_nondeam_5prime = {i: 0 for i in range(1, 16)}
    total_nondeam_3prime = {i: 0 for i in range(1, 16)}
    
    for item in items:
        mutation, proportion = item.split(":")
        
        # Skip CpG/nonCpG entries if present
        if mutation.startswith("CpG") or mutation.startswith("nonCpG"):
            continue
        
        try:
            from_base, to_base, pos = mutation[0], mutation[1], int(mutation[2:])
        except (ValueError, IndexError):
            continue
        
        count = float(proportion)
        
        # 5' positions (1 to 15)
        if 1 <= pos <= 15:
            if from_base == 'C':
                if to_base == 'T':
                    ctp_5prime[pos] += count
                total_c_5prime[pos] += count
            elif from_base == 'G':
                if to_base == 'A':
                    gap_5prime[pos] += count
                total_g_5prime[pos] += count
            
            if from_base == to_base:
                total_nondeam_5prime[pos] += count
            elif not ((from_base == 'C' and to_base == 'T') or (from_base == 'G' and to_base == 'A')):
                other_5prime[pos] += count
                total_nondeam_5prime[pos] += count
        
        # 3' positions (-15 to -1)
        elif -15 <= pos <= -1:
            abs_pos = abs(pos)
            if from_base == 'C':
                if to_base == 'T':
                    ctp_3prime[abs_pos] += count
                total_c_3prime[abs_pos] += count
            elif from_base == 'G':
                if to_base == 'A':
                    gap_3prime[abs_pos] += count
                total_g_3prime[abs_pos] += count
            
            if from_base == to_base:
                total_nondeam_3prime[abs_pos] += count
            elif not ((from_base == 'C' and to_base == 'T') or (from_base == 'G' and to_base == 'A')):
                other_3prime[abs_pos] += count
                total_nondeam_3prime[abs_pos] += count
    
    prop_ct_5prime = {i: ctp_5prime[i] / total_c_5prime[i] if total_c_5prime[i] > 0 else 0 for i in range(1, 16)}
    prop_ga_5prime = {i: gap_5prime[i] / total_g_5prime[i] if total_g_5prime[i] > 0 else 0 for i in range(1, 16)}
    prop_ct_3prime = {i: ctp_3prime[i] / total_c_3prime[i] if total_c_3prime[i] > 0 else 0 for i in range(1, 16)}
    prop_ga_3prime = {i: gap_3prime[i] / total_g_3prime[i] if total_g_3prime[i] > 0 else 0 for i in range(1, 16)}
    prop_other_5prime = {i: other_5prime[i] / total_nondeam_5prime[i] if total_nondeam_5prime[i] > 0 else 0 for i in range(1, 16)}
    prop_other_3prime = {i: other_3prime[i] / total_nondeam_3prime[i] if total_nondeam_3prime[i] > 0 else 0 for i in range(1, 16)}
    
    return prop_ct_5prime, prop_ga_5prime, prop_ct_3prime, prop_ga_3prime, prop_other_5prime, prop_other_3prime


def make_damage_plot(in_subs, tax, plotfile, ymax, per_contig):
    """Create damage plot from subs file."""
    
    if not matplotlib_imported:
        print("Error: Cannot find matplotlib library for plotting. Try: pip install matplotlib")
        sys.exit(1)
    
    if not os.path.exists(in_subs):
        print(f"Error: File {in_subs} does not exist.")
        sys.exit(1)
    
    with open(in_subs, 'r') as f:
        file_content = f.readlines()
    
    if len(file_content) == 0:
        print("Error: Subs file is empty.")
        sys.exit(1)
    
    # Determine which lines to plot
    lines_to_plot = []
    
    if per_contig:
        # Plot all lines
        for line in file_content:
            line = line.strip()
            if line:
                lines_to_plot.append(line)
    elif tax is not None:
        # Find matching line
        tax_pattern = f"^{tax}\\t" if tax.isdigit() else tax
        matched = [line for line in file_content if re.match(tax_pattern, line)]
        
        if len(matched) == 0:
            print(f"Error: No entry found for tax '{tax}'.")
            sys.exit(1)
        elif len(matched) > 1:
            print(f"Error: Multiple matches for '{tax}'. Use a numeric tax ID.")
            sys.exit(1)
        
        lines_to_plot.append(matched[0])
    else:
        # Default: use first line, but error if multiple lines exist
        non_empty_lines = [line for line in file_content if line.strip()]
        if len(non_empty_lines) > 1:
            print(f"Error: Subs file has {len(non_empty_lines)} entries. Please specify --tax or use --per-contig.")
            sys.exit(1)
        lines_to_plot.append(file_content[0])
    
    # Parse all lines
    positions = list(range(-15, 0)) + list(range(1, 16))
    values_all = {key: [] for key in ['Other', 'CT', 'GA']}
    labels = []
    
    for line in lines_to_plot:
        split_line = line.strip().split('\t')
        if len(split_line) < 3:
            continue
        
        tax_id, tax_name, data_part = split_line[0], split_line[1], split_line[2]
        data_items = data_part.split()
        
        if not data_items:
            continue
        
        ctp, gap_5, ctm, gap_3, other_5, other_3 = calculate_damage_for_plot(data_items)
        
        values = {'Other': [0] * 30, 'CT': [0] * 30, 'GA': [0] * 30}
        for i in range(1, 16):
            values['CT'][positions.index(i)] = ctp[i]
            values['GA'][positions.index(i)] = gap_5[i]
            values['Other'][positions.index(i)] = other_5[i]
            
            values['CT'][positions.index(-i)] = ctm[i]
            values['GA'][positions.index(-i)] = gap_3[i]
            values['Other'][positions.index(-i)] = other_3[i]
        
        for key in values:
            values_all[key].append(values[key])
        
        labels.append(tax_name)
    
    if not any(values_all['CT']):
        print("Error: No valid data found in subs file.")
        sys.exit(1)
    
    # Determine y-axis max
    if ymax == 0 or ymax == "0":
        all_vals = []
        for key in ['Other', 'CT', 'GA']:
            for sublist in values_all[key]:
                all_vals.extend(sublist)
        max_y = min(1.0, max(all_vals) * 1.2) if all_vals else 0.5
    else:
        max_y = float(ymax)
    
    # Create plot
    plt.figure(figsize=(10, 5))
    color_palette = {'Other': '#009E73', 'CT': '#F8766D', 'GA': '#56B4E9'}
    
    # 5' plot
    ax1 = plt.subplot(1, 2, 1)
    handles1 = []
    legend_labels = ['Other', 'C to T', 'G to A']
    
    for key in ['Other', 'CT', 'GA']:
        color = color_palette[key]
        for j, file_values in enumerate(values_all[key]):
            alpha = 1.0 if len(values_all[key]) == 1 else 0.7
            line, = ax1.plot(
                [str(pos) for pos in positions if pos > 0],
                [file_values[i] for i, pos in enumerate(positions) if pos > 0],
                color=color, linewidth=2, alpha=alpha
            )
            if j == 0:
                handles1.append(line)
    
    ax1.set_ylim(0, max_y)
    ax1.set_xlabel('Position', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=14)
    ax1.yaxis.set_ticks_position('left')
    ax1.tick_params(right=False, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    
    # 3' plot
    ax2 = plt.subplot(1, 2, 2)
    handles2 = []
    
    for key in ['Other', 'CT', 'GA']:
        color = color_palette[key]
        for j, file_values in enumerate(values_all[key]):
            alpha = 1.0 if len(values_all[key]) == 1 else 0.7
            line, = ax2.plot(
                [str(pos) for pos in positions if pos < 0],
                [file_values[i] for i, pos in enumerate(positions) if pos < 0],
                color=color, linewidth=2, alpha=alpha
            )
            if j == 0:
                handles2.append(line)
    
    ax2.set_ylim(0, max_y)
    ax2.set_xlabel('Position', fontsize=14)
    ax2.yaxis.set_ticks_position('right')
    ax2.tick_params(left=False, labelsize=12)
    ax2.tick_params(axis='x', labelsize=12)
    
    ax2.legend(handles=handles2, labels=legend_labels, loc='upper left', fontsize=12)
    
    # Title
    if len(labels) == 1:
        title = f'Damage Plot for {labels[0]}'
    elif per_contig:
        title = f'Damage Plot ({len(labels)} contigs)'
    else:
        title = 'Damage Plot'
    
    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save
    file_extension = os.path.splitext(plotfile)[1].lower()
    if file_extension == '.pdf':
        plt.savefig(plotfile, format='pdf')
    else:
        if file_extension != '.png':
            print("Warning: Invalid plot file suffix. Saving as PNG.")
        plt.savefig(plotfile)
    
    plt.close()
    print(f"Wrote plot to {plotfile}")


def main():
    args = parse_args()
    make_damage_plot(args.in_subs, args.tax, args.outplot, args.ymax, args.per_contig)


if __name__ == "__main__":
    main()
