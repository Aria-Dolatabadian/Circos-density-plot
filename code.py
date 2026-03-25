import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import warnings
warnings.filterwarnings('ignore')
# ============================================================
# READ DATA FROM XLSX
# ============================================================
xlsx_file = sys.argv[1] if len(sys.argv) > 1 else 'circos_data.xlsx'
print(f"Reading data from: {xlsx_file}")
chr_info = pd.read_excel(xlsx_file, sheet_name='Chromosomes')
heat_data = pd.read_excel(xlsx_file, sheet_name='Heatmap_Data')
print(f"  - {len(chr_info)} chromosomes")
print(f"  - {len(heat_data)} heatmap data points")
tracks = ['A', 'B', 'C', 'D', 'E']
# ============================================================
# CREATE CIRCOS PLOT
# ============================================================
# --- Chromosome colors (for the outer chromosome ring) ---
chr_colors = {
    'Chr1':  '#DC143C',  'Chr2':  '#4169E1',  'Chr3':  '#2E8B57',
    'Chr4':  '#8B008B',  'Chr5':  '#FF8C00',  'Chr6':  '#FFD700',
    'Chr7':  '#8B4513',  'Chr8':  '#FF69B4',  'Chr9':  '#A9A9A9',
    'Chr10': '#000080',  'Chr11': '#FF00FF',  'Chr12': '#32CD32',
    'Chr13': '#00CED1',
}
# --- Genotype colors: each track (A-E) gets one distinct color ---
genotype_colors = {
    'A': '#DC143C',  # Red
    'B': '#4169E1',  # Blue
    'C': '#2E8B57',  # Green
    'D': '#FF8C00',  # Orange
    'E': '#8B008B',  # Purple
}
# --- Heatmap colormap ---
cmap_color_list = [
    '#CC0000', '#E06600', '#FF9900', '#FFCC66', '#FFFFBF',
    '#D5E8A0', '#8DD3C7', '#66B2CC', '#3366CC',
]
cmap = LinearSegmentedColormap.from_list('circos_hm', cmap_color_list, N=256)
norm = Normalize(vmin=5.0, vmax=40.0)
# --- Angular layout ---
gap_deg = 2.0
n_chr = len(chr_info)
total_gap = gap_deg * n_chr
total_data_span = 360.0 - total_gap
total_genome = chr_info['Length_Mb'].sum()
chr_layout = {}
cur_deg = 0.0
for _, row in chr_info.iterrows():
    chrom = row['Chromosome']
    length = row['Length_Mb']
    span = (length / total_genome) * total_data_span
    chr_layout[chrom] = {
        'start': cur_deg, 'span': span,
        'end': cur_deg + span, 'center': cur_deg + span / 2.0,
        'length': length,
    }
    cur_deg += span + gap_deg
# --- Radial layout ---
chr_ring_out = 0.93
chr_ring_h = 0.035
chr_ring_in = chr_ring_out - chr_ring_h
# Ticks are OUTSIDE the chromosome ring
tick_h = 0.022
tick_bot_r = chr_ring_out          # ticks start at outer edge of chr ring
tick_top_r = tick_bot_r + tick_h   # ticks extend outward
heatmap_h = 0.040
bar_h = 0.025
inter_track_gap = 0.005
track_layout = {}
r = chr_ring_in - 0.004  # tracks start just inside the chromosome ring
for t in tracks:
    hm_out = r
    hm_in = r - heatmap_h
    bar_out = hm_in
    bar_in = bar_out - bar_h
    track_layout[t] = {
        'hm_out': hm_out, 'hm_in': hm_in,
        'bar_out': bar_out, 'bar_in': bar_in,
    }
    r = bar_in - inter_track_gap
def deg2rad(d):
    return np.radians(d)
def tangent_rotation(theta_deg):
    rot = 180.0 - theta_deg
    rot = (rot + 180.0) % 360.0 - 180.0
    if rot > 90.0 or rot < -90.0:
        rot += 180.0
        rot = (rot + 180.0) % 360.0 - 180.0
    return rot
# --- Create figure: wider to make room for side legends ---
fig = plt.figure(figsize=(20, 16), facecolor='white')
# Main circos plot in the center-left portion
ax = fig.add_axes([0.05, 0.05, 0.65, 0.90], projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_ylim(0, 1.15)
ax.axis('off')
print("Drawing circos plot...")
# --- Draw each chromosome ---
for chrom_name, layout in chr_layout.items():
    start_rad = deg2rad(layout['start'])
    span_rad = deg2rad(layout['span'])
    center_rad = deg2rad(layout['center'])
    chr_len = layout['length']
    chr_color = chr_colors.get(chrom_name, '#888888')
    # Outer chromosome ring
    ax.bar(center_rad, height=chr_ring_h, width=span_rad,
           bottom=chr_ring_in, color=chr_color, edgecolor='none', linewidth=0)
    # Heatmap data
    mask = heat_data['Chromosome'] == chrom_name
    chr_hm = heat_data.loc[mask].sort_values('Position_Mb')
    n_pts = len(chr_hm)
    # Tracks A-E
    for t_name in tracks:
        tl = track_layout[t_name]
        geno_color = genotype_colors[t_name]
        ax.bar(center_rad, height=tl['bar_out'] - tl['bar_in'],
               width=span_rad, bottom=tl['bar_in'],
               color=geno_color, edgecolor='none', linewidth=0)
        if n_pts > 0:
            seg_w = span_rad / n_pts
            vals = chr_hm[f'Track_{t_name}'].values
            thetas = start_rad + (np.arange(n_pts) + 0.5) * seg_w
            colors = [cmap(norm(v)) for v in vals]
            heights = np.full(n_pts, tl['hm_out'] - tl['hm_in'])
            bottoms = np.full(n_pts, tl['hm_in'])
            widths = np.full(n_pts, seg_w * 1.02)
            ax.bar(thetas, heights, width=widths, bottom=bottoms,
                   color=colors, edgecolor='none', linewidth=0)
    # Tick marks (OUTSIDE the chromosome ring)
    for pos in range(0, int(chr_len) + 1, 5):
        if pos > chr_len:
            break
        frac = pos / chr_len
        angle = start_rad + frac * span_rad
        ax.plot([angle, angle], [tick_bot_r, tick_top_r],
                color='black', lw=0.7, solid_capstyle='butt')
        rot = tangent_rotation(np.degrees(angle))
        ax.text(angle, tick_top_r + 0.010, str(int(pos)),
                ha='center', va='center', fontsize=4, rotation=rot)
    for pos in range(0, int(chr_len) + 1):
        if pos % 5 == 0 or pos > chr_len:
            continue
        frac = pos / chr_len
        angle = start_rad + frac * span_rad
        ax.plot([angle, angle], [tick_bot_r, tick_bot_r + tick_h * 0.5],
                color='black', lw=0.35, solid_capstyle='butt')
    # Chromosome name (outside the ticks)
    lbl_r = tick_top_r + 0.045
    lbl_angle_rad = deg2rad(layout['center'])
    rot = tangent_rotation(layout['center'])
    ax.text(lbl_angle_rad, lbl_r, chrom_name,
            ha='center', va='center',
            fontsize=11, fontweight='bold', fontstyle='italic', rotation=rot)
# Track labels A-E
gap_center_deg = chr_layout[chr_info['Chromosome'].iloc[-1]]['end'] + gap_deg / 2.0
gap_center_rad = deg2rad(gap_center_deg)
for t_name in tracks:
    tl = track_layout[t_name]
    mid_r = (tl['hm_out'] + tl['bar_in']) / 2.0
    ax.text(gap_center_rad, mid_r, t_name,
            ha='center', va='center', fontsize=10, fontweight='bold', color='black')
# ========== LEGENDS ON THE RIGHT SIDE ==========
# --- Density legend (upper right) ---
leg1 = fig.add_axes([0.76, 0.48, 0.20, 0.42])
leg1.axis('off')
legend_values = [5.3, 9.7, 14.0, 18.3, 22.7, 27.0, 31.3, 35.7, 40.0]
box_h = 0.065
box_w = 0.18
spacing = box_h + 0.015
leg1.text(0.05, 0.98, 'Density', fontsize=14, fontweight='bold', va='top', ha='left')
for i, val in enumerate(legend_values):
    y = 0.90 - i * spacing
    c = cmap(norm(val))
    leg1.add_patch(plt.Rectangle(
        (0.05, y - box_h), box_w, box_h,
        facecolor=c, edgecolor='black', linewidth=0.5))
    leg1.text(0.05 + box_w + 0.04, y - box_h / 2.0,
              f'{val}', va='center', ha='left', fontsize=11)
leg1.set_xlim(0, 1)
leg1.set_ylim(-0.05, 1.02)
# --- Genotype legend (lower right) ---
leg2 = fig.add_axes([0.76, 0.12, 0.20, 0.30])
leg2.axis('off')
leg2.text(0.05, 0.98, 'Genotype', fontsize=14, fontweight='bold', va='top', ha='left')
geno_spacing = 0.16
for i, t_name in enumerate(tracks):
    y = 0.85 - i * geno_spacing
    gc = genotype_colors[t_name]
    leg2.add_patch(plt.Rectangle(
        (0.05, y - 0.08), box_w, 0.08,
        facecolor=gc, edgecolor='black', linewidth=0.5))
    leg2.text(0.05 + box_w + 0.04, y - 0.04,
              f'Genotype {t_name}', va='center', ha='left', fontsize=11)
leg2.set_xlim(0, 1)
leg2.set_ylim(-0.05, 1.02)
# --- Save ---
plt.savefig('circos_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('circos_plot.pdf', bbox_inches='tight', facecolor='white')
print("Saved: circos_plot.png, circos_plot.pdf")
plt.close()
print("Done!")
