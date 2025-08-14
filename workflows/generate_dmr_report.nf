process report_dmrs {
    label 'dmr_report'
    label 'process_default'
    publishDir "${params.output_dir}/annotated_dmrs/report", mode: 'copy'

    input:
    path annotated_bed
    path annotation_log

    output:
    path "dmr_summary_report.html", emit: report

    script:
    """
#!/usr/bin/env python
import os
os.environ['MPLCONFIGDIR'] = os.getcwd()
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import base64
import io

# Read annotation summary if available
summary_text = ""
if os.path.exists("${annotation_log}"):
    with open("${annotation_log}", 'r') as f:
        summary_text = f.read().replace('\\n', '<br>')

# Read DMR data
try:
    df = pd.read_csv("${annotated_bed}", sep='\\t')
    
    # Generate plots
    plots = []
    
    # Plot 1: Annotation counts
    if 'annotation' in df.columns:
        # Combine promoter annotations into a single category
        df['annotation'] = df['annotation'].replace(['promoter_plus', 'promoter_minus'], 'promoters')
        fig, ax = plt.subplots(figsize=(6, 4))
        counts = df['annotation'].value_counts().head(10).sort_values(ascending=True)
        # Slightly darker, more vibrant color palette
        color_palette = ['#73a3d4', '#ff9d5c', '#66c37f', '#ff7e79', '#b194e6', '#c8a17e', '#f78ce0', '#b1b1b1', '#f2f17b', '#92e4e1']
        plot_colors = color_palette[:len(counts)]
        bars = ax.barh(counts.index, counts.values, color=plot_colors)
        ax.set_xlabel('Count')
        ax.set_title('DMRs by Genomic Annotation')
        for bar, val in zip(bars, counts.values):
            ax.text(bar.get_width()+0.5, bar.get_y()+bar.get_height()/2, 
                    str(val), va='center')
        plt.tight_layout()
        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        plot1 = base64.b64encode(buf.read()).decode()
        plots.append(plot1)
        plt.close()
    
    # Plot 2: Biotype distribution
    if 'biotype' in df.columns and not df['biotype'].dropna().empty:
        fig, ax = plt.subplots(figsize=(8, 5), subplot_kw=dict(aspect="equal"))
        
        biotype_counts = df['biotype'].value_counts()
        
        # Group small slices into 'Other' to keep the plot clean
        threshold = biotype_counts.sum() * 0.02 # Group if less than 2%
        other_count = biotype_counts[biotype_counts < threshold].sum()
        biotype_counts_main = biotype_counts[biotype_counts >= threshold]
        if other_count > 0:
            biotype_counts_main.loc['Other'] = other_count
        
        # Darker palette of blue, orange, and green
        color_palette = ['#4682B4', '#FF8C00', '#2E8B57', '#5F9EA0', '#FFA500', '#3CB371', '#6495ED', '#FFD700']
        colors = color_palette[:len(biotype_counts_main)]

        wedges, texts, autotexts = ax.pie(biotype_counts_main, 
                                          autopct='%1.1f%%',
                                          startangle=90,
                                          pctdistance=0.85,
                                          colors=colors,
                                          wedgeprops=dict(width=0.4, edgecolor='w'))

        ax.legend(wedges, biotype_counts_main.index,
                  title="Biotypes",
                  loc="center left",
                  bbox_to_anchor=(1, 0, 0.5, 1))

        plt.setp(autotexts, size=10, weight="bold")
        ax.set_title("Distribution of DMRs by Gene Biotype")
        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        plot2 = base64.b64encode(buf.read()).decode()
        plots.append(plot2)
        plt.close()
    
    # Format table
    display_cols = ['chr', 'start', 'end', 'length', 'nSites', 'meanMethy1', 
                    'meanMethy2', 'diff.Methy', 'annotation', 'biotype', 'gene']
    display_cols = [c for c in display_cols if c in df.columns]
    df_display = df[display_cols].drop_duplicates()
    
    # Format numeric columns
    if 'meanMethy1' in df_display.columns:
        df_display['meanMethy1'] = df_display['meanMethy1'].round(3)
    if 'meanMethy2' in df_display.columns:
        df_display['meanMethy2'] = df_display['meanMethy2'].round(3)
    if 'diff.Methy' in df_display.columns:
        df_display['diff.Methy'] = df_display['diff.Methy'].round(3)
    
    table_html = df_display.to_html(index=False, table_id="dmr_table", 
                                    classes="display", escape=False)
    
    # Create plots HTML
    plots_html = '<div style="display:flex;justify-content:center;gap:20px;margin:20px 0;">'
    plots_html = '<div class="plot-container">'
    for plot_b64 in plots:
        plots_html += f'<div class="plot-box"><img src="data:image/png;base64,{plot_b64}" alt="summary plot"></div>'
    plots_html += '</div>'

except Exception as e:
    table_html = f"<p>Error processing data: {str(e)}</p>"
    plots_html = ""

# Generate HTML
html = f'''<!DOCTYPE html>
<html>
<head>
<title>DMR Report</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<style>
    body {{font-family:"Segoe UI", Arial, sans-serif; margin:20px; background:#f5f5f5;}}
    .page-container {{max-width: 1600px; margin: auto;}}
    h1 {{color:#333; text-align:center; border-bottom: 1px solid #ccc; padding-bottom: 10px; margin-bottom: 30px;}}
    h2 {{color:#333; border-bottom: 1px solid #eee; padding-bottom: 8px; margin-top: 0; margin-bottom: 20px; font-size: 1.5em;}}
    .content-box {{
        background: white;
        padding: 25px;
        margin-bottom: 30px;
        border-radius: 8px;
        box-shadow: 0 4px 8px rgba(0,0,0,0.05);
        border: 1px solid #e9ecef;
    }}
    .summary-text {{background:#e8f4f8; padding:15px; border-radius:5px; font-family:monospace;}}
    .plot-container {{display:flex; flex-wrap:wrap; justify-content:center; gap:30px;}}
    .plot-box {{padding:15px; border:1px solid #ddd; border-radius:8px; background-color:#fdfdfd; flex:1; min-width:400px; max-width:48%; text-align:center;}}
    .plot-box img {{max-width:100%; height:auto;}}
    table.dataTable {{width:100%!important; border-collapse: collapse;}}
    table.dataTable th, table.dataTable td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #dddddd; border-left: none; border-right: none;}}
    table.dataTable thead th {{ background-color: #4CAF50; color: white; border-bottom: 2px solid #45a049; }}
    table.dataTable tbody tr:hover {{ background-color: #f1f1f1; }}
</style>
</head>
<body>
<div class="page-container">
    <h1>Differentially Methylated Regions (DMRs) Summary</h1>
    <div class="content-box">
        <h2>Analysis Summary</h2>
        <div class="summary-text">
            {summary_text if summary_text else "No summary available"}
        </div>
    </div>
    <div class="content-box">
        <h2>Summary Plots</h2>
        {plots_html}
    </div>
    <div class="content-box">
        <h2>DMR Table</h2>
        {table_html}
    </div>
</div>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script>
    \$(document).ready(function() {{
        \$('#dmr_table').DataTable({{"pageLength": 25, "order": [[3, "desc"]]}});
    }});
</script>
</body>
</html>'''

with open("dmr_summary_report.html", "w") as f:
    f.write(html)
"""
}
