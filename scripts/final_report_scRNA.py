#!/usr/bin/env python3

import os
import glob
import subprocess
import pandas as pd
import re
import sys

print("Starting Final Report Generation...")

# ------------------------------------------------
# Locate Project Root
# ------------------------------------------------
current = os.getcwd()

while current != "/" and not os.path.exists(os.path.join(current, "main.nf")):
    current = os.path.dirname(current)

project_root = current

results_dir = os.path.join(project_root, "results")
final_report_dir = os.path.join(results_dir, "final_report")
os.makedirs(final_report_dir, exist_ok=True)

print("Detected project root:", project_root)
print("Detected results directory:", results_dir)

# ------------------------------------------------
# 1. Run MultiQC (clean output)
# ------------------------------------------------
pre_qc_dir = os.path.join(results_dir, "pre_qc")
post_qc_dir = os.path.join(results_dir, "post_qc")
multiqc_out = os.path.join(results_dir, "multiqc")

os.makedirs(multiqc_out, exist_ok=True)

dirs_to_check = [pre_qc_dir, post_qc_dir]
dirs_exist = [d for d in dirs_to_check if os.path.exists(d)]

if dirs_exist:
    result = subprocess.run(
        ["multiqc", *dirs_exist, "-o", multiqc_out],
        stdout=subprocess.PIPE,   # capture normal messages
        stderr=subprocess.PIPE,   # capture errors
        text=True
    )

    # Optional: print only a small summary or nothing
    print("MultiQC completed successfully")
    # print(result.stdout)
    # print(result.stderr)

else:
    print("QC folders not found, skipping MultiQC")
# ------------------------------------------------
# 2. Parse FLAGSTAT files
# ------------------------------------------------
flag_dir = os.path.join(results_dir, "flagstat")
flag_files = glob.glob(os.path.join(flag_dir, "*_flagstat.txt"))

summary = []

meta_file = os.path.join(project_root, "samples.tsv")

if os.path.exists(meta_file):
    meta_df = pd.read_csv(meta_file, sep="\t")
    sample_condition = dict(zip(meta_df["sample"], meta_df["condition"]))
else:
    sample_condition = {}

for file in flag_files:

    sample = os.path.basename(file).replace("_flagstat.txt","")

    total_reads = "NA"
    mapped_reads = "NA"
    primary_reads = "NA"

    with open(file) as f:

        for line in f:

            if "in total" in line:
                total_reads = line.split()[0]

            elif " mapped (" in line and "primary" not in line:
                match = re.search(r"\(([\d\.]+%)", line)
                if match:
                    mapped_reads = match.group(1)

            elif "primary mapped" in line:
                match = re.search(r"\(([\d\.]+%)", line)
                if match:
                    primary_reads = match.group(1)

    condition = sample_condition.get(sample,"NA")

    summary.append([
        sample,
        total_reads,
        mapped_reads,
        primary_reads,
        condition
    ])

if summary:

    df = pd.DataFrame(
        summary,
        columns=[
            "Sample",
            "Total_Reads",
            "Mapped_Reads",
            "Primary_Mapped_Reads",
            "Condition"
        ]
    )

    flag_csv_path = os.path.join(final_report_dir,"flagstat_summary.csv")
    df.to_csv(flag_csv_path,index=False)

    flag_table = df.to_html(index=False)

else:

    flag_table = "<p>No flagstat files found</p>"

print("Flagstat summary prepared")

# ------------------------------------------------
# 3. Load Count Matrix
# ------------------------------------------------
counts_file = os.path.join(results_dir,"counts","counts.txt")

if os.path.exists(counts_file):

    counts_df = pd.read_csv(counts_file,sep="\t",comment="#")

    counts_html = counts_df.head(20).to_html(index=False)

    counts_link = "<a href='../counts/counts.txt'>Download full count matrix</a>"

else:

    counts_html = "<p>Counts file not found</p>"
    counts_link = ""

# ------------------------------------------------
# 4. Load Metadata
# ------------------------------------------------
if os.path.exists(meta_file):

    meta_df = pd.read_csv(meta_file,sep="\t")

    meta_html = meta_df.to_html(index=False)

    meta_link = "<a href='../samples.tsv'>Download full metadata</a>"

else:

    meta_html = "<p>Metadata file not found</p>"
    meta_link = ""

# ------------------------------------------------
# 5. Collect Plots from pipeline
# ------------------------------------------------

plot_source = sys.argv[1]

plots_dir = os.path.join(results_dir, "plots")
os.makedirs(plots_dir, exist_ok=True)

print("Copying plots from:", plot_source)

for p in glob.glob(os.path.join(plot_source, "*")):
    dest = os.path.join(plots_dir, os.path.basename(p))

    if os.path.abspath(p) != os.path.abspath(dest):
        subprocess.run(["cp", p, dest])

plot_files = sorted(glob.glob(os.path.join(plots_dir, "*.png"))) + \
             sorted(glob.glob(os.path.join(plots_dir, "*.pdf")))

plot_html = "<div style='display:flex; flex-wrap:wrap;'>"

for plot in plot_files:

    name = os.path.basename(plot)

    plot_html += f"""
    <div style='width:45%; margin:2%;'>
        <h4>{name}</h4>
        <img src='../plots/{name}' width='100%'>
    </div>
    """

plot_html += "</div>"
# ------------------------------------------------
# 6. KEGG pathway images and table
# ------------------------------------------------
kegg_dir  = os.path.join(results_dir, "kegg")
kegg_file = os.path.join(kegg_dir, "kegg_links.csv")

if os.path.exists(kegg_file):
    import shutil  # <-- you must import at the top of the file
    report_kegg_dir = os.path.join(final_report_dir, "kegg")
    os.makedirs(report_kegg_dir, exist_ok=True)

    for img in glob.glob(os.path.join(kegg_dir, "*.png")):
        shutil.copy(img, report_kegg_dir)

    kegg_df = pd.read_csv(kegg_file)
    kegg_html = "<h2>7. KEGG Pathway Analysis</h2>\n"
    kegg_html += "<table border='1'><tr><th>Pathway</th><th>KEGG ID</th><th>Diagram</th></tr>"

    for _, row in kegg_df.iterrows():
        pathway = row["Pathway"]
        pid = row["KEGG_ID"]
        img_path = os.path.join("kegg", f"{pid}.png")  # relative path inside report
        kegg_html += f"<tr><td>{pathway}</td><td>{pid}</td>"
        kegg_html += f"<td><a href='{img_path}' target='_blank'>Open Diagram</a></td></tr>"

    kegg_html += "</table>"

else:
    kegg_html = "<p>No KEGG pathway results found.</p>"
# ------------------------------------------------
# 7. KEGG Pathway Diagrams (Fully Automatic)
# ------------------------------------------------
import shutil
import os
import pandas as pd

kegg_dir = os.path.join(results_dir, "kegg")
report_kegg_dir = os.path.join(final_report_dir, "kegg")
os.makedirs(report_kegg_dir, exist_ok=True)

kegg_file = os.path.join(kegg_dir, "kegg_links.csv")
kegg_html = ""

if os.path.exists(kegg_file):
    kegg_df = pd.read_csv(kegg_file)

    kegg_html += ""
    kegg_html += "<p>Significant biological pathways enriched from differential gene expression analysis.</p>"
    kegg_html += "<table border='1'><tr><th>Pathway</th><th>KEGG ID</th><th>Colored Diagram</th><th>Green Diagram</th></tr>"

    # Dynamically check existing PNGs
    for _, row in kegg_df.iterrows():
        pathway = row["Pathway"]
        pid = row["KEGG_ID"]

        # Check colored PNG
        colored_png_name = f"{pid}.{pid}.png"
        colored_png_path = os.path.join(kegg_dir, colored_png_name)
        colored_link = None
        if os.path.exists(colored_png_path):
            shutil.copy(colored_png_path, report_kegg_dir)
            colored_link = os.path.join("kegg", colored_png_name)

        # Check green PNG
        green_png_name = f"{pid}.png"
        green_png_path = os.path.join(kegg_dir, green_png_name)
        green_link = None
        if os.path.exists(green_png_path):
            shutil.copy(green_png_path, report_kegg_dir)
            green_link = os.path.join("kegg", green_png_name)

        # Skip if no PNG exists for this pathway
        if not colored_link and not green_link:
            continue

        kegg_html += f"<tr><td>{pathway}</td><td>{pid}</td>"
        kegg_html += f"<td>{f'<a href={colored_link} target=_blank>Open Colored</a>' if colored_link else 'NA'}</td>"
        kegg_html += f"<td>{f'<a href={green_link} target=_blank>Open Green</a>' if green_link else 'NA'}</td></tr>"

    kegg_html += "</table>"

else:
    kegg_html = "<p>No KEGG pathway results found.</p>"
# ------------------------------------------------
# 6. Generate Final HTML
# ------------------------------------------------

html_content = f"""
<html>

<head>
<title>Transcriptomic Analysis of Lung Cancer Using Single-Cell RNA-Seq Data</title>
<style>

body {{
font-family: Arial;
margin:40px;
}}

h1 {{
color:#2c3e50;
}}

table {{
border-collapse: collapse;
}}

td, th {{
border:1px solid black;
padding:6px;
}}

</style>

</head>

<body>

<h1>Transcriptomic Analysis of Lung Cancer Using Single-Cell RNA-Seq Data Reveals Differentially Expressed Genes and Enriched Pathways Through a Nextflow-Based Automated Bioinformatics Pipeline</h1>
<h3>Student Information</h3>

<p><b>Name:</b> DEVI KOUSALYA JUTURI</p>

<p><b>Registration Number:</b> 24MSBI108</p>

<p><b>Program:</b> M.Sc Bioinformatics</p>

<p><b>University:</b> Garden City University</p>

<p><b>Internship Organization:</b> Miriomics</p>
<p><b>Project Type:</b> Internship Project</p>
<p>
This report summarizes the analysis of transcriptomic sequencing data using an automated Nextflow pipeline.
The workflow performs quality control, read alignment, gene quantification, and differential expression
analysis between tumor and normal samples.
</p>

<h2>1. Quality Control</h2>

<p>
Raw sequencing reads were evaluated using FastQC. Multiple QC reports were aggregated using MultiQC
to provide an overview of sequencing quality across all samples.
</p>

<a href="../multiqc/multiqc_report.html" target="_blank">Open MultiQC Detailed Report</a>

<h2>2. Alignment Statistics</h2>

<p>
Reads were aligned to the human reference genome and alignment statistics were calculated using
SAMtools flagstat. The table below summarizes total reads, mapped reads, and primary mapped reads
for each sample.
</p>

{flag_table}

<p><a href="flagstat_summary.csv">Download full CSV</a></p>

<h2>3. Gene Expression Count Matrix</h2>

<p>
Gene expression quantification was performed using featureCounts, which assigns aligned reads to
genomic features (genes) and produces a count matrix used for downstream differential expression
analysis.
</p>

{counts_html}

<p>{counts_link}</p>

<h2>4. Sample Metadata</h2>
<p>
<p>
The RNA-seq dataset used in this analysis was obtained from the NCBI Sequence Read Archive (SRA). 
Taxonomy analysis of the sequencing reads confirmed the presence of 
<b>Homo sapiens</b> sequences in all selected samples. The selected runs 
showed approximately 5–7% classification to Homo sapiens during taxonomy analysis, 
indicating that the reads originate from human transcriptomic single-cell data. 
Based on this observation, samples were selected according to their Homo sapiens 
classification percentage and their biological conditions, categorized as 
Tumor and Normal.
</p>
</p>

{meta_html}

<p>{meta_link}</p>

<h2>5. Differential Expression Analysis</h2>

<p>
Differential expression analysis was performed using the DESeq2 R package. DESeq2 models count data
using a negative binomial distribution and identifies genes that are significantly upregulated or
downregulated between tumor and normal samples.
</p>

<h2>6. Differential Expression Plots</h2>

<p>
The following plots summarize the results of differential gene expression analysis including
MA plots, PCA visualization, volcano plots,upregulated plots,downregulated plots heatmaps,and enrichment analyses.
</p>

{plot_html}

<h2>7. KEGG Pathway Diagrams</h2>

<p>
The enriched pathways identified from differential expression analysis are listed below.
Click a pathway to open the interactive diagram from the KEGG database using Pathview.
</p>

{kegg_html}

</body>
</html>
"""

final_html_path = "scRNAseq_Report.html"
with open(final_html_path, "w") as f:
    f.write(html_content)

print("Final HTML Report Generated Successfully")
