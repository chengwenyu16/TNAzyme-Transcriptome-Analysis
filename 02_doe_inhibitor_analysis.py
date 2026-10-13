#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TNAzyme Cleavable Motif Analysis for DOE Inhibitor Targets
Targeted for Code Ocean reproduction pipeline.
"""

import csv
import re
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

# Code Ocean 标准路径设置
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

THRESHOLD = 0.61
MOTIFS = ["ATGCAG","GTGCAG","AAGCAG","AGGCAG","ATGCAA","TTGCAG","ATGCAT","ATGCAC"]
MOTIF_RNA = ["AUG^CAG","GUG^CAG","AAG^CAG","AGG^CAG","AUG^CAA","UUG^CAG","AUG^CAU","AUG^CAC"]

def read_doe(path):
    targets = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            score = float(row["Pred. inhibitor"])
            if score >= THRESHOLD:
                targets[row["Gene"]] = score
    return targets

def read_select_ids(path):
    select_ids = set()
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        i_status = header.index("MANE_status")
        i_nuc    = header.index("Ensembl_nuc")
        for line in f:
            c = line.strip().split("\t")
            if len(c) > max(i_status, i_nuc) and c[i_status] == "MANE Select":
                select_ids.add(c[i_nuc])
    return select_ids

def parse_fasta(path, select_ids):
    genes = {}
    tid = gene = hdr = None
    seq = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if tid and tid in select_ids and gene not in genes:
                    genes[gene] = (hdr, "".join(seq))
                tid  = line.split()[0][1:]
                m    = re.search(r"gene_symbol:(\S+)", line)
                gene = m.group(1) if m else None
                hdr  = line
                seq  = []
            else:
                seq.append(line)
        if tid and tid in select_ids and gene not in genes:
            genes[gene] = (hdr, "".join(seq))
    return genes

def count_motif(seq, motif):
    n = start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1: 
            break
        n += 1
        start = pos + 1
    return n

def scan(seq):
    counts = [count_motif(seq, m) for m in MOTIFS]
    return sum(counts), counts

if __name__ == "__main__":
    targets    = read_doe(DATA_DIR / "doe.csv")
    select_ids = read_select_ids(DATA_DIR / "MANE.GRCh38.v1.5.summary.txt")
    mane       = parse_fasta(DATA_DIR / "MANE.GRCh38.v1.5.select_ensembl_rna.txt", select_ids)

    mapped, unmapped = {}, []
    for gene, score in targets.items():
        if gene in mane:
            hdr, seq = mane[gene]
            total, counts = scan(seq)
            mapped[gene] = (score, total, counts, total > 0)
        else:
            unmapped.append((gene, score))

    n = len(mapped)
    n_cut = sum(1 for v in mapped.values() if v[3])
    n_sites = sum(v[1] for v in mapped.values())

    wb = Workbook()

    hdr_font  = Font(name="Calibri", bold=True, color="FFFFFF", size=11)
    hdr_fill  = PatternFill("solid", fgColor="2E4057")
    hdr_align = Alignment(horizontal="center", vertical="center", wrap_text=True)
    thin_border = Border(left=Side(style="thin"), right=Side(style="thin"),
                         top=Side(style="thin"), bottom=Side(style="thin"))

    def style_header(ws, ncols):
        for col in range(1, ncols + 1):
            c = ws.cell(row=1, column=col)
            c.font = hdr_font
            c.fill = hdr_fill
            c.alignment = hdr_align
            c.border = thin_border

    def autosize(ws, min_w=10, max_w=45):
        for col_cells in ws.columns:
            lengths = [len(str(c.value)) for c in col_cells if c.value is not None]
            w = min(max(max(lengths, default=min_w) + 2, min_w), max_w)
            ws.column_dimensions[col_cells[0].column_letter].width = w

    # Sheet 1: Overview
    ws1 = wb.active
    ws1.title = "Overview"
    overview = [
        ("DNAzyme Cleavable Motif Analysis — DOE Inhibitor Targets", ""),
        ("(MANE Select ONLY — excluding Plus Clinical)", ""),
        ("", ""),
        ("Parameter", "Value"),
        ("Input file", "doe.csv"),
        ("MANE FASTA", "MANE.GRCh38.v1.5.select_ensembl_rna.txt"),
        ("MANE summary", "MANE.GRCh38.v1.5.summary.txt"),
        ("MANE dataset", "MANE Select only (19,363 transcripts)"),
        ("Threshold", f"Pred. inhibitor >= {THRESHOLD}"),
        ("Number of motifs", str(len(MOTIFS))),
        ("", ""),
        ("Total target genes", len(targets)),
        ("Mapped to MANE Select", n),
        ("Mapping rate", f"{n/len(targets)*100:.2f}%"),
        ("Unmapped (not analyzed)", len(unmapped)),
        ("", ""),
        ("Cleavable (>=1 motif site)", n_cut),
        ("Cleavable rate", f"{n_cut/n*100:.1f}%"),
        ("Not cleavable", n - n_cut),
        ("Total cleavage sites", n_sites),
        ("Average sites per gene", f"{n_sites/n:.1f}"),
    ]
    for row in overview:
        ws1.append(row)
    ws1["A1"].font = Font(bold=True, size=14)
    ws1["A2"].font = Font(italic=True, size=11, color="666666")
    style_header(ws1, 2)
    autosize(ws1)

    # Sheet 2: Gene Summary
    ws2 = wb.create_sheet("Gene Summary")
    headers = ["Gene", "Pred. inhibitor score", "Total sites",
               "Cleavable"] + [f"{d} ({r})" for d, r in zip(MOTIFS, MOTIF_RNA)]
    ws2.append(headers)
    for gene in sorted(mapped):
        score, total, counts, has_any = mapped[gene]
        ws2.append([gene, round(score, 4), total, "Yes" if has_any else "No"] + counts)
    style_header(ws2, len(headers))
    ws2.freeze_panes = "A2"
    autosize(ws2)

    # Sheet 3: Motif Summary
    ws3 = wb.create_sheet("Motif Summary")
    ws3.append(["Motif (DNA)", "Motif (RNA)", "Genes with motif",
                "Gene coverage (%)", "Total sites"])
    for i, (d, r) in enumerate(zip(MOTIFS, MOTIF_RNA)):
        g = sum(1 for v in mapped.values() if v[2][i] > 0)
        s = sum(v[2][i] for v in mapped.values())
        ws3.append([d, r, g, f"{g/n*100:.1f}%", s])
    ws3.append([])
    ws3.append(["Total genes analyzed", "", n, "", ""])
    ws3.append(["Cleavable genes", "", n_cut, f"{n_cut/n*100:.1f}%", n_sites])
    style_header(ws3, 5)
    autosize(ws3)

    # Sheet 4: Unmapped Genes
    ws4 = wb.create_sheet("Unmapped Genes")
    ws4.append(["Gene", "Pred. inhibitor score"])
    for gene, score in sorted(unmapped):
        ws4.append([gene, round(score, 4)])
    style_header(ws4, 2)
    ws4.freeze_panes = "A2"
    autosize(ws4)

    # Sheet 5: Methods
    ws5 = wb.create_sheet("Methods")
    methods = [
        ("Methods", ""),
        ("", ""),
        ("1. Target selection",
         f"Genes with 'Pred. inhibitor' score >= {THRESHOLD} in the DOE dataset were selected."),
        ("2. Transcript mapping",
         "Each target gene was mapped to its MANE Select representative transcript (GRCh38 v1.5). "
         "The MANE summary file was used to filter exclusively MANE Select transcripts (19,363), "
         "excluding 74 MANE Plus Clinical transcripts."),
        ("3. Motif scanning",
         "cDNA sequences were scanned for 8 TNAzyme-cleavable motifs: "
         "ATGCAG, GTGCAG, AAGCAG, AGGCAG, ATGCAA, TTGCAG, ATGCAT, ATGCAC. "
         "Overlapping matches were counted."),
        ("4. Cleavability",
         "A target was 'cleavable' if >=1 motif site was found."),
        ("5. Limitation",
         "Genes not in MANE Select (mitochondrial, pseudogenes, etc.) were not analyzed."),
    ]
    for row in methods:
        ws5.append(row)
    ws5["A1"].font = Font(bold=True, size=14)
    autosize(ws5)

    # 输出保存到 results 目录
    output_excel_path = RESULTS_DIR / "Tz2_DOE_inhibitor_MANE_Select_8motifs.xlsx"
    wb.save(output_excel_path)
    print(f"Done: {n} genes mapped, {n_cut} cleavable, {n_sites} sites")
    print(f"Saved result to: {output_excel_path}")