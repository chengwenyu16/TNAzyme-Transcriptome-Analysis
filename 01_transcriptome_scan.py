# -*- coding: utf-8 -*-
"""
Transcriptome Motif Scanner
Targeted for Code Ocean reproduction pipeline.
"""

import gzip
import os
from pathlib import Path
from Bio import SeqIO

# Code Ocean 标准路径设置
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

def count_motif_overlapping(seq, motif):
    n = start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1: 
            break
        n += 1
        start = pos + 1
    return n

def scan_transcriptome(fasta_file, motifs):
    total_transcripts = 0
    targetable_transcripts = 0
    
    transcript_counts = {motif: 0 for motif in motifs}
    site_counts = {motif: 0 for motif in motifs}
    total_sites = 0
    
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_transcripts += 1
            sequence = str(record.seq).upper()
            is_targetable = False
            
            for motif in motifs:
                count = count_motif_overlapping(sequence, motif)
                if count > 0:
                    is_targetable = True
                    transcript_counts[motif] += 1
                    site_counts[motif] += count
                    total_sites += count
                    
            if is_targetable:
                targetable_transcripts += 1
                
    return total_transcripts, targetable_transcripts, transcript_counts, site_counts, total_sites

if __name__ == "__main__":
    fasta_path = DATA_DIR / "gencode.v50.transcripts.fa.gz"
    
    target_motifs = [
        "ATGCAG", "GTGCAG", "AAGCAG", "AGGCAG",
        "ATGCAA", "TTGCAG", "ATGCAT", "ATGCAC"
    ]
    
    total, targeted, motif_transcripts, motif_sites, overall_sites = scan_transcriptome(fasta_path, target_motifs)
    overall_cov = (targeted / total) * 100 if total > 0 else 0.0
    
    summary_lines = [
        "=" * 75,
        "OVERALL SUMMARY",
        "=" * 75,
        f"Total human mRNA transcripts analyzed: {total}",
        f"Transcripts containing >=1 targetable motif: {targeted} ({overall_cov:.2f}%)",
        f"Total cleavage sites across all transcripts: {overall_sites}\n",
        "=" * 75,
        "BREAKDOWN BY MOTIF",
        "=" * 75,
        f"{'Motif':<10} | {'Targeted Transcripts':<22} | {'Coverage (%)':<15} | {'Total Sites (Cuts)':<20}",
        "-" * 75
    ]
    
    for motif in target_motifs:
        t_count = motif_transcripts[motif]
        s_count = motif_sites[motif]
        cov = (t_count / total) * 100 if total > 0 else 0.0
        summary_lines.append(f"{motif:<10} | {t_count:<22} | {cov:<15.2f} | {s_count:<20}")
    summary_lines.append("=" * 75)
    
    output_text = "\n".join(summary_lines)
    print(output_text)
    
    # 将终端打印的摘要同步保存到 results 目录中供下载
    with open(RESULTS_DIR / "transcriptome_scan_summary.txt", "w", encoding="utf-8") as out:
        out.write(output_text)