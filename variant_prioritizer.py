#!/usr/bin/env python3
"""
variant_prioritizer.py

A command-line bioinformatics tool for prioritizing variants from VCF files.
Designed for GitHub portfolio projects / CV-ready presentation.

"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import sys
from dataclasses import dataclass, asdict
from typing import Dict, Generator, Iterable, List, Optional, Tuple


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


# ----------------------------
# Data models
# ----------------------------

@dataclass
class VariantRecord:
    chrom: str
    pos: int
    var_id: str
    ref: str
    alt: str
    qual: str
    flt: str
    gene: str
    consequence: str
    impact: str
    allele_frequency: Optional[float]
    clinical_significance: str
    score: float


# ----------------------------
# File utilities
# ----------------------------

def open_textfile(path: str):
    """Open plain text or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_gene_panel(path: Optional[str]) -> set[str]:
    """Load a gene panel from a text file, one gene symbol per line."""
    genes = set()
    if not path:
        return genes

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene.upper())

    logging.info("Loaded %d genes from panel.", len(genes))
    return genes


# ----------------------------
# VCF parsing
# ----------------------------

def parse_info_field(info_str: str) -> Dict[str, str]:
    """Parse INFO column into a dictionary."""
    info = {}
    if info_str == ".":
        return info

    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
        else:
            info[entry] = "True"
    return info


def extract_annotation_fields(info: Dict[str, str]) -> Tuple[str, str, str]:
    """
    Extract gene, consequence, impact from ANN or CSQ.
    Assumes common SnpEff/VEP-like formatting.
    Returns (gene, consequence, impact).
    """
    # SnpEff ANN format often:
    # Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|...
    if "ANN" in info:
        ann_entries = info["ANN"].split(",")
        if ann_entries:
            parts = ann_entries[0].split("|")
            gene = parts[3] if len(parts) > 3 else ""
            consequence = parts[1] if len(parts) > 1 else ""
            impact = parts[2] if len(parts) > 2 else ""
            return gene, consequence, impact

    # VEP CSQ format varies a lot; this is a best-effort fallback
    # Common early fields:
    # Allele|Consequence|IMPACT|SYMBOL|Gene|...
    if "CSQ" in info:
        csq_entries = info["CSQ"].split(",")
        if csq_entries:
            parts = csq_entries[0].split("|")
            consequence = parts[1] if len(parts) > 1 else ""
            impact = parts[2] if len(parts) > 2 else ""
            gene = parts[3] if len(parts) > 3 else ""
            return gene, consequence, impact

    return "", "", ""


def extract_allele_frequency(info: Dict[str, str]) -> Optional[float]:
    """Try to extract allele frequency from common INFO tags."""
    for key in ("AF", "gnomAD_AF", "ExAC_AF", "TOPMED_AF"):
        if key in info:
            raw = info[key].split(",")[0]
            try:
                return float(raw)
            except ValueError:
                continue
    return None


def extract_clinical_significance(info: Dict[str, str]) -> str:
    """Try to extract clinical significance from common INFO tags."""
    for key in ("CLNSIG", "CLIN_SIG", "CLINICAL_SIGNIFICANCE"):
        if key in info:
            return info[key]
    return ""


def parse_vcf(path: str) -> Generator[Tuple[str, int, str, str, str, str, str, Dict[str, str]], None, None]:
    """
    Yield parsed VCF records:
    (chrom, pos, id, ref, alt, qual, filter, info_dict)
    """
    with open_textfile(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, var_id, ref, alt, qual, flt, info_str = fields[:8]
            info = parse_info_field(info_str)

            for alt_allele in alt.split(","):
                yield chrom, int(pos), var_id, ref, alt_allele, qual, flt, info


# ----------------------------
# Scoring logic
# ----------------------------

def score_variant(
    gene: str,
    impact: str,
    consequence: str,
    af: Optional[float],
    clin_sig: str,
    gene_panel: set[str]
) -> float:
    """Assign a prioritization score to a variant."""
    score = 0.0

    impact_upper = impact.upper()
    consequence_lower = consequence.lower()
    clin_sig_lower = clin_sig.lower()
    gene_upper = gene.upper()

    # Impact score
    if impact_upper == "HIGH":
        score += 6.0
    elif impact_upper == "MODERATE":
        score += 3.5
    elif impact_upper == "LOW":
        score += 1.0

    # Consequence-based bonuses
    severe_terms = [
        "frameshift_variant",
        "stop_gained",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_lost"
    ]
    if any(term in consequence_lower for term in severe_terms):
        score += 4.0
    elif "missense_variant" in consequence_lower:
        score += 2.0
    elif "synonymous_variant" in consequence_lower:
        score -= 1.0

    # Frequency score: rarer variants get more points
    if af is None:
        score += 1.0
    else:
        if af < 0.0001:
            score += 4.0
        elif af < 0.001:
            score += 3.0
        elif af < 0.01:
            score += 1.5
        elif af > 0.05:
            score -= 3.0

    # Gene panel bonus
    if gene_upper and gene_upper in gene_panel:
        score += 5.0

    # Clinical significance
    if "pathogenic" in clin_sig_lower:
        score += 6.0
    elif "likely_pathogenic" in clin_sig_lower:
        score += 4.0
    elif "benign" in clin_sig_lower:
        score -= 4.0

    return round(score, 3)


# ----------------------------
# Prioritization pipeline
# ----------------------------

def prioritize_variants(vcf_path: str, gene_panel: set[str]) -> List[VariantRecord]:
    """Parse VCF and produce scored variant records."""
    results: List[VariantRecord] = []

    total = 0
    for chrom, pos, var_id, ref, alt, qual, flt, info in parse_vcf(vcf_path):
        total += 1

        gene, consequence, impact = extract_annotation_fields(info)
        af = extract_allele_frequency(info)
        clin_sig = extract_clinical_significance(info)

        score = score_variant(
            gene=gene,
            impact=impact,
            consequence=consequence,
            af=af,
            clin_sig=clin_sig,
            gene_panel=gene_panel
        )

        record = VariantRecord(
            chrom=chrom,
            pos=pos,
            var_id=var_id if var_id != "." else "",
            ref=ref,
            alt=alt,
            qual=qual,
            flt=flt,
            gene=gene,
            consequence=consequence,
            impact=impact,
            allele_frequency=af,
            clinical_significance=clin_sig,
            score=score
        )
        results.append(record)

    logging.info("Parsed %d variants.", total)

    results.sort(key=lambda x: x.score, reverse=True)
    return results


def write_tsv(records: Iterable[VariantRecord], out_path: str) -> None:
    """Write prioritized variants to TSV."""
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chrom", "pos", "var_id", "ref", "alt", "qual", "flt",
                "gene", "consequence", "impact",
                "allele_frequency", "clinical_significance", "score"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))

    logging.info("Saved output to: %s", out_path)


# ----------------------------
# CLI
# ----------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prioritize variants from a VCF file using annotation-driven scoring."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF or VCF.GZ file."
    )
    parser.add_argument(
        "--genes",
        required=False,
        default=None,
        help="Optional gene panel text file (one gene symbol per line)."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file with prioritized variants."
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        logging.error("VCF file does not exist: %s", args.vcf)
        return 1

    if args.genes and not os.path.exists(args.genes):
        logging.error("Gene panel file does not exist: %s", args.genes)
        return 1

    gene_panel = load_gene_panel(args.genes)
    prioritized = prioritize_variants(args.vcf, gene_panel)
    write_tsv(prioritized, args.out)

    if prioritized:
        logging.info("Top 5 prioritized variants:")
        for record in prioritized[:5]:
            logging.info(
                "%s:%d %s>%s | gene=%s | impact=%s | score=%.2f",
                record.chrom,
                record.pos,
                record.ref,
                record.alt,
                record.gene or "NA",
                record.impact or "NA",
                record.score
            )
    else:
        logging.warning("No variants found.")

    return 0


if __name__ == "__main__":
    sys.exit(main())