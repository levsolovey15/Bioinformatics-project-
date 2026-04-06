# Variant Prioritizer

A Python CLI tool for prioritizing genomic variants from annotated VCF files.

## Features
- Supports VCF and VCF.GZ
- Parses ANN / CSQ annotations
- Scores variants based on:
  - impact (HIGH, MODERATE, LOW)
  - allele frequency
  - gene panel
  - clinical significance

## Usage

```bash
python variant_prioritizer.py --vcf input.vcf --out results.tsv