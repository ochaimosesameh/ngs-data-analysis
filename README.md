# NGS Data Analysis

**Author:** Ochai Moses Ameh

## Project Overview

This project is a beginner-friendly Python tool for analysing
Next-Generation Sequencing (NGS) reads and their sequencing quality.

## Features

- Reads FASTQ files.
- Counts A, T, G, C, and N bases.
- Calculates read length and GC content.
- Converts Phred+33 characters into quality scores.
- Checks sequence validity.
- Identifies reads requiring quality-control review.
- Generates individual reports and an overall summary.

## Files

- `ngs_analyzer.py`: Main Python program.
- `sample_reads.fastq`: Sample sequencing data.

## How to Run

Run `ngs_analyzer.py` using Python 3.

When prompted, enter the FASTQ filename:
`sample_reads.fastq`

Ensure the FASTQ file is available in the program's working directory.

## Limitations

The quality-control thresholds are illustrative and are not
universal standards for real sequencing datasets.

This project is intended for learning and demonstration.

## Author

Ochai Moses Ameh