# NGS Data Analysis

**Author:** Ochai Moses Ameh

## Project Overview

This project is a beginner-friendly Python tool for analysing illustrative Next-Generation Sequencing (NGS) reads and assessing basic sequencing quality.

It demonstrates fundamental bioinformatics concepts using sample FASTQ data.

## Objectives

- Analyse individual sequencing reads.
- Analyse multiple sequencing reads.
- Calculate read lengths and nucleotide counts.
- Calculate GC content.
- Convert FASTQ quality characters into Phred quality scores.
- Calculate average sequencing quality.
- Identify invalid DNA characters.
- Apply illustrative quality-control checks.
- Generate an overall quality-control summary.

## Features

### 1. Single-Read Analysis

Users can enter a DNA sequencing read and optionally provide its FASTQ quality string.

### 2. Multiple-Read Analysis

The program analyses four built-in example FASTQ reads and generates a report for each read.

### 3. Sequence Statistics

For each read, the program reports:

- Read length
- A, T, G, C and N nucleotide counts
- GC content
- Average Phred quality score, when available

### 4. Quality Control

The program checks read length, GC content and average quality against illustrative thresholds.

It also identifies invalid characters and reports reads that require review.

### 5. Overall Summary

For multiple reads, the program reports:

- Total number of reads
- Valid and invalid reads
- Reads passing quality control
- Reads requiring review
- Average read length
- Average GC content
- Average Phred quality score

## Technologies Used

- Python
- FASTQ format
- Basic sequencing quality-control concepts

## How to Run

1. Open `ngs_analysis.py` in a Python 3 environment.
2. Run the program.
3. Choose option 1 to analyse a single read.
4. Choose option 2 to analyse the built-in FASTQ examples.
5. Choose option 3 to exit.

## Data and Limitations

The multiple-read analysis uses illustrative FASTQ records built into the program. It does not currently load external FASTQ files.

The quality-control thresholds are educational examples, not universal laboratory standards.

This project demonstrates basic quality-control concepts. It does not replace professional sequencing analysis tools such as FastQC.

## Future Improvements

- Support loading external FASTQ files.
- Generate visualisations of sequencing quality.
- Add more comprehensive quality-control statistics.
- Compare results across larger sequencing datasets.

## Author

Ochai Moses Ameh