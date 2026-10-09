"""
NGS Data Analysis
Author: Ochai Moses Ameh

A beginner-friendly tool for analysing illustrative
Next-Generation Sequencing (NGS) reads and their quality.
"""

from statistics import mean

AUTHOR = "Ochai Moses Ameh"

# Illustrative QC thresholds for this learning project.
MIN_READ_LENGTH = 5
MIN_AVERAGE_QUALITY = 20
MIN_GC_CONTENT = 30
MAX_GC_CONTENT = 70


def clean_sequence(sequence):
    """Convert a DNA sequence to uppercase and remove whitespace."""
    return "".join(sequence.split()).upper()


def calculate_gc_content(sequence):
    """Calculate the percentage of G and C bases."""
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def count_bases(sequence):
    """Count the occurrences of A, T, G, C and N."""
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "N": sequence.count("N")
    }


def validate_sequence(sequence):
    """Return characters that are not accepted DNA bases."""
    valid_bases = set("ATGCN")
    return sorted(set(sequence) - valid_bases)


def phred_score_to_quality(character):
    """
    Convert a FASTQ quality character to a Phred score.
    This project uses the standard Phred+33 encoding.
    """
    return ord(character) - 33


def analyze_read(read_id, sequence, quality_string=None):
    """Analyse one sequencing read and return its statistics."""

    sequence = clean_sequence(sequence)
    invalid_bases = validate_sequence(sequence)

    result = {
        "read_id": read_id,
        "sequence": sequence,
        "valid": not invalid_bases and bool(sequence),
        "invalid_bases": invalid_bases,
        "length": len(sequence),
        "base_counts": count_bases(sequence),
        "gc_content": calculate_gc_content(sequence),
        "average_quality": None,
        "qc_status": "Not assessed"
    }

    if not sequence:
        result["error"] = "The sequence is empty."
        result["qc_status"] = "FAIL"
        return result

    if invalid_bases:
        result["error"] = (
            "Invalid characters found: "
            + ", ".join(invalid_bases)
        )
        result["qc_status"] = "FAIL"
        return result

    if quality_string is not None:
        if len(quality_string) != len(sequence):
            result["error"] = (
                "Quality-string length does not match "
                "sequence length."
            )
            result["qc_status"] = "FAIL"
            return result

        quality_scores = [
            phred_score_to_quality(character)
            for character in quality_string
        ]

        if any(score < 0 or score > 93
               for score in quality_scores):
            result["error"] = (
                "Quality string contains invalid "
                "Phred+33 scores."
            )
            result["qc_status"] = "FAIL"
            return result

        result["average_quality"] = mean(quality_scores)

    warnings = []

    if result["length"] < MIN_READ_LENGTH:
        warnings.append("Short read")

    if not (
        MIN_GC_CONTENT
        <= result["gc_content"]
        <= MAX_GC_CONTENT
    ):
        warnings.append("Review GC content")

    if result["average_quality"] is not None:
        if result["average_quality"] < MIN_AVERAGE_QUALITY:
            warnings.append("Low average quality")

    if warnings:
        result["qc_status"] = "REVIEW: " + "; ".join(warnings)
    else:
        result["qc_status"] = "PASS"

    return result


def print_read_report(result):
    """Display the results for one read."""

    print("\n" + "=" * 45)
    print("NGS READ ANALYSIS")
    print("=" * 45)

    print("Read ID:", result["read_id"])
    print("Sequence:", result["sequence"])
    print("Status:", "Valid" if result["valid"] else "Invalid")

    if result["invalid_bases"]:
        print("Invalid characters:", result["invalid_bases"])

    if "error" in result:
        print("Error:", result["error"])

    print("Read length:", result["length"])

    for base, count in result["base_counts"].items():
        print(base + " count:", count)

    print("GC content:", round(result["gc_content"], 2), "%")

    if result["average_quality"] is not None:
        print(
            "Average Phred quality:",
            round(result["average_quality"], 2)
        )
    else:
        print("Average Phred quality: Not available")

    print("QC status:", result["qc_status"])
    print("=" * 45)


def parse_fastq(fastq_text):
    """
    Parse standard four-line FASTQ records.
    Returns a list of (read_id, sequence, quality) tuples.
    """

    lines = [
        line.strip()
        for line in fastq_text.strip().splitlines()
        if line.strip()
    ]

    if len(lines) % 4 != 0:
        raise ValueError(
            "FASTQ data must contain four lines per read."
        )

    reads = []

    for index in range(0, len(lines), 4):
        header = lines[index]
        sequence = lines[index + 1]
        separator = lines[index + 2]
        quality = lines[index + 3]

        if not header.startswith("@"):
            raise ValueError(
                "FASTQ read headers must begin with @."
            )

        if not separator.startswith("+"):
            raise ValueError(
                "The third FASTQ line must begin with +."
            )

        if len(sequence) != len(quality):
            raise ValueError(
                "Sequence and quality lengths do not match "
                "for " + header
            )

        read_id = header[1:].split()[0]

        reads.append((read_id, sequence, quality))

    return reads


def print_summary(results):
    """Summarise all analysed reads."""

    print("\n" + "=" * 45)
    print("OVERALL QUALITY CONTROL SUMMARY")
    print("=" * 45)

    total = len(results)
    valid_results = [r for r in results if r["valid"]]
    invalid_results = [r for r in results if not r["valid"]]

    passing = [
        r for r in valid_results
        if r["qc_status"] == "PASS"
    ]

    requiring_review = [
        r for r in valid_results
        if r["qc_status"].startswith("REVIEW")
    ]

    print("Total reads:", total)
    print("Valid reads:", len(valid_results))
    print("Invalid reads:", len(invalid_results))
    print("Reads passing QC:", len(passing))
    print("Reads requiring review:", len(requiring_review))

    if valid_results:
        print(
            "Average read length:",
            round(mean(r["length"] for r in valid_results), 2)
        )

        print(
            "Average GC content:",
            round(mean(r["gc_content"] for r in valid_results), 2),
            "%"
        )

        quality_values = [
            r["average_quality"]
            for r in valid_results
            if r["average_quality"] is not None
        ]

        if quality_values:
            print(
                "Average Phred quality:",
                round(mean(quality_values), 2)
            )
        else:
            print("Average Phred quality: Not available")

    print("\nNote:")
    print("QC thresholds are illustrative, not universal.")
    print("This program is a learning tool, not a replacement")
    print("for professional tools such as FastQC.")
    print("=" * 45)


def analyze_single_read():
    """Ask the user to analyse one read."""

    print("\nEnter a DNA sequencing read.")
    sequence = input("DNA sequence: ")

    quality = input(
        "FASTQ quality string (press Enter to skip): "
    )

    sequence = clean_sequence(sequence)

    if quality == "":
        result = analyze_read("User_read", sequence)
    else:
        result = analyze_read("User_read", sequence, quality)

    print_read_report(result)


def analyze_multiple_reads():
    """Analyse several illustrative FASTQ reads."""

    example_fastq = """@Read1
ATGCGCTA
+
IIIIIIII
@Read2
GGCC
+
IIII
@Read3
ATATATAT
+
55555555
@Read4
ATBXGCTA
+
IIIIIIII
"""

    try:
        reads = parse_fastq(example_fastq)
    except ValueError as error:
        print("FASTQ parsing error:", error)
        return

    results = []

    for read_id, sequence, quality in reads:
        result = analyze_read(read_id, sequence, quality)
        results.append(result)
        print_read_report(result)

    print_summary(results)


def main():
    """Display the main menu."""

    while True:
        print("\n" + "=" * 45)
        print("NGS DATA ANALYSIS")
        print("Author:", AUTHOR)
        print("=" * 45)

        print("1. Analyse a single sequencing read")
        print("2. Analyse multiple example FASTQ reads")
        print("3. Exit")

        choice = input("Choose an option (1-3): ").strip()

        if choice == "1":
            analyze_single_read()

        elif choice == "2":
            analyze_multiple_reads()

        elif choice == "3":
            print("Thank you for using NGS Data Analysis!")
            break

        else:
            print("Invalid choice. Please enter 1, 2, or 3.")


if __name__ == "__main__":
    main()