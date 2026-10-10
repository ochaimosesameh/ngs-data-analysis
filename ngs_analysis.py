"""
NGS Data Analysis
Author: Ochai Moses Ameh

A beginner-friendly tool for analysing NGS reads
and their sequencing quality.
"""

from statistics import mean

AUTHOR = "Ochai Moses Ameh"

# Illustrative thresholds for this learning project.
MIN_READ_LENGTH = 5
MIN_AVERAGE_QUALITY = 20
MIN_GC_CONTENT = 30
MAX_GC_CONTENT = 70


def clean_sequence(sequence):
    """Convert a DNA sequence to uppercase."""
    return "".join(sequence.split()).upper()


def calculate_gc_content(sequence):
    """Calculate the percentage of G and C bases."""
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def count_bases(sequence):
    """Count A, T, G, C and N bases."""
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "N": sequence.count("N")
    }


def validate_sequence(sequence):
    """Find characters that are not accepted DNA bases."""
    valid_bases = set("ATGCN")
    return sorted(set(sequence) - valid_bases)


def phred_score_to_quality(character):
    """Convert a FASTQ Phred+33 character to a quality score."""
    return ord(character) - 33


def parse_fastq(filename):
    """Read sequencing reads from a FASTQ file."""
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.rstrip("\r\n") for line in file]

    if len(lines) % 4 != 0:
        raise ValueError(
            "FASTQ file must contain four lines per read."
        )

    reads = []

    for index in range(0, len(lines), 4):
        header = lines[index]
        sequence = lines[index + 1].strip()
        separator = lines[index + 2]
        quality = lines[index + 3]

        if not header.startswith("@"):
            raise ValueError(
                "A FASTQ header must begin with @."
            )

        if not separator.startswith("+"):
            raise ValueError(
                "The third FASTQ line must begin with +."
            )

        if len(sequence) != len(quality):
            raise ValueError(
                "Sequence and quality lengths differ for "
                + header
            )

        read_id = header[1:].split()[0]

        reads.append((read_id, sequence, quality))

    if not reads:
        raise ValueError("No sequencing reads found.")

    return reads


def analyze_read(read_id, sequence, quality_string):
    """Analyse one read and return its statistics."""
    sequence = clean_sequence(sequence)
    invalid_bases = validate_sequence(sequence)

    result = {
        "read_id": read_id,
        "sequence": sequence,
        "length": len(sequence),
        "base_counts": count_bases(sequence),
        "gc_content": calculate_gc_content(sequence),
        "average_quality": None,
        "valid": False,
        "qc_status": "FAIL",
        "warnings": []
    }

    if not sequence:
        result["warnings"].append("Empty sequence")
        return result

    if invalid_bases:
        result["warnings"].append(
            "Invalid characters: " + ", ".join(invalid_bases)
        )
        return result

    if len(quality_string) != len(sequence):
        result["warnings"].append(
            "Sequence and quality lengths differ"
        )
        return result

    quality_scores = [
        phred_score_to_quality(character)
        for character in quality_string
    ]

    if any(score < 0 or score > 93 for score in quality_scores):
        result["warnings"].append(
            "Invalid Phred+33 quality score"
        )
        return result

    result["valid"] = True
    result["average_quality"] = mean(quality_scores)

    if result["length"] < MIN_READ_LENGTH:
        result["warnings"].append("Short read")

    if not (
        MIN_GC_CONTENT
        <= result["gc_content"]
        <= MAX_GC_CONTENT
    ):
        result["warnings"].append("Review GC content")

    if result["average_quality"] < MIN_AVERAGE_QUALITY:
        result["warnings"].append("Low average quality")

    if result["warnings"]:
        result["qc_status"] = "REVIEW"
    else:
        result["qc_status"] = "PASS"

    return result


def print_read_report(result):
    """Print the analysis report for one read."""
    print("\n" + "=" * 40)
    print("READ ID:", result["read_id"])
    print("Sequence:", result["sequence"])
    print("Length:", result["length"])

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

    print("Valid sequence:", result["valid"])
    print("QC status:", result["qc_status"])

    if result["warnings"]:
        print("Warnings:", "; ".join(result["warnings"]))

    print("=" * 40)


def print_summary(results):
    """Summarise the results for all reads."""
    valid_reads = [r for r in results if r["valid"]]
    invalid_reads = [r for r in results if not r["valid"]]

    passing = [
        r for r in valid_reads if r["qc_status"] == "PASS"
    ]

    review = [
        r for r in valid_reads if r["qc_status"] == "REVIEW"
    ]

    print("\n" + "=" * 40)
    print("OVERALL QC SUMMARY")
    print("=" * 40)
    print("Total reads:", len(results))
    print("Valid reads:", len(valid_reads))
    print("Invalid reads:", len(invalid_reads))
    print("Reads passing QC:", len(passing))
    print("Reads requiring review:", len(review))

    if valid_reads:
        print(
            "Average read length:",
            round(mean(r["length"] for r in valid_reads), 2)
        )

        print(
            "Average GC content:",
            round(mean(r["gc_content"] for r in valid_reads), 2),
            "%"
        )

        quality_values = [
            r["average_quality"]
            for r in valid_reads
            if r["average_quality"] is not None
        ]

        if quality_values:
            print(
                "Average Phred quality:",
                round(mean(quality_values), 2)
            )

    print("\nQC thresholds are illustrative, not universal.")
    print("=" * 40)


def main():
    """Run the NGS analysis program."""
    print("\nNGS DATA ANALYSIS")
    print("Author:", AUTHOR)

    filename = input(
        "Enter FASTQ filename (default: sample_reads.fastq): "
    ).strip()

    if not filename:
        filename = "sample_reads.fastq"

    try:
        reads = parse_fastq(filename)
    except (OSError, ValueError) as error:
        print("Could not analyse FASTQ file:", error)
        return

    results = []

    for read_id, sequence, quality in reads:
        result = analyze_read(read_id, sequence, quality)
        results.append(result)
        print_read_report(result)

    print_summary(results)


if __name__ == "__main__":
    main()