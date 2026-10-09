"""
NGS Data Analysis
Author: Ochai Moses Ameh

A beginner-friendly program for basic analysis of
Next-Generation Sequencing (NGS) reads.
"""

from statistics import mean


AUTHOR = "Ochai Moses Ameh"


def calculate_gc_content(sequence):
    """Calculate the percentage of G and C bases."""
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def phred_quality_scores(quality_string):
    """Convert FASTQ quality characters to Phred+33 scores."""
    return [ord(character) - 33 for character in quality_string]


def analyze_read(read_id, sequence, quality_string):
    """Analyze one sequencing read."""
    sequence = "".join(sequence.upper().split())
    quality_string = quality_string.strip()

    result = {
        "read_id": read_id,
        "sequence": sequence,
        "length": len(sequence),
        "status": "Valid",
        "error": None,
        "a_count": sequence.count("A"),
        "t_count": sequence.count("T"),
        "g_count": sequence.count("G"),
        "c_count": sequence.count("C"),
        "n_count": sequence.count("N"),
        "gc_content": 0.0,
        "average_quality": None,
        "qc_status": "Review",
    }

    if not sequence:
        result["status"] = "Invalid"
        result["error"] = "Empty sequence"
        return result

    invalid_bases = sorted(set(sequence) - set("ATGCN"))

    if invalid_bases:
        result["status"] = "Invalid"
        result["error"] = (
            "Invalid bases: " + ", ".join(invalid_bases)
        )
        return result

    if len(sequence) != len(quality_string):
        result["status"] = "Invalid"
        result["error"] = (
            "Sequence and quality lengths do not match"
        )
        return result

    scores = phred_quality_scores(quality_string)

    # FASTQ quality characters should represent valid Phred+33 scores.
    if any(score < 0 or score > 93 for score in scores):
        result["status"] = "Invalid"
        result["error"] = "Invalid FASTQ quality character"
        return result

    result["gc_content"] = calculate_gc_content(sequence)
    result["average_quality"] = mean(scores)

    # These are simple educational thresholds, not universal standards.
    if result["average_quality"] >= 30:
        quality_ok = True
    else:
        quality_ok = False

    if 40 <= result["gc_content"] <= 60:
        gc_ok = True
    else:
        gc_ok = False

    if quality_ok and gc_ok:
        result["qc_status"] = "Pass"
    elif not quality_ok:
        result["qc_status"] = "Review quality"
    else:
        result["qc_status"] = "Review GC content"

    return result


def read_fastq_text(fastq_text):
    """Parse FASTQ-formatted text into sequencing reads."""
    lines = [
        line.strip()
        for line in fastq_text.strip().splitlines()
        if line.strip()
    ]

    if not lines:
        raise ValueError("The FASTQ data is empty.")

    if len(lines) % 4 != 0:
        raise ValueError(
            "FASTQ data must contain four lines per read."
        )

    reads = []

    for index in range(0, len(lines), 4):
        header = lines[index]
        sequence = lines[index + 1]
        plus_line = lines[index + 2]
        quality = lines[index + 3]

        if not header.startswith("@"):
            raise ValueError(
                f"Read header on line {index + 1} must start with @."
            )

        if not plus_line.startswith("+"):
            raise ValueError(
                f"Plus line for {header} must start with +."
            )

        read_id = header[1:].split()[0]

        reads.append(
            (read_id, sequence, quality)
        )

    return reads


def print_read_report(result):
    """Print the analysis results for one read."""
    print("\n" + "=" * 45)
    print("NGS READ ANALYSIS")
    print("=" * 45)

    print("Read ID:", result["read_id"])
    print("Status:", result["status"])

    if result["error"]:
        print("Error:", result["error"])
        return

    print("Sequence:", result["sequence"])
    print("Read length:", result["length"])
    print("A count:", result["a_count"])
    print("T count:", result["t_count"])
    print("G count:", result["g_count"])
    print("C count:", result["c_count"])
    print("N count:", result["n_count"])
    print("GC content:", round(result["gc_content"], 2), "%")

    print(
        "Average Phred quality:",
        round(result["average_quality"], 2),
    )

    print("QC status:", result["qc_status"])


def print_summary(results):
    """Print an overall summary of sequencing reads."""
    total = len(results)

    valid_results = [
        result for result in results
        if result["status"] == "Valid"
    ]

    invalid_results = [
        result for result in results
        if result["status"] == "Invalid"
    ]

    passing_results = [
        result for result in valid_results
        if result["qc_status"] == "Pass"
    ]

    low_quality_results = [
        result for result in valid_results
        if result["average_quality"] < 30
    ]

    gc_review_results = [
        result for result in valid_results
        if result["gc_content"] < 40
        or result["gc_content"] > 60
    ]

    print("\n" + "=" * 45)
    print("OVERALL QUALITY CONTROL SUMMARY")
    print("=" * 45)

    print("Total reads:", total)
    print("Valid reads:", len(valid_results))
    print("Invalid reads:", len(invalid_results))
    print("Reads passing educational QC:", len(passing_results))
    print("Reads with average quality below 30:",
          len(low_quality_results))
    print("Reads requiring GC review:", len(gc_review_results))

    if valid_results:
        print(
            "Average read length:",
            round(mean(r["length"] for r in valid_results), 2),
        )

        print(
            "Average GC content:",
            round(mean(r["gc_content"] for r in valid_results), 2),
            "%",
        )

        print(
            "Average Phred quality:",
            round(
                mean(r["average_quality"] for r in valid_results),
                2,
            ),
        )

    print("\nNote:")
    print("The dataset is illustrative.")
    print("QC thresholds are educational, not universal.")
    print("Use established tools such as FastQC for real datasets.")


def get_sample_fastq():
    """Return a small built-in example FASTQ dataset."""
    return """@Read1
ATGCGCTA
+
IIIIIIII
@Read2
GGCCTTAA
+
55555555
@Read3
ATATNNGC
+
IIIIIIII
@Read4
GGCCGGCC
+
IIIIIIII
"""


def analyze_fastq(fastq_text):
    """Analyze every read in FASTQ-formatted text."""
    try:
        reads = read_fastq_text(fastq_text)
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
    """Run the NGS Data Analysis program."""
    print("=" * 45)
    print("NGS DATA ANALYSIS")
    print("Author:", AUTHOR)
    print("=" * 45)

    while True:
        print("\n1. Analyze built-in sample FASTQ data")
        print("2. Enter one sequencing read manually")
        print("3. Exit")

        choice = input("Choose an option (1-3): ").strip()

        if choice == "1":
            analyze_fastq(get_sample_fastq())

        elif choice == "2":
            read_id = input("Enter read ID: ").strip()
            sequence = input("Enter DNA sequence: ").strip()
            quality = input(
                "Enter FASTQ quality string of equal length: "
            ).strip()

            if not read_id:
                read_id = "User_read"

            result = analyze_read(read_id, sequence, quality)
            print_read_report(result)

        elif choice == "3":
            print("Thank you for using NGS Data Analysis!")
            break

        else:
            print("Invalid choice. Please enter 1, 2, or 3.")


if __name__ == "__main__":
    main()