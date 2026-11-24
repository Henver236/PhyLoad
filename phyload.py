
#!/usr/bin/env python3

import argparse
import re
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------



REQUEST_DELAY = 0.4

MARKERS = {
    "ITS": ["ITS"],
    "GAPDH": ["GAPDH"],
    "CHS-1": ["CHS-1", "CHS1"],
    "Histone-3": ["Histone 3", "Histone-3", "Histone3"],
    "Actin": ["Actin"],
    "Tubulin-2": ["Tubulin 2", "Tubulin-2", "TUB2"],
    "Calmodulin": ["Calmodulin"],
}

MISSING_VALUES = {"NO", "NA", "N/A", ""}


# ---------------------------------------------------------------------
# Input table
# ---------------------------------------------------------------------

def load_input_table(input_file):
    """Load and normalize the input TSV file."""

    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    df.columns = [str(col).strip() for col in df.columns]

    # Normalize marker column names
    normalized_columns = {}

    for marker, aliases in MARKERS.items():
        for alias in aliases:
            for col in df.columns:
                if col.strip().lower() == alias.lower():
                    normalized_columns[marker] = col
                    break

            if marker in normalized_columns:
                break

    return df, normalized_columns


# ---------------------------------------------------------------------
# GenBank download
# ---------------------------------------------------------------------

def download_genbank_record(accession, cache):
    """Download a GenBank record and cache successful results."""

    if accession in cache:
        return cache[accession]

    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )

        record = SeqIO.read(handle, "genbank")
        handle.close()

        if not record.seq:
            cache[accession] = None
            return None

        cache[accession] = record

        time.sleep(REQUEST_DELAY)

        return record

    except Exception as error:
        print(f"Download error for {accession}: {error}")

        cache[accession] = None

        time.sleep(REQUEST_DELAY)

        return None


# ---------------------------------------------------------------------
# Organism name extraction
# ---------------------------------------------------------------------

def get_genbank_organism(record):
    """
    Extract the organism name from GenBank.

    Priority:
    1. Organism
    2. Source
    3. Definition
    """

    organism = ""

    # Organism is normally stored in record.annotations
    organism = record.annotations.get("organism", "")

    if organism:
        return organism.strip()

    # Try the SOURCE and ORGANISM information in the GenBank record
    for feature in record.features:

        if feature.type == "source":

            qualifiers = feature.qualifiers

            if "organism" in qualifiers:
                return qualifiers["organism"][0].strip()

            if "source" in qualifiers:
                return qualifiers["source"][0].strip()

    # Last resort: use the GenBank definition
    definition = record.description.strip()

    if definition:
        return definition

    return "Unknown_organism"


def clean_species_name(name):
    """Convert organism names into FASTA-compatible identifiers."""

    name = str(name).strip()

    # Remove dots
    name = name.replace(".", "")

    # Replace spaces, hyphens and other punctuation with underscores
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)

    # Remove repeated and leading/trailing underscores
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")

    return name or "Unknown_organism"


# ---------------------------------------------------------------------
# FASTA output
# ---------------------------------------------------------------------

def write_fasta_record(output_file, accession, organism, marker, sequence):
    """Append one sequence to a FASTA file."""

    organism_id = clean_species_name(organism)

    header = f">{accession}_{organism_id}_{marker}"

    with open(output_file, "a", encoding="utf-8") as handle:

        handle.write(header + "\n")

        sequence = str(sequence)

        for i in range(0, len(sequence), 80):
            handle.write(sequence[i:i + 80] + "\n")


# ---------------------------------------------------------------------
# Marker processing
# ---------------------------------------------------------------------


def process_marker(
    df,
    marker,
    accession_column,
    output_file,
    cache
):
    """Download all accessions for one marker with terminal progress."""

    results = []

    stats = {
        "requested": 0,
        "downloaded": 0,
        "not_found": 0,
        "ignored": 0,
    }

    print(f"\nProcessing marker: {marker}", flush=True)

    for index, row in df.iterrows():

        accession = str(row[accession_column]).strip()

        if accession.upper() in MISSING_VALUES:

            stats["ignored"] += 1

            print(
                f"  Skipping {accession or 'empty value'}... "
                f"--> IGNORED",
                flush=True
            )

            results.append({
                "row": index + 1,
                "marker": marker,
                "accession": accession,
                "organism": "",
                "status": "ignored",
                "fasta_file": ""
            })

            continue

        stats["requested"] += 1

        print(
            f"  Downloading {accession}...",
            end=" ",
            flush=True
        )

        record = download_genbank_record(accession, cache)

        if record is None:

            stats["not_found"] += 1

            print("--> NOT FOUND !", flush=True)

            results.append({
                "row": index + 1,
                "marker": marker,
                "accession": accession,
                "organism": "",
                "status": "not_found",
                "fasta_file": ""
            })

            continue

        organism = get_genbank_organism(record)

        write_fasta_record(
            output_file,
            accession,
            organism,
            marker,
            record.seq
        )

        stats["downloaded"] += 1

        print("--> ok !", flush=True)

        results.append({
            "row": index + 1,
            "marker": marker,
            "accession": accession,
            "organism": organism,
            "status": "downloaded",
            "fasta_file": output_file.name
        })

    return stats, results

# ---------------------------------------------------------------------
# Log output
# ---------------------------------------------------------------------

def write_log(log_file, marker_stats, all_results):

    """Write a detailed text log."""

    with open(log_file, "w", encoding="utf-8") as handle:

        handle.write("GENBANK REFERENCE SEQUENCE DOWNLOAD LOG\n")
        handle.write("=" * 55 + "\n\n")

        total_requested = 0
        total_downloaded = 0
        total_not_found = 0
        total_ignored = 0

        for marker, stats in marker_stats.items():

            handle.write(f"MARKER: {marker}\n")
            handle.write("-" * 40 + "\n")

            handle.write(
                f"Requested:  {stats['requested']}\n"
            )

            handle.write(
                f"Downloaded: {stats['downloaded']}\n"
            )

            handle.write(
                f"Not found:  {stats['not_found']}\n"
            )

            handle.write(
                f"Ignored:    {stats['ignored']}\n\n"
            )

            total_requested += stats["requested"]
            total_downloaded += stats["downloaded"]
            total_not_found += stats["not_found"]
            total_ignored += stats["ignored"]

        handle.write("\nGLOBAL SUMMARY\n")
        handle.write("=" * 55 + "\n")

        handle.write(f"Requested:  {total_requested}\n")
        handle.write(f"Downloaded: {total_downloaded}\n")
        handle.write(f"Not found:  {total_not_found}\n")
        handle.write(f"Ignored:    {total_ignored}\n")

        handle.write("\nACCESSIONS NOT FOUND\n")
        handle.write("=" * 55 + "\n")

        missing = [
            result for result in all_results
            if result["status"] == "not_found"
        ]

        if not missing:
            handle.write("No missing accessions.\n")

        else:
            for result in missing:

                handle.write(
                    f"{result['marker']}: "
                    f"{result['accession']} "
                    f"(input row {result['row']})\n"
                )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():

    parser = argparse.ArgumentParser(
        description="Download GenBank reference sequences from accession numbers."
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Email address required by NCBI Entrez."
    )

    parser.add_argument(
        "--input",
        default=str(INPUT_FILE),
        help="Input TSV file."
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to the input file directory."
    )

    args = parser.parse_args()

    Entrez.email = args.email

    input_file = Path(args.input)

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )

    # Create timestamped output directory
    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")

    if args.outdir:
        base_output_dir = Path(args.outdir)
    else:
        base_output_dir = input_file.parent

    output_dir = base_output_dir / f"ref_locus_{timestamp}"

    output_dir.mkdir(parents=True, exist_ok=True)

    output_dir = output_dir.resolve()


    base_name = input_file.stem

    df, marker_columns = load_input_table(input_file)

    cache = {}

    marker_stats = {}
    all_results = []

    print("\n=======================================================")
    print("GENBANK REFERENCE SEQUENCE DOWNLOAD")
    print("=======================================================")
    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")

    for marker in MARKERS:

        if marker not in marker_columns:
            print(
                f"\nWarning: marker column not found: {marker}",
                flush=True
            )
            continue

        accession_column = marker_columns[marker]

        fasta_file = output_dir / (
            f"{base_name}_{marker}_{timestamp}.fas"
        )

        stats, results = process_marker(
            df,
            marker,
            accession_column,
            fasta_file,
            cache
        )

        marker_stats[marker] = stats
        all_results.extend(results)

    # Write detailed CSV summary
    summary_file = output_dir / (
        f"{base_name}_summary_{timestamp}.csv"
    )

    summary_df = pd.DataFrame(all_results)

    summary_df.to_csv(
        summary_file,
        index=False,
        sep=";",
        encoding="utf-8-sig"
    )

    # Write text log
    log_file = output_dir / (
        f"{base_name}_download_log_{timestamp}.txt"
    )

    write_log(
        log_file,
        marker_stats,
        all_results
    )

    # Display final summary in terminal
    print("\n")
    print("=======================================================")
    print("DOWNLOAD COMPLETED")
    print("=======================================================")

    print(
        f"{'Marker':<15}"
        f"{'Downloaded':>12} | "
        f"{'Not found':>10} | "
        f"{'Ignored':>8}"
    )

    for marker, stats in marker_stats.items():

        print(
            f"{marker:<15}"
            f"{stats['downloaded']:>12} | "
            f"{stats['not_found']:>10} | "
            f"{stats['ignored']:>8}"
        )

    print("\n-------------------------------------------------------")
    print(f"Log file: {log_file}")
    print(f"Summary CSV: {summary_file}")
    print(f"Output directory: {output_dir}")
    print("=======================================================")


if __name__ == "__main__":
    main()
