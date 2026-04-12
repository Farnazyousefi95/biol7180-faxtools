#!/usr/bin/env python3
"""
faxtools — a lightweight FASTA summary and filtering command-line tool.

This script reads one or more FASTA files and reports per-file summary
statistics (number of sequences, total bases, min/mean/max length, and GC
content).  It can optionally filter sequences by a minimum length threshold
and write the passing sequences to a new FASTA file.

Usage examples
--------------
Summarize a single file to stdout::

    python3 faxtools.py data/example-seqs1.fasta

Summarize several files and save a TSV::

    python3 faxtools.py data/example-seqs1.fasta data/example-seqs2.fasta \
        --out summary.tsv

Filter by minimum length and write filtered FASTA::

    python3 faxtools.py data/example-seqs1.fasta --min-len 90 \
        --filtered-out filtered.fasta

Author
------
Farnaz Yousefi — BIOL 7180, Spring 2026
"""

import argparse
import os
import sys
from statistics import mean


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(filepath):
    """Parse a FASTA file and yield (header, sequence) tuples.

    Supports multi-line sequences, ignores blank lines, and raises
    ``ValueError`` if the file does not begin with a '>' header line (after
    skipping any leading blank lines).

    Parameters
    ----------
    filepath : str
        Path to a FASTA-formatted file.

    Yields
    ------
    tuple of (str, str)
        A two-element tuple where the first element is the header line
        (without the leading '>') and the second is the concatenated
        sequence string.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist.
    ValueError
        If the first non-blank line is not a FASTA header ('>').
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    header = None
    sequence_parts = []

    with open(filepath, "r") as fasta_file:
        for line_number, raw_line in enumerate(fasta_file, start=1):
            line = raw_line.strip()

            # Skip blank lines
            if not line:
                continue

            if line.startswith(">"):
                # Yield the previous record (if any)
                if header is not None:
                    yield (header, "".join(sequence_parts))
                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"{filepath}, line {line_number}: expected a '>' "
                        f"header before sequence data"
                    )
                sequence_parts.append(line)

        # Yield the last record
        if header is not None:
            yield (header, "".join(sequence_parts))


# ---------------------------------------------------------------------------
# Sequence-level helpers
# ---------------------------------------------------------------------------

def calculate_gc_percent(sequence):
    """Return the GC percentage of *sequence*, ignoring gaps ('-').

    The calculation is case-insensitive.  If the sequence (after removing
    gaps) is empty, 0.0 is returned.

    Parameters
    ----------
    sequence : str
        A nucleotide sequence string.

    Returns
    -------
    float
        GC content as a percentage (0–100).
    """
    bases = sequence.upper().replace("-", "")
    if not bases:
        return 0.0
    gc_count = bases.count("G") + bases.count("C")
    return (gc_count / len(bases)) * 100.0


def sequence_length(sequence):
    """Return the number of non-gap characters in *sequence*.

    Parameters
    ----------
    sequence : str
        A nucleotide or protein sequence string.

    Returns
    -------
    int
        Length of the sequence after removing gap characters ('-').
    """
    return len(sequence.replace("-", ""))


# ---------------------------------------------------------------------------
# Per-file summary
# ---------------------------------------------------------------------------

def summarize_fasta(filepath):
    """Compute summary statistics for all sequences in a FASTA file.

    Parameters
    ----------
    filepath : str
        Path to a FASTA-formatted file.

    Returns
    -------
    dict
        A dictionary with the following keys:

        - ``file``        – the basename of the input file
        - ``n_seqs``      – number of sequences
        - ``total_bp``    – total bases (excluding gaps)
        - ``min_len``     – shortest sequence length
        - ``mean_len``    – mean sequence length (rounded to 1 decimal)
        - ``max_len``     – longest sequence length
        - ``gc_percent``  – overall GC% (rounded to 2 decimals)

    Raises
    ------
    ValueError
        If the file contains zero sequences.
    """
    lengths = []
    gc_numerator = 0   # total G + C bases across all sequences
    total_bases = 0

    for _header, seq in parse_fasta(filepath):
        bases = seq.upper().replace("-", "")
        seq_len = len(bases)
        lengths.append(seq_len)
        total_bases += seq_len
        gc_numerator += bases.count("G") + bases.count("C")

    if not lengths:
        raise ValueError(f"No sequences found in {filepath}")

    overall_gc = (gc_numerator / total_bases) * 100.0 if total_bases else 0.0

    return {
        "file": os.path.basename(filepath),
        "n_seqs": len(lengths),
        "total_bp": total_bases,
        "min_len": min(lengths),
        "mean_len": round(mean(lengths), 1),
        "max_len": max(lengths),
        "gc_percent": round(overall_gc, 2),
    }


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_sequences(filepaths, min_len):
    """Yield (header, sequence) tuples that meet a minimum-length threshold.

    Parameters
    ----------
    filepaths : list of str
        Paths to FASTA-formatted files.
    min_len : int
        Minimum sequence length (non-gap characters).  Sequences shorter
        than this value are excluded.

    Yields
    ------
    tuple of (str, str)
        (header, sequence) for every sequence whose non-gap length is
        >= *min_len*.
    """
    for filepath in filepaths:
        for header, seq in parse_fasta(filepath):
            if sequence_length(seq) >= min_len:
                yield (header, seq)


def write_fasta(records, output_path, line_width=80):
    """Write (header, sequence) records to a FASTA file.

    Parameters
    ----------
    records : iterable of (str, str)
        Each element is a (header, sequence) tuple.
    output_path : str
        Destination file path.
    line_width : int, optional
        Maximum number of characters per sequence line (default 80).
    """
    with open(output_path, "w") as out_fh:
        for header, seq in records:
            out_fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                out_fh.write(seq[i:i + line_width] + "\n")


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------

TSV_COLUMNS = ["file", "n_seqs", "total_bp", "min_len", "mean_len",
               "max_len", "gc_percent"]


def format_summary_tsv(summaries):
    """Format a list of summary dicts as a TSV string (with header row).

    Parameters
    ----------
    summaries : list of dict
        Each dict should contain the keys listed in ``TSV_COLUMNS``.

    Returns
    -------
    str
        Tab-separated text with a header row and one data row per file.
    """
    lines = ["\t".join(TSV_COLUMNS)]
    for summary in summaries:
        row = "\t".join(str(summary[col]) for col in TSV_COLUMNS)
        lines.append(row)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def build_argument_parser():
    """Create and return the ``argparse.ArgumentParser`` for faxtools.

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarize one or more FASTA files (sequence counts, length "
            "statistics, GC%%) and optionally filter sequences by minimum "
            "length."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python3 faxtools.py data/example-seqs1.fasta\n"
            "  python3 faxtools.py *.fasta --out summary.tsv\n"
            "  python3 faxtools.py *.fasta --min-len 100 "
            "--filtered-out long_seqs.fasta\n"
        ),
    )
    parser.add_argument(
        "fasta_files",
        nargs="+",
        metavar="FASTA",
        help="One or more FASTA files to process.",
    )
    parser.add_argument(
        "--out",
        metavar="TSV",
        default=None,
        help="Write the summary table to this TSV file (default: stdout).",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=None,
        metavar="N",
        help="Minimum sequence length for filtering (requires --filtered-out).",
    )
    parser.add_argument(
        "--filtered-out",
        metavar="FASTA",
        default=None,
        help="Write sequences passing --min-len to this FASTA file.",
    )
    return parser


def main(argv=None):
    """Entry point for the faxtools command-line tool.

    Parameters
    ----------
    argv : list of str or None
        Command-line arguments (defaults to ``sys.argv[1:]``).
    """
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    # --- Validate argument combinations --------------------------------
    if args.min_len is not None and args.filtered_out is None:
        parser.error("--min-len requires --filtered-out")
    if args.filtered_out is not None and args.min_len is None:
        parser.error("--filtered-out requires --min-len")

    # --- Summarize each FASTA file -------------------------------------
    summaries = []
    for filepath in args.fasta_files:
        try:
            summary = summarize_fasta(filepath)
            summaries.append(summary)
        except (FileNotFoundError, ValueError) as err:
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)

    tsv_text = format_summary_tsv(summaries)

    if args.out:
        with open(args.out, "w") as tsv_fh:
            tsv_fh.write(tsv_text)
        print(f"Summary written to {args.out}", file=sys.stderr)
    else:
        print(tsv_text, end="")

    # --- Optional filtering --------------------------------------------
    if args.min_len is not None:
        filtered_records = list(
            filter_sequences(args.fasta_files, args.min_len)
        )
        write_fasta(filtered_records, args.filtered_out)
        print(
            f"Wrote {len(filtered_records)} sequence(s) with length >= "
            f"{args.min_len} to {args.filtered_out}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
