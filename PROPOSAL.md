# Project Proposal — faxtools (FASTA Toolkit)

**Student:** Farnaz Yousefi
**Course:** BIOL 7180 — Scripting for Biologists (Spring 2026)
**Project type:** Solo
**GitHub repo:** https://github.com/Farnazyousefi95/biol7180-faxtools

## Title

**faxtools** — a lightweight FASTA summary + filtering command-line tool

## Goal

Develop a small, reproducible command-line script that summarizes FASTA files and optionally filters sequences by minimum length, producing outputs suitable for quality control (QC), reporting, and downstream analysis.

## Motivation / biological relevance

FASTA is a standard format for biological sequences (e.g., amplicons, reference genes, and assembled contigs). Before downstream analyses, researchers commonly need quick QC such as number of sequences, sequence-length statistics, total bases, and GC content. faxtools will provide these summaries in a consistent, scriptable format.

## Inputs and outputs

**Inputs:** One or more FASTA file paths (supports multi-line sequences).

**Outputs:**

1. **TSV summary** (to stdout or saved using `--out`) with per-file metrics: number of sequences, total bases, min/mean/max length, and GC%.
2. **Optional filtered FASTA** (using `--min-len` and `--filtered-out`) containing only sequences with length ≥ the specified threshold.

## Planned functionality (MVP)

* Robust FASTA parsing (multi-line sequences; ignores blank lines; validates header placement)
* Per-file statistics (n_seqs, total_bp, min/mean/max length, GC%)
* Simple command-line interface with clear help/usage examples
* Unit tests verifying counts and filtering behavior using small example FASTA files
* Tutorial-style README demonstrating how to run the script and interpret outputs

## Testing and reproducibility plan

* Include small example FASTA files in a `data/` directory
* Provide automated tests runnable with a single command (e.g., `python3 -m unittest -v`)
* Document exact run commands and expected output format in the README to support reproducibility

## Timeline

* **Week 1:** Implement FASTA parsing + TSV summary output
* **Week 2:** Add filtering option + finalize command-line arguments/help text
* **Week 3:** Add tests + finalize tutorial documentation and examples
* **Final:** Polish documentation and prepare a workshop-style walkthrough demonstration

## Deliverables

A public GitHub repository containing code, tests, example input files, and documentation sufficient for classmates to run the tool and reproduce example outputs without live assistance.

